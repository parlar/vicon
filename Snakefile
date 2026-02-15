"""
ViCon — Viral Consensus Pipeline (Snakemake)

Iterative reference mapping → de novo assembly (Trinity, SPAdes, coronaSPAdes,
MEGAHIT, IVA) → iVar consensus → assembly-guided polishing → phylogenetic
analysis (IQ-TREE) with homoplasy metrics. Runs dedup and nodedup lines
in parallel for each sample.
"""

configfile: "config.yaml"

from pathlib import Path


# ── Shared utilities ──────────────────────────────────────────

def read_fasta(path):
    """Parse a FASTA file → dict of {name: sequence (uppercased)}."""
    seqs = {}
    name = None
    parts = []
    with open(str(path)) as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith(">"):
                if name:
                    seqs[name] = "".join(parts).upper()
                name = ln[1:].split()[0]
                parts = []
            else:
                parts.append(ln)
    if name:
        seqs[name] = "".join(parts).upper()
    return seqs


from collections import Counter
from io import StringIO

try:
    from Bio import Phylo
except ImportError:
    Phylo = None  # deferred; only needed at rule execution time


# ── Fitch parsimony utilities ─────────────────────────────────

def fitch_bottom_up(clade, leaf_states, node_cache):
    """Fitch bottom-up pass for one site. Returns total changes."""
    if clade.is_terminal():
        name = clade.name
        if name in leaf_states:
            node_cache[id(clade)] = {leaf_states[name]}
        else:
            node_cache[id(clade)] = set("ACGT")  # missing → any state
        return 0

    changes = 0
    child_sets = []
    for child in clade.clades:
        changes += fitch_bottom_up(child, leaf_states, node_cache)
        child_sets.append(node_cache[id(child)])

    result = child_sets[0].copy()
    for cs in child_sets[1:]:
        inter = result & cs
        if inter:
            result = inter
        else:
            result = result | cs
            changes += 1

    node_cache[id(clade)] = result
    return changes


def fitch_top_down(clade, node_cache, assigned, parent_state=None):
    """Assign a single state per node, top-down."""
    state_set = node_cache.get(id(clade), set())
    if parent_state is not None and parent_state in state_set:
        assigned[id(clade)] = parent_state
    elif state_set:
        assigned[id(clade)] = min(state_set)  # deterministic pick
    else:
        assigned[id(clade)] = parent_state  # fallback
    for child in clade.clades:
        fitch_top_down(child, node_cache, assigned, assigned[id(clade)])


def build_parent_map(tree):
    """Build dict mapping node id → parent clade."""
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parents[id(child)] = clade
    return parents


# ── Configuration ──────────────────────────────────────────────

FASTQ_DIR  = config["fastq_dir"]
REF_DIR    = config["ref_dir"]
OUTDIR     = config["outdir"]
SAMPLE_DIR = OUTDIR + "/samples"
COMMON_DIR = OUTDIR + "/common"
ALIGNER    = config.get("aligner", "bwa-mem2")
MAX_ITERATIONS = config.get("max_iterations", 10)
MIN_ITERATIONS = config.get("min_iterations", 2)
THREADS    = config.get("threads", 8)
LINES      = config.get("lines", ["dedup", "nodedup"])
ENFORCE_INFRAME = config.get("enforce_inframe", True)
TRINITY_MEMORY  = config.get("trinity_memory", 4)
CHECK_STOPS     = config.get("check_stops", True)
ASSEMBLERS      = config.get("assemblers", ["trinity", "spades", "coronaspades", "megahit", "iva"])
SPADES_MEMORY   = config.get("spades_memory", 8)
MEGAHIT_MEMORY  = config.get("megahit_memory", 8)
HUMAN_BT2_INDEX    = config["human_bt2_index"]
BANKVOLE_BT2_INDEX = config["bankvole_bt2_index"]
GENBANK_DIR        = config.get("genbank_dir", "../genbank_reference")
SEGMENTS           = config.get("segments", ["L", "M", "S"])
_default_refs      = {s: f"{GENBANK_DIR}/{s}.full_len.cds.clean.uniq.fa" for s in SEGMENTS}
GENBANK_REFS       = {s: config.get("genbank_refs", {}).get(s, _default_refs[s]) for s in SEGMENTS}
AUTO_SELECT_REFS   = config.get("auto_select_refs", False)
DENOVO_SEED        = config.get("denovo_seed", False)
EXTEND_REFERENCE   = config.get("extend_reference", False)
MAX_GENBANK_REFS   = config.get("max_genbank_refs", 50)


# ── Sample detection ───────────────────────────────────────────

def detect_samples():
    fastq_dir = Path(FASTQ_DIR)
    samples = {}
    for r1 in sorted(fastq_dir.glob("*_R1_*.fastq.gz")):
        name = r1.name.split("_")[0]
        r2 = fastq_dir / r1.name.replace("_R1_", "_R2_")
        if r2.exists():
            samples[name] = {"r1": str(r1), "r2": str(r2)}
    return samples

ALL_SAMPLES = detect_samples()

if "samples" in config and config["samples"]:
    SAMPLES = {s: ALL_SAMPLES[s] for s in config["samples"] if s in ALL_SAMPLES}
else:
    SAMPLES = ALL_SAMPLES

if not SAMPLES:
    raise ValueError(f"No paired-end samples found in {FASTQ_DIR}")


# ── Input helper functions ─────────────────────────────────────

def get_raw_r1(wildcards):
    return SAMPLES[wildcards.sample]["r1"]

def get_raw_r2(wildcards):
    return SAMPLES[wildcards.sample]["r2"]

def get_r1(wildcards):
    return f"{SAMPLE_DIR}/{wildcards.sample}/host_filtered/clean_R1.fastq.gz"

def get_r2(wildcards):
    return f"{SAMPLE_DIR}/{wildcards.sample}/host_filtered/clean_R2.fastq.gz"

def _analysis_bam(base, line):
    """Return dedup BAM path for dedup line, raw BAM for nodedup."""
    if line == "dedup":
        return f"{base}/mapped.dedup.bam"
    return f"{base}/mapped.bam"

def get_assembly_fasta(wildcards):
    """Return the assembly FASTA path for a given assembler."""
    paths = {
        "trinity":      f"{SAMPLE_DIR}/{wildcards.sample}/{wildcards.line}/trinity/assembly/Trinity.fasta",
        "spades":       f"{SAMPLE_DIR}/{wildcards.sample}/{wildcards.line}/spades/contigs.fasta",
        "coronaspades": f"{SAMPLE_DIR}/{wildcards.sample}/{wildcards.line}/coronaspades/contigs.fasta",
        "megahit":      f"{SAMPLE_DIR}/{wildcards.sample}/{wildcards.line}/megahit/final.contigs.fa",
        "iva":          f"{SAMPLE_DIR}/{wildcards.sample}/{wildcards.line}/iva/contigs.fasta",
    }
    return paths[wildcards.assembler]

def get_ref_input(wildcards):
    """Return reference FASTA: auto-selected or pre-existing."""
    if AUTO_SELECT_REFS:
        return f"{SAMPLE_DIR}/{wildcards.sample}/autoref/best_refs.fasta"
    return f"{REF_DIR}/{wildcards.sample}_best_refs.fasta"

def get_iterative_ref(wildcards):
    """Return starting reference for iterative mapping: seeded or original."""
    if DENOVO_SEED:
        return f"{SAMPLE_DIR}/{wildcards.sample}/denovo_seed/seeded_ref.fasta"
    return f"{SAMPLE_DIR}/{wildcards.sample}/reference/ref.fasta"

def get_iterative_ref_fai(wildcards):
    """Return .fai for the iterative mapping starting reference."""
    if DENOVO_SEED:
        return f"{SAMPLE_DIR}/{wildcards.sample}/denovo_seed/seeded_ref.fasta.fai"
    return f"{SAMPLE_DIR}/{wildcards.sample}/reference/ref.fasta.fai"

def get_genbank_ref(wildcards):
    """Return GenBank reference FASTA for a given segment."""
    return GENBANK_REFS[wildcards.segment]

def get_validated_genbank_ref(wildcards):
    """Return CDS-validated GenBank reference FASTA for a given segment."""
    return COMMON_DIR + f"/genbank_validated/{wildcards.segment}.validated.fa"


# ── Wildcard constraints ──────────────────────────────────────

wildcard_constraints:
    sample    = r"[^/]+",
    line      = r"dedup|nodedup",
    assembler = r"trinity|spades|coronaspades|megahit|iva",
    segment   = "|".join(SEGMENTS),


# ═══════════════════════════════════════════════════════════════
# TARGET
# ═══════════════════════════════════════════════════════════════

rule all:
    input:
        expand(
            "{sd}/{sample}/{line}/summary.tsv",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/homoplasy/homoplasy_report.tsv",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/iqtree/{segment}/iqtree.treefile",
            sd=SAMPLE_DIR,
            line=LINES,
            segment=SEGMENTS,
            sample=SAMPLES.keys(),
        ),
        COMMON_DIR + "/iqtree/iqtree_summary.tsv",
        COMMON_DIR + "/iqtree/individual_placement_homoplasy.tsv",
        expand(
            "{sd}/{sample}/{line}/best_consensus/best_consensus.fasta",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/best_consensus/validation.tsv",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/best_consensus/read_support.tsv",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/report.html",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/{line}/README.md",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{sd}/{sample}/README.md",
            sd=SAMPLE_DIR,
            sample=SAMPLES.keys(),
        ),
        COMMON_DIR + "/overview_report.html",


# ── Module includes (ordered by pipeline stage) ──────────────
include: "rules/preprocessing.smk"
include: "rules/reference.smk"
include: "rules/denovo_seed.smk"
include: "rules/iterative_mapping.smk"
include: "rules/assembly.smk"
include: "rules/consensus.smk"
include: "rules/alignment_and_homoplasy.smk"
include: "rules/phylogenetics.smk"
include: "rules/homoplasy_summary.smk"
include: "rules/best_consensus.smk"
include: "rules/reporting.smk"
