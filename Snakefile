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


# ═══════════════════════════════════════════════════════════════
# ADAPTER TRIMMING
# ═══════════════════════════════════════════════════════════════

rule adapter_trim:
    """Trim adapters and low-quality bases from raw reads using fastp."""
    input:
        r1 = get_raw_r1,
        r2 = get_raw_r2,
    output:
        r1   = SAMPLE_DIR + "/{sample}/trimmed/trimmed_R1.fastq.gz",
        r2   = SAMPLE_DIR + "/{sample}/trimmed/trimmed_R2.fastq.gz",
        html = SAMPLE_DIR + "/{sample}/trimmed/fastp.html",
        json = SAMPLE_DIR + "/{sample}/trimmed/fastp.json",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.r1})

        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --detect_adapter_for_pe \
            --html {output.html} --json {output.json} \
            --thread {threads}
        """


# ═══════════════════════════════════════════════════════════════
# HOST REMOVAL
# ═══════════════════════════════════════════════════════════════

rule host_removal:
    """Remove human and bank vole host reads using bowtie2 --very-sensitive."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/trimmed/trimmed_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/trimmed/trimmed_R2.fastq.gz",
    output:
        r1 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R2.fastq.gz",
    params:
        human_idx    = HUMAN_BT2_INDEX,
        bankvole_idx = BANKVOLE_BT2_INDEX,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.r1})
        mkdir -p "$DIR"

        # Step 1: Remove human reads
        bowtie2 --very-sensitive \
            -x {params.human_idx} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            --un-conc-gz "$DIR/tmp_no_human_R%.fastq.gz" \
            -S /dev/null 2> "$DIR/human_bowtie2.log"

        # Step 2: Remove bank vole reads from human-filtered output
        bowtie2 --very-sensitive \
            -x {params.bankvole_idx} \
            -1 "$DIR/tmp_no_human_R1.fastq.gz" \
            -2 "$DIR/tmp_no_human_R2.fastq.gz" \
            --threads {threads} \
            --un-conc-gz "$DIR/clean_R%.fastq.gz" \
            -S /dev/null 2> "$DIR/bankvole_bowtie2.log"

        # Clean up intermediate files
        rm -f "$DIR/tmp_no_human_R1.fastq.gz" "$DIR/tmp_no_human_R2.fastq.gz"
        """


# ═══════════════════════════════════════════════════════════════
# AUTOMATIC REFERENCE SELECTION (optional)
# ═══════════════════════════════════════════════════════════════

rule autoref_map:
    """Map host-filtered reads against all GenBank refs for one segment."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R2.fastq.gz",
        ref = get_genbank_ref,
    output:
        bam = SAMPLE_DIR + "/{sample}/autoref/{segment}/mapped.bam",
        bai = SAMPLE_DIR + "/{sample}/autoref/{segment}/mapped.bam.bai",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        minimap2 -a -x sr -t {threads} {input.ref} {input.r1} {input.r2} \
          | samtools sort -@ 4 -o {output.bam}
        samtools index {output.bam}
        """


rule autoref_select:
    """Select best reference per segment by completeness (fraction covered >2x)."""
    input:
        bam = SAMPLE_DIR + "/{sample}/autoref/{segment}/mapped.bam",
        bai = SAMPLE_DIR + "/{sample}/autoref/{segment}/mapped.bam.bai",
        ref = get_genbank_ref,
    output:
        stats = SAMPLE_DIR + "/{sample}/autoref/{segment}/ref_stats.tsv",
        best_fa = SAMPLE_DIR + "/{sample}/autoref/{segment}/best.fasta",
    run:
        import subprocess
        from collections import defaultdict

        # 1. Reference lengths + mapped read counts from idxstats
        result = subprocess.run(
            ["samtools", "idxstats", str(input.bam)],
            capture_output=True, text=True, check=True,
        )
        ref_info = {}
        for line in result.stdout.strip().split("\n"):
            parts = line.split("\t")
            if parts[0] == "*":
                continue
            ref_info[parts[0]] = {
                "length": int(parts[1]),
                "mapped_reads": int(parts[2]),
            }

        # 2. Completeness: fraction of positions with depth > 2
        #    Stream samtools depth to avoid loading all output into memory
        covered = defaultdict(int)
        proc = subprocess.Popen(
            ["samtools", "depth", "-a", str(input.bam)],
            stdout=subprocess.PIPE, text=True,
        )
        for line in proc.stdout:
            parts = line.split("\t")
            if int(parts[2]) > 2:
                covered[parts[0]] += 1
        proc.wait()

        # 3. Score each reference and write stats
        stats = []
        for name, info in ref_info.items():
            length = info["length"]
            mapped = info["mapped_reads"]
            cov_bases = covered.get(name, 0)
            completeness = cov_bases / length if length > 0 else 0.0
            stats.append({
                "reference": name,
                "length": length,
                "mapped_reads": mapped,
                "bases_covered_gt2": cov_bases,
                "completeness": round(completeness, 6),
            })

        # Primary sort: completeness desc, secondary: mapped_reads desc
        stats.sort(
            key=lambda x: (x["completeness"], x["mapped_reads"]),
            reverse=True,
        )

        import os
        os.makedirs(os.path.dirname(str(output.stats)), exist_ok=True)
        with open(str(output.stats), "w") as f:
            f.write("reference\tlength\tmapped_reads\t"
                    "bases_covered_gt2\tcompleteness\n")
            for s in stats:
                f.write(f"{s['reference']}\t{s['length']}\t"
                        f"{s['mapped_reads']}\t{s['bases_covered_gt2']}\t"
                        f"{s['completeness']}\n")

        # 4. Extract best reference sequence from the GenBank FASTA
        best_name = stats[0]["reference"] if stats else None
        seqs = {}
        headers = {}
        cur_name = None
        cur_seq = []
        with open(str(input.ref)) as f:
            for ln in f:
                ln = ln.rstrip("\n")
                if ln.startswith(">"):
                    if cur_name:
                        seqs[cur_name] = "".join(cur_seq)
                    cur_name = ln[1:].split()[0]
                    headers[cur_name] = ln
                    cur_seq = []
                else:
                    cur_seq.append(ln)
            if cur_name:
                seqs[cur_name] = "".join(cur_seq)

        with open(str(output.best_fa), "w") as f:
            if best_name and best_name in seqs:
                chosen = best_name
            else:
                # Fallback: longest sequence
                chosen = max(seqs, key=lambda k: len(seqs[k]))
            f.write(headers[chosen] + "\n")
            seq = seqs[chosen]
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


rule autoref_combine:
    """Combine best references for all segments into a single per-sample FASTA."""
    input:
        segments = expand(
            SAMPLE_DIR + "/{{sample}}/autoref/{segment}/best.fasta",
            segment=SEGMENTS,
        ),
    output:
        ref = SAMPLE_DIR + "/{sample}/autoref/best_refs.fasta",
    shell:
        "cat {input.segments} > {output.ref}"


# ═══════════════════════════════════════════════════════════════
# GENBANK REFERENCE SUBSAMPLING (Mash-based diversity selection)
# ═══════════════════════════════════════════════════════════════

rule subsample_genbank:
    """Subsample GenBank references to N diverse representatives using Mash distances."""
    input:
        fasta = get_genbank_ref,
    output:
        fasta = COMMON_DIR + "/genbank_subsampled/{segment}.subsampled.fa",
    params:
        max_refs = MAX_GENBANK_REFS,
    threads: 1
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})
        python3 {workflow.basedir}/scripts/mash_subsample.py \
            {input.fasta} {output.fasta} {params.max_refs}
        """


# ═══════════════════════════════════════════════════════════════
# REFERENCE PREPARATION
# ═══════════════════════════════════════════════════════════════

rule rename_ref:
    """Match reference contigs to configured segment names and rename to {segment}_{sample}."""
    input:
        ref = get_ref_input,
    output:
        ref = SAMPLE_DIR + "/{sample}/reference/ref.fasta",
    params:
        segments = SEGMENTS,
    run:
        import os

        seqs = read_fasta(str(input.ref))

        # Match each contig to a segment: contig name must start with the segment name
        matched = {}
        for seg in params.segments:
            for contig_name, seq in seqs.items():
                if contig_name == seg or contig_name.startswith(seg + "_"):
                    matched[seg] = seq
                    break

        missing = [s for s in params.segments if s not in matched]
        if missing:
            raise ValueError(
                f"Sample {wildcards.sample}: could not match contigs to segments "
                f"{missing}. Reference contig names: {list(seqs.keys())}. "
                f"Contig names must start with the segment name."
            )

        os.makedirs(os.path.dirname(str(output.ref)), exist_ok=True)
        with open(str(output.ref), "w") as out:
            for seg in params.segments:
                seq = matched[seg]
                out.write(f">{seg}_{wildcards.sample}\n")
                for j in range(0, len(seq), 80):
                    out.write(seq[j:j+80] + "\n")


rule index_ref:
    """Index the renamed per-sample reference (shared by both lines)."""
    input:
        ref = SAMPLE_DIR + "/{sample}/reference/ref.fasta",
    output:
        fai = SAMPLE_DIR + "/{sample}/reference/ref.fasta.fai",
    params:
        aligner = ALIGNER,
    shell:
        """
        set -euo pipefail
        if [ "{params.aligner}" = "bwa-mem2" ] && [ ! -f "{input.ref}.bwt.2bit.64" ]; then
            bwa-mem2 index {input.ref}
        fi
        samtools faidx {input.ref}
        """


# ═══════════════════════════════════════════════════════════════
# DE NOVO SEED (optional: assemble first, then seed iterative mapping)
# ═══════════════════════════════════════════════════════════════

rule denovo_seed_assembly:
    """Quick de novo assembly with MEGAHIT to seed iterative mapping."""
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        contigs = SAMPLE_DIR + "/{sample}/denovo_seed/megahit/final.contigs.fa",
    params:
        outdir = SAMPLE_DIR + "/{sample}/denovo_seed/megahit",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -t {threads} \
            -o {params.outdir}
        """


rule denovo_seed_scaffold:
    """Scaffold de novo contigs against reference to build a seeded reference."""
    input:
        contigs = SAMPLE_DIR + "/{sample}/denovo_seed/megahit/final.contigs.fa",
        ref     = SAMPLE_DIR + "/{sample}/reference/ref.fasta",
        fai     = SAMPLE_DIR + "/{sample}/reference/ref.fasta.fai",
    output:
        seeded = SAMPLE_DIR + "/{sample}/denovo_seed/seeded_ref.fasta",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.seeded})

        # Align de novo contigs to the reference
        minimap2 -ax asm5 -t {threads} {input.ref} {input.contigs} \
            | samtools sort -o "$DIR/contigs_vs_ref.bam"
        samtools index "$DIR/contigs_vs_ref.bam"

        # Call variants from the alignment
        bcftools mpileup -Ou -f {input.ref} "$DIR/contigs_vs_ref.bam" \
            | bcftools call -Ou -c --ploidy 1 \
            | bcftools norm -f {input.ref} -Oz -o "$DIR/seed_variants.vcf.gz"
        bcftools index "$DIR/seed_variants.vcf.gz"

        # Apply variants to the reference to produce seeded reference
        bcftools consensus -f {input.ref} "$DIR/seed_variants.vcf.gz" \
            > {output.seeded}
        """


rule index_seeded_ref:
    """Index the seeded reference."""
    input:
        ref = SAMPLE_DIR + "/{sample}/denovo_seed/seeded_ref.fasta",
    output:
        fai = SAMPLE_DIR + "/{sample}/denovo_seed/seeded_ref.fasta.fai",
    params:
        aligner = ALIGNER,
    shell:
        """
        set -euo pipefail
        if [ "{params.aligner}" = "bwa-mem2" ] && [ ! -f "{input.ref}.bwt.2bit.64" ]; then
            bwa-mem2 index {input.ref}
        fi
        samtools faidx {input.ref}
        """


# ═══════════════════════════════════════════════════════════════
# ITERATIVE MAPPING (convergence-based)
# ═══════════════════════════════════════════════════════════════

rule iterative_mapping:
    """Iteratively map reads until mapped reads and completeness stop improving."""
    input:
        ref = get_iterative_ref,
        fai = get_iterative_ref_fai,
        r1  = get_r1,
        r2  = get_r2,
    output:
        consensus = SAMPLE_DIR + "/{sample}/{line}/iterative/final_consensus.fasta",
        bam       = SAMPLE_DIR + "/{sample}/{line}/iterative/analysis.bam",
        bai       = SAMPLE_DIR + "/{sample}/{line}/iterative/analysis.bam.bai",
        stats     = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
    params:
        sample          = "{sample}",
        line            = "{line}",
        aligner         = ALIGNER,
        max_iterations  = MAX_ITERATIONS,
        min_iterations  = MIN_ITERATIONS,
        enforce_inframe = ENFORCE_INFRAME,
        extend_ref      = EXTEND_REFERENCE,
        outdir          = SAMPLE_DIR + "/{sample}/{line}/iterative",
    threads: THREADS
    shell:
        r"""
        set -euo pipefail

        OUTBASE="{params.outdir}"
        LINE="{params.line}"
        ALIGNER="{params.aligner}"
        MAX_ITER={params.max_iterations}
        MIN_ITER={params.min_iterations}
        ENFORCE="{params.enforce_inframe}"
        EXTEND_REF="{params.extend_ref}"
        SAMPLE="{params.sample}"

        REF="{input.ref}"
        PREV_MAPPED=0
        PREV_COMPLETENESS="0.00"
        ROUND=0
        LAST_CONSENSUS=""
        LAST_BAM=""

        mkdir -p "$OUTBASE"
        printf "round\tmapped_reads\tcompleteness\n" > {output.stats}

        while true; do
            ROUND=$((ROUND + 1))
            ROUND_DIR="$OUTBASE/round_$ROUND"
            mkdir -p "$ROUND_DIR"

            # Copy reference used this round
            cp "$REF" "$ROUND_DIR/reference.fasta"

            # ── Map reads (index + align + add RG + index BAM) ──
            bash {workflow.basedir}/scripts/map_reads.sh \
                "$REF" {input.r1} {input.r2} {threads} \
                "{params.sample}" "$ROUND_DIR/mapped.bam" "{params.aligner}"

            # ── Dedup if dedup line ──
            if [ "$LINE" = "dedup" ]; then
                picard MarkDuplicates \
                    I="$ROUND_DIR/mapped.bam" \
                    O="$ROUND_DIR/mapped.dedup.bam" \
                    M="$ROUND_DIR/dup_metrics.txt" \
                    REMOVE_DUPLICATES=true \
                    ASSUME_SORTED=true \
                    VALIDATION_STRINGENCY=SILENT
                samtools index "$ROUND_DIR/mapped.dedup.bam"
                ANALYSIS_BAM="$ROUND_DIR/mapped.dedup.bam"
            else
                ANALYSIS_BAM="$ROUND_DIR/mapped.bam"
            fi

            # ── Compute stats: mapped reads and completeness (≥2x) ──
            MAPPED=$(samtools view -c -F 4 "$ANALYSIS_BAM")
            TOTAL_BASES=$(awk '{{sum+=$2}} END{{print sum+0}}' "${{REF}}.fai")
            COVERED=$(samtools depth -a -J "$ANALYSIS_BAM" \
                | awk '$3 >= 2 {{n++}} END {{print n+0}}')
            if [ "$TOTAL_BASES" -gt 0 ]; then
                COMPLETENESS=$(awk "BEGIN {{printf \"%.2f\", 100.0 * $COVERED / $TOTAL_BASES}}")
            else
                COMPLETENESS="0.00"
            fi

            printf "%s\t%s\t%s\n" "$ROUND" "$MAPPED" "$COMPLETENESS" >> {output.stats}
            echo "Round $ROUND: mapped=$MAPPED, completeness=$COMPLETENESS%"

            # ── Call consensus ──
            bcftools mpileup -Ou -f "$REF" "$ANALYSIS_BAM" \
                | bcftools call -Ou -c --ploidy 1 \
                | bcftools norm -f "$REF" -Oz -o "$ROUND_DIR/variants.raw.vcf.gz"
            bcftools index "$ROUND_DIR/variants.raw.vcf.gz"

            if [ "$ENFORCE" = "True" ]; then
                bcftools filter \
                    -e 'INDEL && (abs(strlen(REF)-strlen(ALT)) % 3)!=0' \
                    -Oz -o "$ROUND_DIR/variants.inframe.vcf.gz" \
                    "$ROUND_DIR/variants.raw.vcf.gz"
                bcftools index "$ROUND_DIR/variants.inframe.vcf.gz"

                python3 {workflow.basedir}/scripts/filter_stop_variants.py \
                    "$REF" "$ROUND_DIR/variants.inframe.vcf.gz" \
                    "$ROUND_DIR/variants.filtered.vcf"
                bgzip -c "$ROUND_DIR/variants.filtered.vcf" > "$ROUND_DIR/variants.vcf.gz"
                rm "$ROUND_DIR/variants.filtered.vcf"
                bcftools index "$ROUND_DIR/variants.vcf.gz"

                bcftools consensus -f "$REF" "$ROUND_DIR/variants.vcf.gz" \
                    > "$ROUND_DIR/consensus.fasta"
            else
                bcftools consensus -f "$REF" "$ROUND_DIR/variants.raw.vcf.gz" \
                    > "$ROUND_DIR/consensus.fasta"
            fi

            # ── Optionally extend consensus using soft-clipped reads ──
            if [ "$EXTEND_REF" = "True" ]; then
                python3 {workflow.basedir}/scripts/extend_consensus.py \
                    "$ROUND_DIR/consensus.fasta" "$ANALYSIS_BAM" \
                    "$ROUND_DIR/consensus_extended.fasta"
                if [ -s "$ROUND_DIR/consensus_extended.fasta" ]; then
                    mv "$ROUND_DIR/consensus_extended.fasta" "$ROUND_DIR/consensus.fasta"
                fi
            fi

            LAST_CONSENSUS="$ROUND_DIR/consensus.fasta"
            LAST_BAM="$ANALYSIS_BAM"

            # ── Check convergence (only after min_iterations) ──
            if [ "$ROUND" -ge "$MIN_ITER" ]; then
                IMPROVED=$(python3 -c "print(1 if $MAPPED > $PREV_MAPPED and float($COMPLETENESS) > float($PREV_COMPLETENESS) else 0)")
                if [ "$IMPROVED" -eq 0 ]; then
                    echo "Converged after $ROUND rounds"
                    break
                fi
            fi

            if [ "$ROUND" -ge "$MAX_ITER" ]; then
                echo "Reached max iterations ($MAX_ITER)"
                break
            fi

            PREV_MAPPED=$MAPPED
            PREV_COMPLETENESS=$COMPLETENESS
            REF="$ROUND_DIR/consensus.fasta"
        done

        # Copy final outputs
        cp "$LAST_CONSENSUS" {output.consensus}
        cp "$LAST_BAM" {output.bam}
        cp "${{LAST_BAM}}.bai" {output.bai}
        """


rule dedup:
    """Remove PCR duplicates with Picard MarkDuplicates (generic rule)."""
    input:
        bam = "{path}/mapped.bam",
    output:
        bam     = "{path}/mapped.dedup.bam",
        bai     = "{path}/mapped.dedup.bam.bai",
        metrics = "{path}/dup_metrics.txt",
    shell:
        """
        set -euo pipefail
        picard MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT
        samtools index {output.bam}
        """


# ═══════════════════════════════════════════════════════════════
# DE NOVO ASSEMBLY (Trinity, SPAdes, MEGAHIT)
# ═══════════════════════════════════════════════════════════════

rule extract_assembly_reads:
    """Extract paired FASTQ from the final iterative-mapping BAM for assembly."""
    input:
        bam = SAMPLE_DIR + "/{sample}/{line}/iterative/analysis.bam",
    output:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    shell:
        """
        set -euo pipefail
        samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            {input.bam}
        """


rule trinity:
    """Run Trinity de novo assembly via Apptainer (use --use-singularity)."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = SAMPLE_DIR + "/{sample}/{line}/trinity/assembly/Trinity.fasta",
    params:
        memory = TRINITY_MEMORY,
    threads: THREADS
    singularity: "docker://trinityrnaseq/trinityrnaseq"
    shell:
        """
        set -euo pipefail
        TRINITY_OUT=$(dirname {output.fasta})/trinity_assembly

        rm -rf "$TRINITY_OUT"

        Trinity \
            --seqType fq \
            --left {input.r1} \
            --right {input.r2} \
            --max_memory {params.memory}G \
            --CPU {threads} \
            --output "$TRINITY_OUT"

        cp "$TRINITY_OUT.Trinity.fasta" {output.fasta}
        """


rule spades:
    """Run SPAdes in rnaviral mode."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = SAMPLE_DIR + "/{sample}/{line}/spades/contigs.fasta",
    params:
        outdir = SAMPLE_DIR + "/{sample}/{line}/spades",
        memory = SPADES_MEMORY,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        spades.py --rnaviral \
            -1 {input.r1} \
            -2 {input.r2} \
            -t {threads} \
            -m {params.memory} \
            -o {params.outdir}
        """


rule megahit:
    """Run MEGAHIT de novo assembly."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = SAMPLE_DIR + "/{sample}/{line}/megahit/final.contigs.fa",
    params:
        outdir = SAMPLE_DIR + "/{sample}/{line}/megahit",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -t {threads} \
            -o {params.outdir}
        """


rule coronaspades:
    """Run coronaSPAdes for coronavirus-like genome assembly."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = SAMPLE_DIR + "/{sample}/{line}/coronaspades/contigs.fasta",
    params:
        outdir = SAMPLE_DIR + "/{sample}/{line}/coronaspades",
        memory = SPADES_MEMORY,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        coronaspades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -t {threads} \
            -m {params.memory} \
            -o {params.outdir}
        # coronaSPAdes produces raw_contigs.fasta instead of contigs.fasta
        if [ ! -f {output.fasta} ] && [ -f {params.outdir}/raw_contigs.fasta ]; then
            cp {params.outdir}/raw_contigs.fasta {output.fasta}
        fi
        """


rule iva:
    """Run IVA (Iterative Virus Assembler)."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = SAMPLE_DIR + "/{sample}/{line}/iva/contigs.fasta",
    params:
        outdir = SAMPLE_DIR + "/{sample}/{line}/iva",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        iva \
            -f {input.r1} \
            -r {input.r2} \
            --threads {threads} \
            {params.outdir}
        """


# ═══════════════════════════════════════════════════════════════
# iVAR PER-CONTIG CONSENSUS
# ═══════════════════════════════════════════════════════════════

rule ivar_map:
    """Remap original reads to the iterative-mapping final consensus for iVar."""
    input:
        ref = SAMPLE_DIR + "/{sample}/{line}/iterative/final_consensus.fasta",
        r1  = get_r1,
        r2  = get_r2,
    output:
        bam = SAMPLE_DIR + "/{sample}/{line}/ivar/mapped.bam",
        bai = SAMPLE_DIR + "/{sample}/{line}/ivar/mapped.bam.bai",
    params:
        sample  = "{sample}",
        aligner = ALIGNER,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        bash {workflow.basedir}/scripts/map_reads.sh \
            {input.ref} {input.r1} {input.r2} {threads} \
            {params.sample} {output.bam} {params.aligner}
        cp "{input.ref}" "$(dirname {output.bam})/reference.fasta"
        """


rule ivar_consensus:
    """Run iVar per-contig consensus calling."""
    input:
        bam = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/ivar", wc.line
        ),
        ref = SAMPLE_DIR + "/{sample}/{line}/iterative/final_consensus.fasta",
    output:
        consensus = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
    shell:
        r"""
        set -euo pipefail
        TMP_DIR="$(dirname {output.consensus})/ivar_tmp"
        mkdir -p "$TMP_DIR"

        samtools faidx {input.ref}
        samtools idxstats {input.bam} > "$TMP_DIR/idxstats.txt"

        > {output.consensus}

        while IFS=$'\t' read -r CONTIG LENGTH MAPPED UNMAPPED; do
            if [ "$CONTIG" != "*" ] && [ "$MAPPED" -gt 0 ]; then
                samtools faidx {input.ref} "$CONTIG" > "$TMP_DIR/$CONTIG.fasta"
                samtools view -b {input.bam} "$CONTIG" > "$TMP_DIR/$CONTIG.bam"
                samtools index "$TMP_DIR/$CONTIG.bam"

                samtools mpileup -aa -A -d 0 -Q 0 -r "$CONTIG" \
                    -f "$TMP_DIR/$CONTIG.fasta" "$TMP_DIR/$CONTIG.bam" \
                    | ivar consensus -m 2 -p "$TMP_DIR/$CONTIG"

                if [ -f "$TMP_DIR/$CONTIG.fa" ]; then
                    cat "$TMP_DIR/$CONTIG.fa" >> {output.consensus}
                else
                    SEQ_LEN=$(awk 'NR>1' "$TMP_DIR/$CONTIG.fasta" | tr -d '\n' | wc -c)
                    python3 -c "print('>${{CONTIG}}'); print('N' * $SEQ_LEN)" \
                        >> {output.consensus}
                fi
            fi
        done < "$TMP_DIR/idxstats.txt"

        rm -rf "$TMP_DIR"
        """


# ═══════════════════════════════════════════════════════════════
# FINAL MAPPING
# ═══════════════════════════════════════════════════════════════

rule final_map:
    """Map reads to the iVar final consensus."""
    input:
        ref = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        r1  = get_r1,
        r2  = get_r2,
    output:
        bam = SAMPLE_DIR + "/{sample}/{line}/final/mapped.bam",
        bai = SAMPLE_DIR + "/{sample}/{line}/final/mapped.bam.bai",
    params:
        sample  = "{sample}",
        aligner = ALIGNER,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        bash {workflow.basedir}/scripts/map_reads.sh \
            {input.ref} {input.r1} {input.r2} {threads} \
            {params.sample} {output.bam} {params.aligner}
        cp "{input.ref}" "$(dirname {output.bam})/reference.fasta"
        """


# ═══════════════════════════════════════════════════════════════
# ASSEMBLY-GUIDED POLISHING
# ═══════════════════════════════════════════════════════════════

rule assembly_guided_polish:
    """Polish final consensus using a de novo assembly."""
    input:
        assembly  = get_assembly_fasta,
        consensus = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
    output:
        polished = SAMPLE_DIR + "/{sample}/{line}/{assembler}_polish/{assembler}_guided_polished.fasta",
    threads: THREADS
    shell:
        """
        set -euo pipefail
        DIR=$(dirname {output.polished})

        samtools faidx {input.consensus}

        minimap2 -ax asm5 -t {threads} {input.consensus} {input.assembly} \
            | samtools sort -o "$DIR/assembly_vs_consensus.bam"
        samtools index "$DIR/assembly_vs_consensus.bam"

        bcftools mpileup -Ou -f {input.consensus} "$DIR/assembly_vs_consensus.bam" \
            | bcftools call -Ou -c --ploidy 1 \
            | bcftools norm -f {input.consensus} -Oz -o "$DIR/variants.vcf.gz"
        bcftools index "$DIR/variants.vcf.gz"

        bcftools consensus -f {input.consensus} "$DIR/variants.vcf.gz" \
            > {output.polished}
        """


# ═══════════════════════════════════════════════════════════════
# ASSEMBLY-TO-CONSENSUS COMPARISON
# ═══════════════════════════════════════════════════════════════

rule assembly_vs_consensus:
    """Compare a de novo assembly to the final consensus."""
    input:
        assembly  = get_assembly_fasta,
        consensus = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
    output:
        paf = SAMPLE_DIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.paf",
        bam = SAMPLE_DIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.bam",
        vcf = SAMPLE_DIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.vcf.gz",
    threads: THREADS
    shell:
        """
        set -euo pipefail

        samtools faidx {input.consensus}

        minimap2 -ax asm5 --paf-no-hit -t {threads} \
            {input.consensus} {input.assembly} \
            | samtools sort -o {output.bam}
        samtools index {output.bam}

        samtools view -h {output.bam} | paftools.js sam2paf - > {output.paf}

        bcftools mpileup -Ou -f {input.consensus} {output.bam} \
            | bcftools call -Ou -c --ploidy 1 \
            | bcftools norm -f {input.consensus} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """


# ═══════════════════════════════════════════════════════════════
# FINAL ALIGNMENT VS INITIAL REFERENCE
# ═══════════════════════════════════════════════════════════════

rule final_alignment:
    """Align all consensus sequences against the initial reference, per contig."""
    input:
        ref      = SAMPLE_DIR + "/{sample}/reference/ref.fasta",
        iter_stats = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        ivar     = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
    output:
        alignment = SAMPLE_DIR + "/{sample}/{line}/final_alignment/aligned.fasta",
    params:
        assemblers  = ASSEMBLERS,
        sample_dir  = SAMPLE_DIR,
    threads: THREADS
    run:
        import subprocess, os, glob, re, tempfile

        # Discover round directories (created by iterative_mapping)
        iter_base = os.path.join(
            str(params.sample_dir), wildcards.sample, wildcards.line, "iterative"
        )
        round_dirs = sorted(
            [d for d in glob.glob(os.path.join(iter_base, "round_*")) if os.path.isdir(d)],
            key=lambda d: int(re.search(r"round_(\d+)", d).group(1)),
        )

        # Collect all FASTAs with labels
        steps = [("reference", str(input.ref))]
        for rd in round_dirs:
            n = re.search(r"round_(\d+)", rd).group(1)
            cons_path = os.path.join(rd, "consensus.fasta")
            if os.path.exists(cons_path):
                steps.append((f"round_{n}", cons_path))
        steps.append(("ivar_consensus", str(input.ivar)))
        for i, assembler in enumerate(params.assemblers):
            steps.append((f"{assembler}_polished", str(input.polished[i])))

        # Read all sequences
        all_seqs = {}
        for label, path in steps:
            all_seqs[label] = read_fasta(path)

        # Get contig names from reference
        ref_contigs = list(all_seqs["reference"].keys())

        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        # Clear output file
        open(str(output.alignment), "w").close()

        # For each contig, collect sequences and align with MAFFT
        for contig in ref_contigs:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fa", delete=False
            ) as tmp_in:
                for label, _ in steps:
                    seqs = all_seqs[label]
                    seq = None
                    for cname, cseq in seqs.items():
                        if contig in cname or cname in contig:
                            seq = cseq
                            break
                    if seq:
                        tmp_in.write(f">{contig}__{label}\n{seq}\n")
                tmp_in_path = tmp_in.name

            # Run MAFFT
            result = subprocess.run(
                ["mafft", "--auto", "--thread", str(threads), tmp_in_path],
                capture_output=True, text=True,
            )

            # Append to output
            with open(str(output.alignment), "a") as out:
                out.write(result.stdout)

            os.unlink(tmp_in_path)

        print(f"Alignment written to {output.alignment}")


# ═══════════════════════════════════════════════════════════════
# HOMOPLASY ANALYSIS
# ═══════════════════════════════════════════════════════════════

rule homoplasy_analysis:
    """Compute delta score (Holland et al. 2002) and homoplasy metrics from MAFFT alignment."""
    input:
        alignment = SAMPLE_DIR + "/{sample}/{line}/final_alignment/aligned.fasta",
    output:
        report    = SAMPLE_DIR + "/{sample}/{line}/homoplasy/homoplasy_report.tsv",
        distances = SAMPLE_DIR + "/{sample}/{line}/homoplasy/pairwise_distances.tsv",
    run:
        import os
        from itertools import combinations

        def read_aligned_fasta(path):
            seqs = []
            name = None
            parts = []
            with open(path) as f:
                for ln in f:
                    ln = ln.strip()
                    if ln.startswith(">"):
                        if name:
                            seqs.append((name, "".join(parts).upper()))
                        name = ln[1:].strip()
                        parts = []
                    else:
                        parts.append(ln)
            if name:
                seqs.append((name, "".join(parts).upper()))
            return seqs

        def group_by_contig(seqs):
            """Group sequences by contig (part before __)."""
            contigs = {}
            for name, seq in seqs:
                if "__" in name:
                    contig, label = name.split("__", 1)
                else:
                    contig, label = name, name
                contigs.setdefault(contig, {})[label] = seq
            return contigs

        def hamming_prop(s1, s2):
            """Proportion of differing positions (excluding gap-vs-gap)."""
            diffs = 0
            compared = 0
            for a, b in zip(s1, s2):
                if a == "-" and b == "-":
                    continue
                compared += 1
                if a != b:
                    diffs += 1
            return diffs / compared if compared > 0 else 0.0

        def delta_score(n, dist):
            """Holland et al. 2002 delta score for treelikeness.
            0 = perfectly tree-like (no homoplasy).
            1 = maximally conflicting signal.
            """
            if n < 4:
                return float("nan")
            deltas = []
            for i, j, k, l in combinations(range(n), 4):
                sums = sorted([
                    dist[i][j] + dist[k][l],
                    dist[i][k] + dist[j][l],
                    dist[i][l] + dist[j][k],
                ])
                denom = sums[2] - sums[0]
                deltas.append((sums[2] - sums[1]) / denom if denom > 0 else 0.0)
            return sum(deltas) / len(deltas)

        # ── Parse alignment ──
        all_seqs = read_aligned_fasta(str(input.alignment))
        contigs = group_by_contig(all_seqs)

        os.makedirs(os.path.dirname(str(output.report)), exist_ok=True)

        report_rows = []
        dist_rows = []

        for contig, seqs_dict in contigs.items():
            labels = list(seqs_dict.keys())
            sequences = [seqs_dict[l] for l in labels]
            n = len(labels)
            aln_len = len(sequences[0]) if sequences else 0

            # Pairwise Hamming distances
            dist = [[0.0] * n for _ in range(n)]
            for i, j in combinations(range(n), 2):
                d = hamming_prop(sequences[i], sequences[j])
                dist[i][j] = d
                dist[j][i] = d
                seg = contig.split("_", 1)[0] if "_" in contig else contig
                dist_rows.append((seg, contig, labels[i], labels[j], f"{d:.6f}"))

            # Variable and parsimony-informative sites
            var_sites = 0
            pi_sites = 0
            for pos in range(aln_len):
                bases = [s[pos] for s in sequences if s[pos] != "-"]
                if not bases:
                    continue
                unique = set(bases)
                if len(unique) > 1:
                    var_sites += 1
                    counts = {}
                    for b in bases:
                        counts[b] = counts.get(b, 0) + 1
                    if sum(1 for c in counts.values() if c >= 2) >= 2:
                        pi_sites += 1

            # Delta score
            ds = delta_score(n, dist)

            # Mean pairwise distance
            pair_dists = [dist[i][j] for i, j in combinations(range(n), 2)]
            mean_dist = sum(pair_dists) / len(pair_dists) if pair_dists else 0.0

            # Extract segment name (L, M, S) from contig like "L_B1"
            segment = contig.split("_", 1)[0] if "_" in contig else contig

            report_rows.append({
                "segment": segment,
                "contig": contig,
                "num_sequences": n,
                "alignment_length": aln_len,
                "variable_sites": var_sites,
                "parsimony_informative_sites": pi_sites,
                "delta_score": f"{ds:.4f}" if ds == ds else "NA",
                "mean_pairwise_distance": f"{mean_dist:.6f}",
            })

        # ── Write report ──
        hdr = ["sample", "line", "segment", "contig", "num_sequences",
               "alignment_length", "variable_sites",
               "parsimony_informative_sites",
               "delta_score", "mean_pairwise_distance"]
        with open(str(output.report), "w") as f:
            f.write("\t".join(hdr) + "\n")
            for r in report_rows:
                f.write("\t".join([
                    wildcards.sample, wildcards.line,
                    r["segment"], r["contig"],
                    str(r["num_sequences"]),
                    str(r["alignment_length"]), str(r["variable_sites"]),
                    str(r["parsimony_informative_sites"]),
                    r["delta_score"], r["mean_pairwise_distance"],
                ]) + "\n")

        # ── Write pairwise distances ──
        with open(str(output.distances), "w") as f:
            f.write("\t".join(["segment", "contig", "method_1", "method_2", "hamming_distance"]) + "\n")
            for row in dist_rows:
                f.write("\t".join(row) + "\n")


# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

rule summary:
    """Write per-segment summary with coverage, completeness, N counts, and stop-codon check."""
    input:
        ref        = SAMPLE_DIR + "/{sample}/reference/ref.fasta",
        iter_stats = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        consensus  = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished   = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
        comparison = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.paf",
            assembler=ASSEMBLERS,
        ),
        alignment  = SAMPLE_DIR + "/{sample}/{line}/final_alignment/aligned.fasta",
        final_bam  = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        ivar_bam = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/ivar", wc.line
        ),
    output:
        summary = SAMPLE_DIR + "/{sample}/{line}/summary.tsv",
    params:
        check_stops = CHECK_STOPS,
        assemblers  = ASSEMBLERS,
        sample_dir  = SAMPLE_DIR,
    run:
        import subprocess, glob, os, re

        GENETIC_CODE = {
            "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
            "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
            "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
            "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
            "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
            "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
            "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
            "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
            "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
            "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
            "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
            "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
            "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
            "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
            "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
            "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
        }

        def get_coverage(bam_path):
            """Run samtools coverage → dict of contig → {numreads, coverage, meandepth}."""
            result = subprocess.run(
                ["samtools", "coverage", str(bam_path)],
                capture_output=True, text=True,
            )
            stats = {}
            for ln in result.stdout.strip().split("\n"):
                if ln.startswith("#"):
                    continue
                f = ln.split("\t")
                stats[f[0]] = {
                    "numreads": int(f[3]),
                    "coverage": float(f[5]),
                    "meandepth": float(f[6]),
                }
            return stats

        def find_match(contig, name_dict):
            """Find matching key by exact match or substring containment."""
            if contig in name_dict:
                return name_dict[contig]
            for key, val in name_dict.items():
                if contig in key or key in contig:
                    return val
            return None

        def has_internal_stop(seq):
            aa = [GENETIC_CODE.get(seq[i:i+3], "X")
                  for i in range(0, len(seq) - 2, 3)]
            return "*" in aa[:-1]

        # ── Reference contigs (defines segment order) ──
        ref_seqs = read_fasta(str(input.ref))
        ref_contigs = list(ref_seqs.keys())

        # ── Discover round directories (created by iterative_mapping) ──
        iter_base = os.path.join(
            str(params.sample_dir), wildcards.sample, wildcards.line, "iterative"
        )
        round_dirs = sorted(
            [d for d in glob.glob(os.path.join(iter_base, "round_*")) if os.path.isdir(d)],
            key=lambda d: int(re.search(r"round_(\d+)", d).group(1)),
        )

        # ── BAM coverage per step ──
        bam_steps = []
        for rd in round_dirs:
            n = re.search(r"round_(\d+)", rd).group(1)
            if wildcards.line == "dedup":
                bam_path = os.path.join(rd, "mapped.dedup.bam")
            else:
                bam_path = os.path.join(rd, "mapped.bam")
            if os.path.exists(bam_path):
                bam_steps.append((f"round_{n}", bam_path))
        bam_steps.append(("ivar", str(input.ivar_bam)))
        bam_steps.append(("final", str(input.final_bam)))

        bam_cov = {}
        for step, path in bam_steps:
            bam_cov[step] = get_coverage(path)

        # ── Consensus N counts per step ──
        cons_steps = []
        for rd in round_dirs:
            n = re.search(r"round_(\d+)", rd).group(1)
            cons_path = os.path.join(rd, "consensus.fasta")
            if os.path.exists(cons_path):
                cons_steps.append((f"round_{n}", cons_path))
        cons_steps.append(("ivar", str(input.consensus)))
        for i, assembler in enumerate(params.assemblers):
            cons_steps.append((f"{assembler}_polished", str(input.polished[i])))

        cons_seqs = {}
        for step, path in cons_steps:
            cons_seqs[step] = read_fasta(path)

        # ── Build header ──
        header = ["sample", "line", "segment", "segment_length"]
        for step, _ in bam_steps:
            header.extend([
                f"{step}_mapped_reads",
                f"{step}_mean_depth",
                f"{step}_completeness",
            ])
        for step, _ in cons_steps:
            header.extend([f"{step}_N", f"{step}_N_frac"])
        header.append("internal_stop")

        # ── Build one row per segment ──
        rows = []
        for contig in ref_contigs:
            seg_len = len(ref_seqs[contig])
            row = [wildcards.sample, wildcards.line, contig, str(seg_len)]

            # BAM stats
            for step, _ in bam_steps:
                hit = find_match(contig, bam_cov[step])
                if hit:
                    row.append(str(hit["numreads"]))
                    row.append(f"{hit['meandepth']:.1f}")
                    row.append(f"{hit['coverage']:.1f}")
                else:
                    row.extend(["0", "0.0", "0.0"])

            # Consensus N stats
            for step, _ in cons_steps:
                seqs = cons_seqs[step]
                hit = find_match(contig, seqs)
                if hit:
                    n_count = hit.count("N")
                    n_frac = 100.0 * n_count / len(hit) if len(hit) > 0 else 0
                    row.append(str(n_count))
                    row.append(f"{n_frac:.1f}")
                else:
                    row.extend(["NA", "NA"])

            # Internal stop codon check (on ivar consensus)
            ivar_seq = find_match(contig, cons_seqs.get("ivar", {}))
            stop = 0
            if ivar_seq:
                stop = int(has_internal_stop(ivar_seq))
                if stop:
                    print(f"[WARN] {contig}: internal stop codon in ivar consensus")
            row.append(str(stop))

            rows.append(row)

        with open(str(output.summary), "w") as out:
            out.write("\t".join(header) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")


# ═══════════════════════════════════════════════════════════════
# IQ-TREE TWO-PHASE PHYLOGENETIC ANALYSIS
# Phase 1: backbone tree from GenBank references (per segment)
# Phase 2: place each sample's consensus sequences onto backbone
# ═══════════════════════════════════════════════════════════════

rule iqtree_backbone:
    """Phase 1: Build backbone ML tree from GenBank references only (once per segment)."""
    input:
        genbank = COMMON_DIR + "/genbank_subsampled/{segment}.subsampled.fa",
    output:
        alignment = COMMON_DIR + "/iqtree/backbone/{segment}/aligned.fasta",
        treefile  = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.treefile",
        iqtree    = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.iqtree",
        log       = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.log",
    threads: THREADS
    run:
        import subprocess, os

        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        genbank_seqs = read_fasta(str(input.genbank))

        if len(genbank_seqs) < 4:
            with open(str(output.alignment), "w") as f:
                for name, seq in genbank_seqs.items():
                    f.write(f">{name}\n{seq}\n")
            with open(str(output.treefile), "w") as f:
                f.write("(no_tree);\n")
            with open(str(output.iqtree), "w") as f:
                f.write("# Too few sequences for IQ-TREE backbone\n")
            with open(str(output.log), "w") as f:
                f.write("# Too few sequences for IQ-TREE backbone\n")
            print(f"[WARN] {wildcards.segment}: only {len(genbank_seqs)} GenBank refs, skipping backbone")
            return

        # Write unaligned GenBank FASTA
        unaligned_path = os.path.join(out_dir, "unaligned.fasta")
        with open(unaligned_path, "w") as f:
            for name, seq in genbank_seqs.items():
                f.write(f">{name}\n{seq}\n")

        print(f"[INFO] {wildcards.segment}: {len(genbank_seqs)} GenBank refs for backbone")

        # MAFFT alignment
        result = subprocess.run(
            ["mafft", "--auto", "--thread", str(threads), unaligned_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr}")
        with open(str(output.alignment), "w") as f:
            f.write(result.stdout)

        # IQ-TREE full ML search with ModelFinder
        prefix = os.path.join(out_dir, "backbone")
        cmd = [
            "iqtree",
            "-s", str(output.alignment),
            "-m", "MFP",          # ModelFinder Plus
            "-bb", "1000",        # ultrafast bootstrap
            "-alrt", "1000",      # SH-aLRT branch test
            "-T", str(threads),
            "--prefix", prefix,
            "-redo",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"IQ-TREE backbone failed: {result.stderr}")

        print(f"[INFO] {wildcards.segment}: backbone tree completed → {output.treefile}")


rule iqtree_place_sample:
    """Phase 2: Place a sample's consensus sequences onto the fixed backbone tree."""
    input:
        backbone_tree = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.treefile",
        backbone_aln  = COMMON_DIR + "/iqtree/backbone/{segment}/aligned.fasta",
        backbone_iq   = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.iqtree",
        ivar = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
        iter_stats = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
    output:
        alignment  = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/aligned.fasta",
        treefile   = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.treefile",
        iqtree     = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.iqtree",
        log        = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.log",
        placements = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/consensus_placements.tsv",
    params:
        assemblers = ASSEMBLERS,
        sample_dir = SAMPLE_DIR,
    threads: THREADS
    run:
        import subprocess, os, glob, re, shutil

        def extract_segment(fasta_path, segment, sample):
            """Extract the sequence for a given segment from a multi-segment FASTA."""
            seqs = read_fasta(fasta_path)
            target = f"{segment}_{sample}"
            if target in seqs:
                return seqs[target]
            for name, seq in seqs.items():
                if f"_{segment}_{sample}" in name or name.startswith(f"{segment}_"):
                    return seq
            return None

        def run_single_placement(taxon_name, seq, backbone_aln, backbone_tree,
                                 model, sub_dir, n_threads):
            """Place one consensus sequence on the backbone tree.

            Returns (iq_path, tree_path, aln_path, log_path) or None on failure.
            """
            os.makedirs(sub_dir, exist_ok=True)

            # Write single-sequence FASTA
            query_fasta = os.path.join(sub_dir, "query.fasta")
            with open(query_fasta, "w") as fh:
                fh.write(f">{taxon_name}\n{seq}\n")

            # MAFFT --add: add this one sequence to the backbone alignment
            res = subprocess.run(
                ["mafft", "--add", query_fasta, "--keeplength",
                 "--thread", str(n_threads), backbone_aln],
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                print(f"[WARN] MAFFT --add failed for {taxon_name}: {res.stderr}")
                return None

            aln_path = os.path.join(sub_dir, "aligned.fasta")
            with open(aln_path, "w") as fh:
                fh.write(res.stdout)

            # IQ-TREE with constraint backbone
            prefix = os.path.join(sub_dir, "iqtree")
            cmd = [
                "iqtree",
                "-s", aln_path,
                "-g", backbone_tree,
                "-m", model,
                "-T", str(n_threads),
                "--prefix", prefix,
                "-redo",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                print(f"[WARN] IQ-TREE failed for {taxon_name}: {res.stderr}")
                return None

            return (
                os.path.join(sub_dir, "iqtree.iqtree"),
                os.path.join(sub_dir, "iqtree.treefile"),
                aln_path,
                os.path.join(sub_dir, "iqtree.log"),
            )

        seg = wildcards.segment
        line = wildcards.line
        sample = wildcards.sample
        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        # Check backbone is valid
        with open(str(input.backbone_tree)) as f:
            tree_str = f.read().strip()
        if tree_str.startswith("(no_tree)") or not tree_str:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# No backbone tree available\n")
            with open(str(output.placements), "w") as f:
                f.write("# No backbone tree available\n")
            print(f"[WARN] {seg}/{line}/{sample}: no backbone tree, skipping placement")
            return

        # Parse model from backbone IQ-TREE report
        best_model = "MFP"
        with open(str(input.backbone_iq)) as f:
            for bline in f:
                m = re.search(r"Best-fit model according to BIC:\s+(.+)", bline)
                if m:
                    best_model = m.group(1).strip()
                    break

        # ── Collect sample consensus sequences ──
        # Dict of method_label -> (taxon_name, sequence)
        collected = {}

        # Reference used for iterative mapping
        ref_path = os.path.join(
            str(params.sample_dir), sample, line, "iterative", "round_1", "reference.fasta"
        )
        if os.path.exists(ref_path):
            seq = extract_segment(ref_path, seg, sample)
            if seq and seq.replace("N", ""):
                collected["reference"] = (f"{sample}_reference", seq)

        # Iterative mapping rounds
        iter_base = os.path.join(str(params.sample_dir), sample, line, "iterative")
        round_dirs = sorted(
            [d for d in glob.glob(os.path.join(iter_base, "round_*")) if os.path.isdir(d)],
            key=lambda d: int(re.search(r"round_(\d+)", d).group(1)),
        )
        for rd in round_dirs:
            n = re.search(r"round_(\d+)", rd).group(1)
            cons_path = os.path.join(rd, "consensus.fasta")
            if os.path.exists(cons_path):
                seq = extract_segment(cons_path, seg, sample)
                if seq and seq.replace("N", ""):
                    collected[f"round_{n}"] = (f"{sample}_round_{n}", seq)

        # iVar consensus
        if os.path.exists(str(input.ivar)):
            seq = extract_segment(str(input.ivar), seg, sample)
            if seq and seq.replace("N", ""):
                collected["ivar"] = (f"{sample}_ivar", seq)

        # Polished consensus per assembler
        for i, assembler in enumerate(params.assemblers):
            pol_path = str(input.polished[i])
            if os.path.exists(pol_path):
                seq = extract_segment(pol_path, seg, sample)
                if seq and seq.replace("N", ""):
                    collected[assembler] = (f"{sample}_{assembler}", seq)

        if len(collected) == 0:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# No sample sequences for this segment\n")
            with open(str(output.placements), "w") as f:
                f.write("# No sample sequences for this segment\n")
            print(f"[WARN] {seg}/{line}/{sample}: no sample sequences, skipping")
            return

        print(f"[INFO] {seg}/{line}/{sample}: {len(collected)} consensus sequences "
              f"to place individually on backbone")

        # ── Individual placements ──
        placements_dir = os.path.join(out_dir, "placements")
        placement_results = []  # (method, taxon_name, iq_path, tree_path, aln_path)

        for method, (taxon_name, seq) in collected.items():
            sub_dir = os.path.join(placements_dir, method)
            paths = run_single_placement(
                taxon_name, seq,
                str(input.backbone_aln), str(input.backbone_tree),
                best_model, sub_dir, threads,
            )
            if paths:
                iq_path, tree_path, aln_path, log_path = paths
                placement_results.append((method, taxon_name, iq_path, tree_path, aln_path))
                print(f"[INFO] {seg}/{line}/{sample}/{method}: individual placement done")
            else:
                print(f"[WARN] {seg}/{line}/{sample}/{method}: placement failed, skipping")

        # ── Populate top-level outputs for backward compatibility ──
        # Prefer the ivar placement; fall back to last successful placement
        main_method = None
        for pref in ["ivar"]:
            for m, tn, iq, tr, al in placement_results:
                if m == pref:
                    main_method = (m, tn, iq, tr, al)
                    break
            if main_method:
                break
        if not main_method and placement_results:
            main_method = placement_results[-1]

        if main_method:
            _, _, main_iq, main_tree, main_aln = main_method
            main_log = main_iq.replace("iqtree.iqtree", "iqtree.log")
            shutil.copy2(main_aln, str(output.alignment))
            shutil.copy2(main_tree, str(output.treefile))
            shutil.copy2(main_iq, str(output.iqtree))
            if os.path.exists(main_log):
                shutil.copy2(main_log, str(output.log))
            else:
                with open(str(output.log), "w") as f:
                    f.write("# See placements/ subdirectory\n")
        else:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# All individual placements failed\n")

        # ── Write consensus_placements.tsv ──
        placement_header = [
            "sample", "line", "segment", "method", "taxon_name",
            "tree_path", "alignment_path", "iqtree_path",
        ]
        with open(str(output.placements), "w") as f:
            f.write("\t".join(placement_header) + "\n")
            for method, taxon_name, iq_path, tree_path, aln_path in placement_results:
                f.write("\t".join([
                    sample, line, seg, method, taxon_name,
                    tree_path, aln_path, iq_path,
                ]) + "\n")

        print(f"[INFO] {seg}/{line}/{sample}: all placements completed "
              f"({len(placement_results)} individual)")


rule iqtree_homoplasy_summary:
    """Parse per-sample IQ-TREE outputs and compute CI/RI/HI via Fitch parsimony."""
    input:
        iqtree_files = expand(
            SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.iqtree",
            line=LINES,
            segment=SEGMENTS,
            sample=SAMPLES.keys(),
        ),
        tree_files = expand(
            SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.treefile",
            line=LINES,
            segment=SEGMENTS,
            sample=SAMPLES.keys(),
        ),
        aln_files = expand(
            SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/aligned.fasta",
            line=LINES,
            segment=SEGMENTS,
            sample=SAMPLES.keys(),
        ),
        placement_files = expand(
            SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/consensus_placements.tsv",
            line=LINES,
            segment=SEGMENTS,
            sample=SAMPLES.keys(),
        ),
    output:
        summary           = COMMON_DIR + "/iqtree/iqtree_summary.tsv",
        taxon_report      = COMMON_DIR + "/iqtree/taxon_homoplasy.tsv",
        individual_report = COMMON_DIR + "/iqtree/individual_placement_homoplasy.tsv",
    run:
        import os, re

        # ── Combined analysis: global metrics + per-taxon ──
        def analyse_tree(tree_path, aln_path):
            """Return (global_metrics_dict, per_taxon_stats_dict)."""
            na_global = {
                "parsimony_score": "NA", "ci": "NA", "ri": "NA", "hi": "NA",
            }
            with open(tree_path) as f:
                tree_str = f.read().strip()
            if tree_str.startswith("(no_tree)") or tree_str.startswith("# No") or not tree_str:
                return na_global, {}

            tree = Phylo.read(StringIO(tree_str), "newick")
            aln = read_fasta(aln_path)
            if not aln:
                return na_global, {}

            terminals = [t.name for t in tree.get_terminals()]
            parent_map = build_parent_map(tree)
            aln_len = len(next(iter(aln.values())))

            # ── Precompute pairwise patristic distances ──
            patristic = {}
            for i, n1 in enumerate(terminals):
                for n2 in terminals[i + 1:]:
                    d = tree.distance(n1, n2)
                    patristic[(n1, n2)] = d
                    patristic[(n2, n1)] = d
            all_pat_dists = sorted(patristic.values())
            median_pat = all_pat_dists[len(all_pat_dists) // 2] if all_pat_dists else 0.0

            # Per-taxon accumulators
            taxon_stats = {
                name: {
                    "terminal_changes": 0,
                    "convergent": 0,
                    "autapomorphic": 0,
                    "near_convergent": 0,
                    "distant_convergent": 0,
                    "conv_distances": [],
                }
                for name in terminals
            }
            total_variable = 0

            # Global accumulators
            total_parsimony = 0
            total_min = 0
            total_max = 0

            for pos in range(aln_len):
                states = {}
                for name in terminals:
                    if name in aln:
                        s = aln[name][pos]
                        if s not in ("-", "N", "?"):
                            states[name] = s

                unique = set(states.values())
                if len(unique) <= 1:
                    continue

                total_variable += 1

                # --- Fitch bottom-up ---
                node_cache = {}
                changes = fitch_bottom_up(tree.root, states, node_cache)
                total_parsimony += changes
                total_min += len(unique) - 1
                counts = Counter(states.values())
                total_max += len(states) - max(counts.values())

                # --- Fitch top-down assignment ---
                assigned = {}
                fitch_top_down(tree.root, node_cache, assigned)

                # --- Per-taxon: check terminal branch ---
                for leaf in tree.get_terminals():
                    if leaf.name not in states:
                        continue
                    parent = parent_map.get(id(leaf))
                    if parent is None:
                        continue
                    leaf_state = states[leaf.name]
                    parent_state = assigned.get(id(parent))
                    if leaf_state != parent_state:
                        taxon_stats[leaf.name]["terminal_changes"] += 1
                        # Find other taxa sharing the same derived state
                        others = [
                            n for n, s in states.items()
                            if n != leaf.name and s == leaf_state
                        ]
                        if others:
                            taxon_stats[leaf.name]["convergent"] += 1
                            nearest_dist = min(
                                patristic.get((leaf.name, n), float("inf"))
                                for n in others
                            )
                            taxon_stats[leaf.name]["conv_distances"].append(nearest_dist)
                            if nearest_dist < median_pat:
                                taxon_stats[leaf.name]["near_convergent"] += 1
                            else:
                                taxon_stats[leaf.name]["distant_convergent"] += 1
                        else:
                            taxon_stats[leaf.name]["autapomorphic"] += 1

            # Global metrics
            ci = total_min / total_parsimony if total_parsimony > 0 else float("nan")
            ri_denom = total_max - total_min
            ri = (total_max - total_parsimony) / ri_denom if ri_denom > 0 else float("nan")
            hi = 1.0 - ci if ci == ci else float("nan")

            global_metrics = {
                "parsimony_score": str(total_parsimony),
                "ci": f"{ci:.4f}" if ci == ci else "NA",
                "ri": f"{ri:.4f}" if ri == ri else "NA",
                "hi": f"{hi:.4f}" if hi == hi else "NA",
            }

            # Summarise per-taxon stats
            for name in taxon_stats:
                ts = taxon_stats[name]
                ts["total_variable"] = total_variable
                tc = ts["terminal_changes"]
                ts["homoplasy_frac"] = ts["convergent"] / tc if tc > 0 else 0.0
                dists = ts["conv_distances"]
                ts["mean_conv_dist"] = sum(dists) / len(dists) if dists else 0.0
                del ts["conv_distances"]

            return global_metrics, taxon_stats

        # ── IQ-TREE report parser ──
        def parse_iqtree(path):
            info = {
                "num_sequences": "NA",
                "alignment_length": "NA",
                "constant_sites": "NA",
                "constant_sites_pct": "NA",
                "variable_sites": "NA",
                "parsimony_informative": "NA",
                "parsimony_informative_pct": "NA",
                "singleton_sites": "NA",
                "best_model": "NA",
                "log_likelihood": "NA",
                "tree_length": "NA",
            }

            with open(path) as f:
                content = f.read()
            if content.startswith("# Too few") or content.startswith("# No"):
                return info

            # v2: "Number of sequences: N" / v3: "Input data: N sequences with M nucleotide sites"
            m = re.search(r"Number of sequences:\s+(\d+)", content)
            if m:
                info["num_sequences"] = m.group(1)
            else:
                m = re.search(r"Input data:\s+(\d+)\s+sequences\s+with\s+(\d+)\s+nucleotide", content)
                if m:
                    info["num_sequences"] = m.group(1)
                    info["alignment_length"] = m.group(2)

            # v2: "Number of columns: N"
            m = re.search(r"Number of columns:\s+(\d+)", content)
            if m:
                info["alignment_length"] = m.group(1)

            m = re.search(r"Number of constant sites:\s+(\d+)\s+\(=\s+([\d.]+)%", content)
            if m:
                info["constant_sites"] = m.group(1)
                info["constant_sites_pct"] = m.group(2)

            m = re.search(r"Number of parsimony informative sites:\s+(\d+)", content)
            if m:
                info["parsimony_informative"] = m.group(1)
                try:
                    pct = 100.0 * int(m.group(1)) / int(info["alignment_length"])
                    info["parsimony_informative_pct"] = f"{pct:.1f}"
                except (ValueError, TypeError, ZeroDivisionError):
                    pass

            m = re.search(r"Number of singleton sites:\s+(\d+)", content)
            if m:
                info["singleton_sites"] = m.group(1)

            # Variable sites = parsimony informative + singleton, or alignment_length - constant
            try:
                info["variable_sites"] = str(
                    int(info["parsimony_informative"]) + int(info["singleton_sites"])
                )
            except (ValueError, TypeError):
                try:
                    info["variable_sites"] = str(
                        int(info["alignment_length"]) - int(info["constant_sites"])
                    )
                except (ValueError, TypeError):
                    pass

            m = re.search(r"Best-fit model according to BIC:\s+(.+)", content)
            if m:
                info["best_model"] = m.group(1).strip()

            m = re.search(r"Log-likelihood of the tree:\s+([-\d.]+)", content)
            if m:
                info["log_likelihood"] = m.group(1)

            m = re.search(r"Total tree length.*?:\s+([\d.]+)", content)
            if m:
                info["tree_length"] = m.group(1)

            return info

        # ── Build summary from per-sample trees ──
        # Group inputs by (line, segment, sample)
        path_map = {}
        for iq_path in input.iqtree_files:
            parts = str(iq_path).split("/")
            # Path: .../samples/{sample}/{line}/iqtree/{segment}/iqtree.iqtree
            segment = parts[-2]
            # parts[-3] == "iqtree"
            line = parts[-4]
            sample_name = parts[-5]
            base_dir = os.path.dirname(str(iq_path))
            path_map[(line, segment, sample_name)] = base_dir

        summary_rows = []
        all_taxon_rows = []

        for (line, segment, sample_name), base_dir in sorted(path_map.items()):
            iq_path = os.path.join(base_dir, "iqtree.iqtree")
            tree_path = os.path.join(base_dir, "iqtree.treefile")
            aln_path = os.path.join(base_dir, "aligned.fasta")

            info = parse_iqtree(iq_path)
            global_metrics, taxon_stats = analyse_tree(tree_path, aln_path)
            info.update(global_metrics)
            summary_rows.append((line, segment, sample_name, info))

            # ── Per-taxon rows: only report sample taxa (not GenBank backbone) ──
            genbank_agg = {
                "terminal_changes": 0, "convergent": 0,
                "autapomorphic": 0, "near_convergent": 0,
                "distant_convergent": 0, "conv_dist_sum": 0.0,
                "conv_dist_n": 0, "count": 0,
            }

            for taxon, ts in sorted(taxon_stats.items()):
                # Identify sample taxa by prefix (e.g. "B1_ivar", "B1_round_3")
                is_sample = taxon.startswith(f"{sample_name}_")
                if not is_sample:
                    # Aggregate GenBank backbone stats
                    genbank_agg["terminal_changes"] += ts["terminal_changes"]
                    genbank_agg["convergent"] += ts["convergent"]
                    genbank_agg["autapomorphic"] += ts["autapomorphic"]
                    genbank_agg["near_convergent"] += ts["near_convergent"]
                    genbank_agg["distant_convergent"] += ts["distant_convergent"]
                    if ts["convergent"] > 0:
                        genbank_agg["conv_dist_sum"] += ts["mean_conv_dist"] * ts["convergent"]
                        genbank_agg["conv_dist_n"] += ts["convergent"]
                    genbank_agg["count"] += 1
                    continue

                # Parse method from taxon name (e.g. "B1_ivar" → method="ivar",
                # "B1_round_3" → method="round_3", "B1_reference" → method="reference")
                method = taxon[len(sample_name) + 1:]  # strip "SAMPLE_" prefix

                tc = ts["terminal_changes"]
                all_taxon_rows.append({
                    "line": line, "segment": segment,
                    "taxon": taxon, "sample": sample_name, "method": method,
                    "terminal_changes": tc,
                    "convergent_changes": ts["convergent"],
                    "near_convergent": ts["near_convergent"],
                    "distant_convergent": ts["distant_convergent"],
                    "autapomorphic_changes": ts["autapomorphic"],
                    "total_variable_sites": ts["total_variable"],
                    "homoplasy_fraction": f"{ts['homoplasy_frac']:.4f}",
                    "mean_convergent_distance": f"{ts['mean_conv_dist']:.4f}",
                })

            # GenBank backbone aggregate row for comparison
            gc = genbank_agg["count"]
            if gc > 0:
                gtc = genbank_agg["terminal_changes"]
                gconv = genbank_agg["convergent"]
                gauto = genbank_agg["autapomorphic"]
                gnear = genbank_agg["near_convergent"]
                gdist = genbank_agg["distant_convergent"]
                gfrac = gconv / gtc if gtc > 0 else 0.0
                gmean_d = (genbank_agg["conv_dist_sum"] / genbank_agg["conv_dist_n"]
                           if genbank_agg["conv_dist_n"] > 0 else 0.0)
                all_taxon_rows.append({
                    "line": line, "segment": segment,
                    "taxon": f"backbone_mean (n={gc})",
                    "sample": sample_name, "method": "backbone",
                    "terminal_changes": round(gtc / gc, 1),
                    "convergent_changes": round(gconv / gc, 1),
                    "near_convergent": round(gnear / gc, 1),
                    "distant_convergent": round(gdist / gc, 1),
                    "autapomorphic_changes": round(gauto / gc, 1),
                    "total_variable_sites": taxon_stats[
                        next((k for k in taxon_stats if not k.startswith(f"{sample_name}_")), None)
                    ]["total_variable"] if any(not k.startswith(f"{sample_name}_") for k in taxon_stats) else 0,
                    "homoplasy_fraction": f"{gfrac:.4f}",
                    "mean_convergent_distance": f"{gmean_d:.4f}",
                })

        # ── Write global summary ──
        hdr = [
            "line", "segment", "sample", "num_sequences", "alignment_length",
            "constant_sites", "constant_sites_pct",
            "variable_sites", "parsimony_informative",
            "parsimony_informative_pct", "singleton_sites",
            "best_model", "log_likelihood", "tree_length",
            "parsimony_score", "consistency_index", "retention_index",
            "homoplasy_index",
        ]

        os.makedirs(os.path.dirname(str(output.summary)), exist_ok=True)
        with open(str(output.summary), "w") as f:
            f.write("\t".join(hdr) + "\n")
            for line, segment, sample_name, info in summary_rows:
                f.write("\t".join([
                    line, segment, sample_name,
                    info["num_sequences"], info["alignment_length"],
                    info["constant_sites"], info["constant_sites_pct"],
                    info["variable_sites"], info["parsimony_informative"],
                    info["parsimony_informative_pct"], info["singleton_sites"],
                    info["best_model"], info["log_likelihood"],
                    info["tree_length"], info["parsimony_score"],
                    info["ci"], info["ri"], info["hi"],
                ]) + "\n")

        # ── Write per-taxon homoplasy ──
        taxon_hdr = [
            "line", "segment", "taxon", "sample", "method",
            "terminal_changes", "convergent_changes",
            "near_convergent", "distant_convergent",
            "autapomorphic_changes", "total_variable_sites",
            "homoplasy_fraction", "mean_convergent_distance",
        ]
        with open(str(output.taxon_report), "w") as f:
            f.write("\t".join(taxon_hdr) + "\n")
            for r in all_taxon_rows:
                f.write("\t".join([
                    str(r[c]) for c in taxon_hdr
                ]) + "\n")

        # ── Analyse individual consensus placements ──
        individual_rows = []

        for placement_tsv in input.placement_files:
            with open(str(placement_tsv)) as f:
                header_line = f.readline()
                if header_line.startswith("#"):
                    continue  # skip empty/no-backbone files
                for row_line in f:
                    fields = row_line.strip().split("\t")
                    if len(fields) < 8:
                        continue
                    p_sample, p_line, p_seg, p_method, p_taxon = fields[:5]
                    p_tree, p_aln, p_iq = fields[5], fields[6], fields[7]

                    if not os.path.exists(p_iq) or not os.path.exists(p_tree):
                        continue

                    # Parse IQ-TREE report for this individual placement
                    indiv_info = parse_iqtree(p_iq)

                    # Compute Fitch parsimony for this individual placement
                    indiv_global, indiv_taxon = analyse_tree(p_tree, p_aln)

                    # Extract the single sample taxon's stats
                    ts = indiv_taxon.get(p_taxon, {})
                    tc = ts.get("terminal_changes", 0)
                    conv = ts.get("convergent", 0)
                    near = ts.get("near_convergent", 0)
                    dist_conv = ts.get("distant_convergent", 0)
                    auto = ts.get("autapomorphic", 0)
                    total_var = ts.get("total_variable", 0)
                    hfrac = ts.get("homoplasy_frac", 0.0)
                    mean_cd = ts.get("mean_conv_dist", 0.0)

                    individual_rows.append({
                        "sample": p_sample,
                        "line": p_line,
                        "segment": p_seg,
                        "method": p_method,
                        "taxon": p_taxon,
                        "num_sequences": indiv_info.get("num_sequences", "NA"),
                        "alignment_length": indiv_info.get("alignment_length", "NA"),
                        "tree_length": indiv_info.get("tree_length", "NA"),
                        "log_likelihood": indiv_info.get("log_likelihood", "NA"),
                        "parsimony_score": indiv_global.get("parsimony_score", "NA"),
                        "consistency_index": indiv_global.get("ci", "NA"),
                        "retention_index": indiv_global.get("ri", "NA"),
                        "homoplasy_index": indiv_global.get("hi", "NA"),
                        "terminal_changes": tc,
                        "convergent_changes": conv,
                        "near_convergent": near,
                        "distant_convergent": dist_conv,
                        "autapomorphic_changes": auto,
                        "total_variable_sites": total_var,
                        "homoplasy_fraction": f"{hfrac:.4f}",
                        "mean_convergent_distance": f"{mean_cd:.4f}",
                    })

        # ── Write individual placement homoplasy report ──
        indiv_hdr = [
            "sample", "line", "segment", "method", "taxon",
            "num_sequences", "alignment_length", "tree_length", "log_likelihood",
            "parsimony_score", "consistency_index", "retention_index", "homoplasy_index",
            "terminal_changes", "convergent_changes",
            "near_convergent", "distant_convergent",
            "autapomorphic_changes", "total_variable_sites",
            "homoplasy_fraction", "mean_convergent_distance",
        ]
        with open(str(output.individual_report), "w") as f:
            f.write("\t".join(indiv_hdr) + "\n")
            for r in individual_rows:
                f.write("\t".join([str(r[c]) for c in indiv_hdr]) + "\n")

        print(f"IQ-TREE summary written to {output.summary}")
        print(f"Per-taxon homoplasy written to {output.taxon_report}")
        print(f"Individual placement homoplasy written to {output.individual_report}")


rule select_best_consensus:
    """Rank consensus sequences by phylogenetic homoplasy, select best, and build hybrid."""
    input:
        placements = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/iqtree/{segment}/consensus_placements.tsv",
            segment=SEGMENTS,
        ),
    output:
        best_fasta   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/best_consensus.fasta",
        hybrid_fasta = SAMPLE_DIR + "/{sample}/{line}/best_consensus/hybrid_consensus.fasta",
        ranking      = SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
        site_report  = SAMPLE_DIR + "/{sample}/{line}/best_consensus/site_homoplasies.tsv",
    run:
        import os

        sample = wildcards.sample
        line = wildcards.line
        out_dir = os.path.dirname(str(output.best_fasta))
        os.makedirs(out_dir, exist_ok=True)

        # ── Parse all placement manifests ──
        # segment -> list of (method, taxon_name, tree_path, aln_path, iq_path)
        seg_placements = {}
        for tsv_path in input.placements:
            with open(str(tsv_path)) as f:
                header_line = f.readline()
                if header_line.startswith("#"):
                    continue
                for row_line in f:
                    fields = row_line.strip().split("\t")
                    if len(fields) < 8:
                        continue
                    _, _, seg, method, taxon = fields[:5]
                    tree_p, aln_p, iq_p = fields[5], fields[6], fields[7]
                    if not os.path.exists(tree_p) or not os.path.exists(aln_p):
                        continue
                    seg_placements.setdefault(seg, []).append(
                        (method, taxon, tree_p, aln_p, iq_p)
                    )

        # ── Part A: site-level Fitch analysis ──
        # site_records: list of dicts for site_homoplasies.tsv
        # method_counts: {(segment, method) -> {convergent, autapomorphic, total_changes, variable}}
        site_records = []
        method_counts = {}

        for seg, entries in sorted(seg_placements.items()):
            for method, taxon, tree_path, aln_path, iq_path in entries:
                key = (seg, method)
                counts = {"convergent": 0, "autapomorphic": 0,
                          "total_changes": 0, "variable": 0}

                with open(tree_path) as f:
                    tree_str = f.read().strip()
                if (tree_str.startswith("(no_tree)") or tree_str.startswith("# No")
                        or not tree_str):
                    method_counts[key] = counts
                    continue

                tree = Phylo.read(StringIO(tree_str), "newick")
                aln = read_fasta(aln_path)
                if not aln:
                    method_counts[key] = counts
                    continue

                terminals = [t.name for t in tree.get_terminals()]
                parent_map = build_parent_map(tree)
                aln_len = len(next(iter(aln.values())))

                # Precompute patristic distances
                patristic = {}
                for i, n1 in enumerate(terminals):
                    for n2 in terminals[i + 1:]:
                        d = tree.distance(n1, n2)
                        patristic[(n1, n2)] = d
                        patristic[(n2, n1)] = d
                all_pat_dists = sorted(patristic.values())
                median_pat = (all_pat_dists[len(all_pat_dists) // 2]
                              if all_pat_dists else 0.0)

                # Find the sample leaf
                sample_leaf = None
                for leaf in tree.get_terminals():
                    if leaf.name == taxon:
                        sample_leaf = leaf
                        break
                if sample_leaf is None:
                    method_counts[key] = counts
                    continue

                parent_clade = parent_map.get(id(sample_leaf))

                for pos in range(aln_len):
                    states = {}
                    for name in terminals:
                        if name in aln:
                            s = aln[name][pos]
                            if s not in ("-", "N", "?"):
                                states[name] = s

                    unique = set(states.values())
                    if len(unique) <= 1:
                        continue
                    counts["variable"] += 1

                    if taxon not in states:
                        continue

                    # Fitch bottom-up + top-down
                    node_cache = {}
                    fitch_bottom_up(tree.root, states, node_cache)
                    assigned = {}
                    fitch_top_down(tree.root, node_cache, assigned)

                    leaf_state = states[taxon]
                    parent_state = assigned.get(id(parent_clade)) if parent_clade else None

                    if leaf_state == parent_state or parent_state is None:
                        continue

                    counts["total_changes"] += 1

                    # Classify the change
                    others = [n for n, s in states.items()
                              if n != taxon and s == leaf_state]
                    if others:
                        counts["convergent"] += 1
                        classification = "convergent"
                        n_partners = len(others)
                    else:
                        counts["autapomorphic"] += 1
                        classification = "autapomorphic"
                        n_partners = 0

                    site_records.append({
                        "sample": sample, "line": line, "segment": seg,
                        "method": method, "position": pos + 1,
                        "sample_state": leaf_state,
                        "ancestral_state": parent_state,
                        "classification": classification,
                        "convergent_partners": n_partners,
                    })

                method_counts[key] = counts

        # ── Write site_homoplasies.tsv ──
        site_hdr = [
            "sample", "line", "segment", "method", "position",
            "sample_state", "ancestral_state", "classification",
            "convergent_partners",
        ]
        with open(str(output.site_report), "w") as f:
            f.write("\t".join(site_hdr) + "\n")
            for r in site_records:
                f.write("\t".join([str(r[c]) for c in site_hdr]) + "\n")

        # ── Part B: ranking ──
        # Group by segment, sort by (convergent ASC, autapomorphic ASC, total_changes ASC)
        ranking_rows = []
        best_per_seg = {}  # segment -> best method

        for seg in SEGMENTS:
            seg_methods = [(m, c) for (s, m), c in method_counts.items() if s == seg]
            if not seg_methods:
                continue
            seg_methods.sort(key=lambda x: (
                x[1]["convergent"], x[1]["autapomorphic"], x[1]["total_changes"]
            ))
            for rank, (method, c) in enumerate(seg_methods, 1):
                ranking_rows.append({
                    "sample": sample, "line": line, "segment": seg,
                    "rank": rank, "method": method,
                    "homoplastic_sites": c["convergent"],
                    "autapomorphic_sites": c["autapomorphic"],
                    "total_changes": c["total_changes"],
                    "total_variable_sites": c["variable"],
                })
                if rank == 1:
                    best_per_seg[seg] = method

        # Write consensus_ranking.tsv
        rank_hdr = [
            "sample", "line", "segment", "rank", "method",
            "homoplastic_sites", "autapomorphic_sites",
            "total_changes", "total_variable_sites",
        ]
        with open(str(output.ranking), "w") as f:
            f.write("\t".join(rank_hdr) + "\n")
            for r in ranking_rows:
                f.write("\t".join([str(r[c]) for c in rank_hdr]) + "\n")

        # ── Part C: best consensus selection ──
        # For each segment, get the raw sequence from the best method's query.fasta
        best_seqs = {}
        for seg, entries in seg_placements.items():
            best_method = best_per_seg.get(seg)
            if not best_method:
                continue
            for method, taxon, tree_p, aln_p, iq_p in entries:
                if method == best_method:
                    # Read the raw query sequence
                    query_fasta = os.path.join(os.path.dirname(tree_p), "query.fasta")
                    if os.path.exists(query_fasta):
                        seqs = read_fasta(query_fasta)
                        if seqs:
                            best_seqs[seg] = (taxon, next(iter(seqs.values())))
                    break

        with open(str(output.best_fasta), "w") as f:
            for seg in SEGMENTS:
                if seg in best_seqs:
                    taxon, seq = best_seqs[seg]
                    f.write(f">{seg}_{sample} best_method={best_per_seg[seg]}\n{seq}\n")

        print(f"[INFO] {line}/{sample}: best consensus methods: "
              + ", ".join(f"{s}={best_per_seg.get(s, 'NA')}" for s in SEGMENTS))

        # ── Part D: hybrid consensus ──
        # For each segment, collect aligned sample sequences from all methods,
        # then build a phylogenetically-informed majority-rule consensus.
        # Positions flagged as homoplastic for a method are deprioritised.
        homoplastic_positions = {}  # (segment, method) -> set of 0-based positions
        for r in site_records:
            if r["classification"] == "convergent":
                key = (r["segment"], r["method"])
                homoplastic_positions.setdefault(key, set()).add(r["position"] - 1)

        hybrid_seqs = {}
        for seg, entries in sorted(seg_placements.items()):
            # Collect aligned sample sequences
            method_aligned = {}  # method -> aligned sequence string
            for method, taxon, tree_p, aln_p, iq_p in entries:
                aln = read_fasta(aln_p)
                if taxon in aln:
                    method_aligned[method] = aln[taxon]

            if not method_aligned:
                continue

            aln_len = len(next(iter(method_aligned.values())))
            hybrid_aligned = []

            for pos in range(aln_len):
                bases = {}  # method -> base
                for method, seq in method_aligned.items():
                    b = seq[pos]
                    if b not in ("-", "N", "?"):
                        bases[method] = b

                if not bases:
                    # All methods have gap/N at this position
                    # Use most common character (including N/gap)
                    all_chars = [method_aligned[m][pos] for m in method_aligned]
                    hybrid_aligned.append(Counter(all_chars).most_common(1)[0][0])
                    continue

                unique_bases = set(bases.values())
                if len(unique_bases) == 1:
                    # All agree
                    hybrid_aligned.append(next(iter(unique_bases)))
                    continue

                # Disagreement: prefer bases from non-homoplastic methods
                non_homo_bases = []
                for method, base in bases.items():
                    if pos not in homoplastic_positions.get((seg, method), set()):
                        non_homo_bases.append(base)

                if non_homo_bases:
                    # Majority rule among non-homoplastic methods
                    hybrid_aligned.append(Counter(non_homo_bases).most_common(1)[0][0])
                else:
                    # All are homoplastic at this position; plain majority rule
                    hybrid_aligned.append(
                        Counter(bases.values()).most_common(1)[0][0]
                    )

            # Strip gaps to get raw hybrid sequence
            hybrid_raw = "".join(b for b in hybrid_aligned if b != "-")
            hybrid_seqs[seg] = hybrid_raw

        with open(str(output.hybrid_fasta), "w") as f:
            for seg in SEGMENTS:
                if seg in hybrid_seqs:
                    f.write(f">{seg}_{sample} hybrid_consensus\n{hybrid_seqs[seg]}\n")

        print(f"[INFO] {line}/{sample}: hybrid consensus written → {output.hybrid_fasta}")


rule sample_report:
    """Generate a comprehensive HTML report consolidating all sample results."""
    input:
        summary      = SAMPLE_DIR + "/{sample}/{line}/summary.tsv",
        conv_stats   = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        final_bam    = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        homoplasy    = SAMPLE_DIR + "/{sample}/{line}/homoplasy/homoplasy_report.tsv",
        pairwise     = SAMPLE_DIR + "/{sample}/{line}/homoplasy/pairwise_distances.tsv",
        ranking      = SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
        site_homo    = SAMPLE_DIR + "/{sample}/{line}/best_consensus/site_homoplasies.tsv",
        best_fasta   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/best_consensus.fasta",
        validation   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/validation.tsv",
        read_support = SAMPLE_DIR + "/{sample}/{line}/best_consensus/read_support.tsv",
        placements   = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/iqtree/{segment}/consensus_placements.tsv",
            segment=SEGMENTS,
        ),
    output:
        report = SAMPLE_DIR + "/{sample}/{line}/report.html",
    params:
        segments = SEGMENTS,
    shell:
        """
        python {workflow.basedir}/scripts/sample_report.py \
            --sample {wildcards.sample} --line {wildcards.line} \
            --summary {input.summary} \
            --convergence {input.conv_stats} \
            --final-bam {input.final_bam} \
            --homoplasy {input.homoplasy} \
            --pairwise {input.pairwise} \
            --ranking {input.ranking} \
            --site-homoplasy {input.site_homo} \
            --best-fasta {input.best_fasta} \
            --validation {input.validation} \
            --read-support {input.read_support} \
            --placements {input.placements} \
            --segments {params.segments} \
            --output {output.report}
        """


rule sample_readme:
    """Generate a README describing the directory structure for a sample/line."""
    input:
        report = SAMPLE_DIR + "/{sample}/{line}/report.html",
    output:
        readme = SAMPLE_DIR + "/{sample}/{line}/README.md",
    run:
        assemblers = ASSEMBLERS
        content = f"""\
# {wildcards.sample} / {wildcards.line}

Output directory for sample **{wildcards.sample}**, analysis line **{wildcards.line}**.

## Directory structure

### Top-level files

| File | Description |
|------|-------------|
| `report.html` | Comprehensive HTML report with coverage plots, convergence graphs, consensus quality, homoplasy analysis, and phylogenetic placement |
| `summary.tsv` | Per-segment coverage statistics (mapped reads, mean depth, completeness) across all pipeline stages |
| `README.md` | This file |

### `iterative/` — Iterative reference mapping

Reads are mapped against the auto-selected reference and the consensus is
refined over multiple rounds until convergence.

| File | Description |
|------|-------------|
| `convergence_stats.tsv` | Per-round mapped reads and genome completeness |
| `final_consensus.fasta` | Consensus sequence after the final iteration |
| `analysis.bam` | BAM from the final iteration (used downstream) |
| `round_N/mapped.bam` | Alignments for iteration N |
| `round_N/reference.fasta` | Reference used in iteration N |
| `round_N/consensus.fasta` | Consensus called after iteration N |
| `round_N/variants.vcf.gz` | Filtered variant calls for iteration N |
| `round_N/variants.raw.vcf.gz` | Unfiltered variant calls |
| `round_N/variants.inframe.vcf.gz` | In-frame filtered variants |

### `ivar/` — iVar consensus calling

Reads are re-mapped to the iterative consensus and iVar generates a new
consensus using a frequency-based approach.

| File | Description |
|------|-------------|
| `mapped.bam` | Alignments to the iterative consensus |
| `reference.fasta` | Reference used (iterative consensus) |
| `final_consensus.fasta` | iVar consensus sequence |

### `final/` — Final mapping

Reads are mapped one last time to the iVar consensus to produce the
definitive alignment used for depth plots and read-support analysis.

| File | Description |
|------|-------------|
| `mapped.bam` | Final alignments |
| `reference.fasta` | iVar consensus used as reference |

### `assembly_reads/` — Mapped read extraction

| File | Description |
|------|-------------|
| `mapped_R1.fastq.gz` | Forward reads that mapped during iterative mapping |
| `mapped_R2.fastq.gz` | Reverse reads that mapped during iterative mapping |

### De novo assemblers

Each assembler has two directories: the raw assembly output and a
`*_polish/` directory where contigs are aligned to the iVar consensus
and polished.

"""
        for asm in assemblers:
            content += f"""\
#### `{asm}/` — {asm} raw assembly output

Raw contigs/scaffolds produced by {asm}. Contains assembler-specific
files (contigs, scaffolds, logs, graphs).

#### `{asm}_polish/` — Assembly-guided polishing

| File | Description |
|------|-------------|
| `{asm}_guided_polished.fasta` | Polished consensus: iVar consensus corrected by assembly evidence |
| `assembly_vs_consensus.bam` | Assembly contigs aligned to the iVar consensus |
| `variants.vcf.gz` | Variants between assembly and iVar consensus |

"""
        content += """\
### `final_alignment/` — Multi-method alignment

| File | Description |
|------|-------------|
| `aligned.fasta` | MAFFT alignment of all consensus methods (iVar, iterative rounds, assembler-polished, reference) for the sample |

### `homoplasy/` — Homoplasy analysis

| File | Description |
|------|-------------|
| `homoplasy_report.tsv` | Delta score (treelikeness) and parsimony metrics per segment |
| `pairwise_distances.tsv` | Pairwise patristic distances between consensus methods |

### `best_consensus/` — Consensus ranking and selection

| File | Description |
|------|-------------|
| `consensus_ranking.tsv` | All consensus methods ranked by homoplastic site count per segment |
| `best_consensus.fasta` | The top-ranked consensus for each segment |
| `hybrid_consensus.fasta` | Hybrid consensus combining best bases from multiple methods |
| `site_homoplasies.tsv` | Per-site homoplasy classification (convergent / autapomorphic) |
| `read_support.tsv` | Allele frequencies from the final BAM at disputed (homoplastic) sites |
| `validation.tsv` | Comparison of best/hybrid consensus vs. original rank-1 method |
| `validation/` | IQ-TREE runs for best and hybrid consensus placement |

### `iqtree/` — Phylogenetic placement

Per-segment subdirectories (`L/`, `M/`, `S/`), each containing:

| File | Description |
|------|-------------|
| `aligned.fasta` | Backbone alignment with GenBank references |
| `iqtree.treefile` | Maximum-likelihood tree |
| `iqtree.iqtree` | IQ-TREE report (model, log-likelihood, tree length) |
| `sample_seqs.fasta` | Sample consensus sequences used as queries |
| `consensus_placements.tsv` | Placement summary for each consensus method |
| `placements/<method>/` | Per-method IQ-TREE placement run |
"""
        with open(str(output.readme), "w") as f:
            f.write(content)

        print(f"[INFO] {wildcards.line}/{wildcards.sample}: README written → {output.readme}")


rule sample_root_readme:
    """Generate a README for the sample root directory."""
    input:
        readmes = expand(
            SAMPLE_DIR + "/{{sample}}/{line}/README.md",
            line=LINES,
        ),
    output:
        readme = SAMPLE_DIR + "/{sample}/README.md",
    run:
        lines_str = ", ".join(f"`{l}`" for l in LINES)
        content = f"""\
# {wildcards.sample}

Output directory for sample **{wildcards.sample}**.

## Directory structure

| Directory | Description |
|-----------|-------------|
| `host_filtered/` | Host-depleted reads (human and bank vole removed via bowtie2) |
| `autoref/` | Automatic reference selection — per-segment best GenBank reference |
| `reference/` | Combined reference FASTA (L + M + S best references concatenated) |
| `dedup/` | Analysis line with PCR duplicate removal (see `dedup/README.md`) |
| `nodedup/` | Analysis line without duplicate removal (see `nodedup/README.md`) |

### `host_filtered/`

| File | Description |
|------|-------------|
| `clean_R1.fastq.gz` | Forward reads after host removal |
| `clean_R2.fastq.gz` | Reverse reads after host removal |
| `human_bowtie2.log` | Bowtie2 log for human genome filtering |
| `bankvole_bowtie2.log` | Bowtie2 log for bank vole genome filtering |

### `autoref/`

Per-segment subdirectories (`L/`, `M/`, `S/`), each containing:

| File | Description |
|------|-------------|
| `best.fasta` | Best-matching GenBank reference for this segment |
| `ref_stats.tsv` | Mapping statistics for all candidate references |
| `mapped.bam` | Alignments used for reference evaluation |

Top-level file:

| File | Description |
|------|-------------|
| `best_refs.fasta` | Combined best references (all segments) |

### `reference/`

| File | Description |
|------|-------------|
| `ref.fasta` | Combined reference FASTA used as starting point for iterative mapping |
| `ref.fasta.fai` | samtools faidx index |

### Analysis lines ({lines_str})

Each analysis line runs the full pipeline independently. The `dedup` line
removes PCR duplicates before consensus calling; the `nodedup` line retains
all reads. See the README.md inside each line directory for a detailed
description of all subdirectories and files.
"""
        with open(str(output.readme), "w") as f:
            f.write(content)

        print(f"[INFO] {wildcards.sample}: root README written → {output.readme}")


rule validate_consensus:
    """Place best and hybrid consensus back on backbone tree and measure homoplasy."""
    input:
        best_fasta   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/best_consensus.fasta",
        hybrid_fasta = SAMPLE_DIR + "/{sample}/{line}/best_consensus/hybrid_consensus.fasta",
        ranking      = SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
        backbone_tree = expand(
            COMMON_DIR + "/iqtree/backbone/{segment}/backbone.treefile", segment=SEGMENTS,
        ),
        backbone_aln = expand(
            COMMON_DIR + "/iqtree/backbone/{segment}/aligned.fasta", segment=SEGMENTS,
        ),
        backbone_iq = expand(
            COMMON_DIR + "/iqtree/backbone/{segment}/backbone.iqtree", segment=SEGMENTS,
        ),
    output:
        validation = SAMPLE_DIR + "/{sample}/{line}/best_consensus/validation.tsv",
    threads: THREADS
    run:
        import subprocess, os, re

        sample = wildcards.sample
        line = wildcards.line
        out_dir = os.path.dirname(str(output.validation))
        val_dir = os.path.join(out_dir, "validation")
        os.makedirs(val_dir, exist_ok=True)

        # Build segment -> backbone paths mapping
        backbone = {}
        for i, seg in enumerate(SEGMENTS):
            tree_p = str(input.backbone_tree[i])
            aln_p = str(input.backbone_aln[i])
            iq_p = str(input.backbone_iq[i])
            with open(tree_p) as f:
                tree_str = f.read().strip()
            if tree_str.startswith("(no_tree)") or not tree_str:
                continue
            # Parse best model
            best_model = "MFP"
            with open(iq_p) as f:
                for bline in f:
                    m = re.search(r"Best-fit model according to BIC:\s+(.+)", bline)
                    if m:
                        best_model = m.group(1).strip()
                        break
            backbone[seg] = {"tree": tree_p, "aln": aln_p, "model": best_model}

        # Read ranking to get original best-method metrics for comparison
        rank1_metrics = {}  # segment -> {convergent, autapomorphic, total_changes}
        with open(str(input.ranking)) as f:
            header = f.readline()
            for row_line in f:
                fields = row_line.strip().split("\t")
                if len(fields) >= 9 and fields[3] == "1":
                    seg = fields[2]
                    rank1_metrics[seg] = {
                        "method": fields[4],
                        "convergent": int(fields[5]),
                        "autapomorphic": int(fields[6]),
                        "total_changes": int(fields[7]),
                        "variable": int(fields[8]),
                    }

        def extract_segment_seq(fasta_path, segment, sample_name):
            """Extract a segment's sequence from a multi-segment FASTA."""
            seqs = read_fasta(fasta_path)
            target = f"{segment}_{sample_name}"
            if target in seqs:
                return seqs[target]
            for name, seq in seqs.items():
                if segment in name:
                    return seq
            return None

        def place_and_analyse(taxon_name, seq, seg, sub_dir):
            """Place one sequence on backbone and compute site-level Fitch metrics."""
            bb = backbone.get(seg)
            if not bb:
                return None
            os.makedirs(sub_dir, exist_ok=True)

            # Write query
            query_path = os.path.join(sub_dir, "query.fasta")
            with open(query_path, "w") as fh:
                fh.write(f">{taxon_name}\n{seq}\n")

            # MAFFT --add
            res = subprocess.run(
                ["mafft", "--add", query_path, "--keeplength",
                 "--thread", str(threads), bb["aln"]],
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                return None
            aln_path = os.path.join(sub_dir, "aligned.fasta")
            with open(aln_path, "w") as fh:
                fh.write(res.stdout)

            # IQ-TREE
            prefix = os.path.join(sub_dir, "iqtree")
            res = subprocess.run(
                ["iqtree", "-s", aln_path, "-g", bb["tree"],
                 "-m", bb["model"], "-T", str(threads),
                 "--prefix", prefix, "-redo"],
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                return None

            tree_path = prefix + ".treefile"

            # Site-level Fitch analysis
            with open(tree_path) as f:
                tree_str = f.read().strip()
            if not tree_str or tree_str.startswith("(no_tree)"):
                return None

            tree = Phylo.read(StringIO(tree_str), "newick")
            aln = read_fasta(aln_path)
            if not aln:
                return None

            terminals = [t.name for t in tree.get_terminals()]
            parent_map = build_parent_map(tree)

            sample_leaf = None
            for leaf in tree.get_terminals():
                if leaf.name == taxon_name:
                    sample_leaf = leaf
                    break
            if not sample_leaf:
                return None

            parent_clade = parent_map.get(id(sample_leaf))
            aln_len = len(next(iter(aln.values())))

            counts = {"convergent": 0, "autapomorphic": 0,
                      "total_changes": 0, "variable": 0}

            for pos in range(aln_len):
                states = {}
                for name in terminals:
                    if name in aln:
                        s = aln[name][pos]
                        if s not in ("-", "N", "?"):
                            states[name] = s
                if len(set(states.values())) <= 1:
                    continue
                counts["variable"] += 1
                if taxon_name not in states:
                    continue

                node_cache = {}
                fitch_bottom_up(tree.root, states, node_cache)
                assigned = {}
                fitch_top_down(tree.root, node_cache, assigned)

                leaf_state = states[taxon_name]
                parent_state = assigned.get(id(parent_clade)) if parent_clade else None
                if leaf_state == parent_state or parent_state is None:
                    continue
                counts["total_changes"] += 1
                others = [n for n, s in states.items()
                          if n != taxon_name and s == leaf_state]
                if others:
                    counts["convergent"] += 1
                else:
                    counts["autapomorphic"] += 1

            return counts

        # ── Run validation for best and hybrid ──
        validation_rows = []

        for seg in SEGMENTS:
            for cons_type, fasta_path in [("best", str(input.best_fasta)),
                                          ("hybrid", str(input.hybrid_fasta))]:
                seq = extract_segment_seq(fasta_path, seg, sample)
                if not seq or not seq.replace("N", ""):
                    continue

                taxon = f"{sample}_{cons_type}"
                sub_dir = os.path.join(val_dir, f"{cons_type}_{seg}")
                counts = place_and_analyse(taxon, seq, seg, sub_dir)

                if counts:
                    validation_rows.append({
                        "sample": sample, "line": line, "segment": seg,
                        "consensus_type": cons_type, "taxon": taxon,
                        "homoplastic_sites": counts["convergent"],
                        "autapomorphic_sites": counts["autapomorphic"],
                        "total_changes": counts["total_changes"],
                        "total_variable_sites": counts["variable"],
                    })

            # Add rank-1 original method for comparison
            r1 = rank1_metrics.get(seg)
            if r1:
                validation_rows.append({
                    "sample": sample, "line": line, "segment": seg,
                    "consensus_type": f"original_best ({r1['method']})",
                    "taxon": f"{sample}_{r1['method']}",
                    "homoplastic_sites": r1["convergent"],
                    "autapomorphic_sites": r1["autapomorphic"],
                    "total_changes": r1["total_changes"],
                    "total_variable_sites": r1["variable"],
                })

        # Write validation.tsv
        val_hdr = [
            "sample", "line", "segment", "consensus_type", "taxon",
            "homoplastic_sites", "autapomorphic_sites",
            "total_changes", "total_variable_sites",
        ]
        with open(str(output.validation), "w") as f:
            f.write("\t".join(val_hdr) + "\n")
            for r in validation_rows:
                f.write("\t".join([str(r[c]) for c in val_hdr]) + "\n")

        print(f"[INFO] {line}/{sample}: validation completed → {output.validation}")


rule read_support:
    """Check read-level allele frequencies at sites flagged as homoplastic."""
    input:
        site_homo  = SAMPLE_DIR + "/{sample}/{line}/best_consensus/site_homoplasies.tsv",
        final_bam  = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        ivar_cons  = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        placements = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/iqtree/{segment}/consensus_placements.tsv",
            segment=SEGMENTS,
        ),
    output:
        support = SAMPLE_DIR + "/{sample}/{line}/best_consensus/read_support.tsv",
    run:
        import subprocess, os, re
        from collections import defaultdict

        sample = wildcards.sample
        line = wildcards.line

        # ── Read site_homoplasies.tsv ──
        sites = []
        with open(str(input.site_homo)) as f:
            header = f.readline()
            if header.startswith("#"):
                # No data
                with open(str(output.support), "w") as out:
                    out.write("# No homoplastic sites to validate\n")
                return
            for row_line in f:
                fields = row_line.strip().split("\t")
                if len(fields) >= 9:
                    sites.append({
                        "sample": fields[0], "line": fields[1],
                        "segment": fields[2], "method": fields[3],
                        "position": int(fields[4]),  # 1-based alignment position
                        "sample_state": fields[5],
                        "ancestral_state": fields[6],
                        "classification": fields[7],
                        "convergent_partners": fields[8],
                    })

        if not sites:
            with open(str(output.support), "w") as out:
                out.write("# No homoplastic sites to validate\n")
            return

        # ── Build alignment position → ungapped position mapping ──
        # For each (segment, method), read the aligned sequence and build
        # aln_pos (1-based) → ungapped_pos (1-based) mapping
        placement_dirs = {}  # (segment, method) -> base directory
        for tsv_path in input.placements:
            with open(str(tsv_path)) as f:
                hdr = f.readline()
                if hdr.startswith("#"):
                    continue
                for row_line in f:
                    fields = row_line.strip().split("\t")
                    if len(fields) >= 8:
                        seg, method, taxon = fields[2], fields[3], fields[4]
                        aln_p = fields[6]
                        placement_dirs[(seg, method)] = {
                            "aln": aln_p, "taxon": taxon,
                        }

        aln_to_ungapped = {}  # (segment, method) -> {aln_pos: ungapped_pos}
        for (seg, method), info in placement_dirs.items():
            aln_path = info["aln"]
            taxon = info["taxon"]
            if not os.path.exists(aln_path):
                continue
            aln_seqs = read_fasta(aln_path)
            if taxon not in aln_seqs:
                continue
            aligned_seq = aln_seqs[taxon]
            mapping = {}
            ungapped = 0
            for i, base in enumerate(aligned_seq):
                if base != "-":
                    ungapped += 1
                    mapping[i + 1] = ungapped  # 1-based
            aln_to_ungapped[(seg, method)] = mapping

        # ── Get contig names from iVar consensus ──
        ivar_seqs = read_fasta(str(input.ivar_cons))
        seg_to_contig = {}
        for contig_name in ivar_seqs:
            for seg in SEGMENTS:
                if seg in contig_name:
                    seg_to_contig[seg] = contig_name
                    break

        # ── Run samtools mpileup per contig ──
        pileup_data = {}  # contig -> {pos: {A:n, C:n, G:n, T:n, del:n, depth:n}}
        for seg, contig in seg_to_contig.items():
            result = subprocess.run(
                ["samtools", "mpileup", "-a", "-A", "-Q", "0",
                 "-r", contig, str(input.final_bam)],
                capture_output=True, text=True,
            )
            contig_pileup = {}
            for pline in result.stdout.strip().split("\n"):
                if not pline:
                    continue
                parts = pline.split("\t")
                if len(parts) < 5:
                    continue
                pos = int(parts[1])
                ref_base = parts[2].upper()
                depth = int(parts[3])
                bases_str = parts[4].upper()

                counts = {"A": 0, "C": 0, "G": 0, "T": 0, "del": 0}
                i = 0
                while i < len(bases_str):
                    c = bases_str[i]
                    if c in (".", ","):
                        counts[ref_base] = counts.get(ref_base, 0) + 1
                    elif c in ("A", "C", "G", "T"):
                        counts[c] += 1
                    elif c == "*" or c == "#":
                        counts["del"] += 1
                    elif c == "^":
                        i += 1  # skip quality char
                    elif c in ("+", "-"):
                        # Indel: read the length digits then skip that many bases
                        j = i + 1
                        while j < len(bases_str) and bases_str[j].isdigit():
                            j += 1
                        indel_len = int(bases_str[i+1:j]) if j > i + 1 else 0
                        i = j + indel_len - 1
                        if c == "-":
                            counts["del"] += 1
                    i += 1

                counts["depth"] = depth
                counts["ref_base"] = ref_base
                contig_pileup[pos] = counts
            pileup_data[seg] = contig_pileup

        # ── Annotate each site with read support ──
        support_rows = []
        for site in sites:
            seg = site["segment"]
            method = site["method"]
            aln_pos = site["position"]

            # Map alignment position to ungapped position
            mapping = aln_to_ungapped.get((seg, method), {})
            ungapped_pos = mapping.get(aln_pos)
            if ungapped_pos is None:
                continue

            # Look up pileup
            contig_pileup = pileup_data.get(seg, {})
            pile = contig_pileup.get(ungapped_pos, {})
            if not pile:
                continue

            depth = pile.get("depth", 0)
            ref_base = pile.get("ref_base", "N")
            cons_base = site["sample_state"].upper()
            cons_count = pile.get(cons_base, 0)
            cons_freq = cons_count / depth if depth > 0 else 0.0

            if cons_freq > 0.8:
                verdict = "supported"
            elif cons_freq >= 0.5:
                verdict = "weak"
            else:
                verdict = "contradicted"

            support_rows.append({
                "sample": sample, "line": line,
                "segment": seg, "method": method,
                "aln_position": aln_pos,
                "contig_position": ungapped_pos,
                "sample_state": site["sample_state"],
                "ancestral_state": site["ancestral_state"],
                "classification": site["classification"],
                "depth": depth,
                "ref_base": ref_base,
                "A_count": pile.get("A", 0),
                "C_count": pile.get("C", 0),
                "G_count": pile.get("G", 0),
                "T_count": pile.get("T", 0),
                "del_count": pile.get("del", 0),
                "consensus_allele_freq": f"{cons_freq:.3f}",
                "read_support_verdict": verdict,
            })

        # Write output
        sup_hdr = [
            "sample", "line", "segment", "method",
            "aln_position", "contig_position",
            "sample_state", "ancestral_state", "classification",
            "depth", "ref_base",
            "A_count", "C_count", "G_count", "T_count", "del_count",
            "consensus_allele_freq", "read_support_verdict",
        ]
        with open(str(output.support), "w") as f:
            f.write("\t".join(sup_hdr) + "\n")
            for r in support_rows:
                f.write("\t".join([str(r[c]) for c in sup_hdr]) + "\n")

        print(f"[INFO] {line}/{sample}: read support for {len(support_rows)} sites → {output.support}")


rule cross_sample_report:
    """Generate a global HTML overview comparing all samples."""
    input:
        summaries = expand(
            SAMPLE_DIR + "/{sample}/{line}/summary.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
        rankings = expand(
            SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
        validations = expand(
            SAMPLE_DIR + "/{sample}/{line}/best_consensus/validation.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
    output:
        report = COMMON_DIR + "/overview_report.html",
    shell:
        """
        python {workflow.basedir}/scripts/cross_sample_report.py \
            --summaries {input.summaries} \
            --rankings {input.rankings} \
            --validations {input.validations} \
            --output {output.report}
        """
