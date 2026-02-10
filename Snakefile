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


# ── Configuration ──────────────────────────────────────────────

FASTQ_DIR  = config["fastq_dir"]
REF_DIR    = config["ref_dir"]
OUTDIR     = config["outdir"]
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
SEGMENTS           = ["L", "M", "S"]
AUTO_SELECT_REFS   = config.get("auto_select_refs", False)
DENOVO_SEED        = config.get("denovo_seed", False)
EXTEND_REFERENCE   = config.get("extend_reference", False)


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
    return f"{OUTDIR}/{wildcards.sample}/host_filtered/clean_R1.fastq.gz"

def get_r2(wildcards):
    return f"{OUTDIR}/{wildcards.sample}/host_filtered/clean_R2.fastq.gz"

def _analysis_bam(base, line):
    """Return dedup BAM path for dedup line, raw BAM for nodedup."""
    if line == "dedup":
        return f"{base}/mapped.dedup.bam"
    return f"{base}/mapped.bam"

def get_assembly_fasta(wildcards):
    """Return the assembly FASTA path for a given assembler."""
    paths = {
        "trinity":      f"{OUTDIR}/{wildcards.sample}/{wildcards.line}/trinity/assembly/Trinity.fasta",
        "spades":       f"{OUTDIR}/{wildcards.sample}/{wildcards.line}/spades/contigs.fasta",
        "coronaspades": f"{OUTDIR}/{wildcards.sample}/{wildcards.line}/coronaspades/contigs.fasta",
        "megahit":      f"{OUTDIR}/{wildcards.sample}/{wildcards.line}/megahit/final.contigs.fa",
        "iva":          f"{OUTDIR}/{wildcards.sample}/{wildcards.line}/iva/contigs.fasta",
    }
    return paths[wildcards.assembler]

def get_ref_input(wildcards):
    """Return reference FASTA: auto-selected or pre-existing."""
    if AUTO_SELECT_REFS:
        return f"{OUTDIR}/{wildcards.sample}/autoref/best_refs.fasta"
    return f"{REF_DIR}/{wildcards.sample}_best_refs.fasta"

def get_iterative_ref(wildcards):
    """Return starting reference for iterative mapping: seeded or original."""
    if DENOVO_SEED:
        return f"{OUTDIR}/{wildcards.sample}/denovo_seed/seeded_ref.fasta"
    return f"{OUTDIR}/{wildcards.sample}/reference/ref.fasta"

def get_iterative_ref_fai(wildcards):
    """Return .fai for the iterative mapping starting reference."""
    if DENOVO_SEED:
        return f"{OUTDIR}/{wildcards.sample}/denovo_seed/seeded_ref.fasta.fai"
    return f"{OUTDIR}/{wildcards.sample}/reference/ref.fasta.fai"


# ── Wildcard constraints ──────────────────────────────────────

wildcard_constraints:
    sample    = r"[^/]+",
    line      = r"dedup|nodedup",
    assembler = r"trinity|spades|coronaspades|megahit|iva",
    segment   = r"L|M|S",


# ═══════════════════════════════════════════════════════════════
# TARGET
# ═══════════════════════════════════════════════════════════════

rule all:
    input:
        expand(
            "{outdir}/{sample}/{line}/summary.tsv",
            outdir=OUTDIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{outdir}/{sample}/{line}/homoplasy/homoplasy_report.tsv",
            outdir=OUTDIR,
            sample=SAMPLES.keys(),
            line=LINES,
        ),
        expand(
            "{outdir}/iqtree/{line}/{segment}/iqtree.treefile",
            outdir=OUTDIR,
            line=LINES,
            segment=SEGMENTS,
        ),
        OUTDIR + "/iqtree/iqtree_summary.tsv",


# ═══════════════════════════════════════════════════════════════
# HOST REMOVAL
# ═══════════════════════════════════════════════════════════════

rule host_removal:
    """Remove human and bank vole host reads using bowtie2 --very-sensitive."""
    input:
        r1 = get_raw_r1,
        r2 = get_raw_r2,
    output:
        r1 = OUTDIR + "/{sample}/host_filtered/clean_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/host_filtered/clean_R2.fastq.gz",
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
        r1 = OUTDIR + "/{sample}/host_filtered/clean_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/host_filtered/clean_R2.fastq.gz",
        ref = GENBANK_DIR + "/{segment}.full_len.cds.clean.uniq.fa",
    output:
        bam = OUTDIR + "/{sample}/autoref/{segment}/mapped.bam",
        bai = OUTDIR + "/{sample}/autoref/{segment}/mapped.bam.bai",
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
        bam = OUTDIR + "/{sample}/autoref/{segment}/mapped.bam",
        bai = OUTDIR + "/{sample}/autoref/{segment}/mapped.bam.bai",
        ref = GENBANK_DIR + "/{segment}.full_len.cds.clean.uniq.fa",
    output:
        stats = OUTDIR + "/{sample}/autoref/{segment}/ref_stats.tsv",
        best_fa = OUTDIR + "/{sample}/autoref/{segment}/best.fasta",
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
    """Combine best L, M, S references into a single per-sample FASTA."""
    input:
        L = OUTDIR + "/{sample}/autoref/L/best.fasta",
        M = OUTDIR + "/{sample}/autoref/M/best.fasta",
        S = OUTDIR + "/{sample}/autoref/S/best.fasta",
    output:
        ref = OUTDIR + "/{sample}/autoref/best_refs.fasta",
    shell:
        "cat {input.L} {input.M} {input.S} > {output.ref}"


# ═══════════════════════════════════════════════════════════════
# REFERENCE PREPARATION
# ═══════════════════════════════════════════════════════════════

rule rename_ref:
    """Rename reference contigs by size: L (largest), M (middle), S (smallest)."""
    input:
        ref = get_ref_input,
    output:
        ref = OUTDIR + "/{sample}/reference/ref.fasta",
    run:
        seqs = []
        name = None
        parts = []
        with open(str(input.ref)) as f:
            for ln in f:
                ln = ln.rstrip("\n")
                if ln.startswith(">"):
                    if name:
                        seqs.append((name, "".join(parts)))
                    name = ln[1:].split()[0]
                    parts = []
                else:
                    parts.append(ln)
        if name:
            seqs.append((name, "".join(parts)))

        # Sort by length descending → L, M, S
        seqs.sort(key=lambda x: len(x[1]), reverse=True)
        labels = ["L", "M", "S"]

        import os
        os.makedirs(os.path.dirname(str(output.ref)), exist_ok=True)
        with open(str(output.ref), "w") as out:
            for i, (orig_name, seq) in enumerate(seqs):
                if i < len(labels):
                    new_name = f"{labels[i]}_{wildcards.sample}"
                else:
                    new_name = f"{orig_name}_{wildcards.sample}"
                out.write(f">{new_name}\n")
                for j in range(0, len(seq), 80):
                    out.write(seq[j:j+80] + "\n")


rule index_ref:
    """Index the renamed per-sample reference (shared by both lines)."""
    input:
        ref = OUTDIR + "/{sample}/reference/ref.fasta",
    output:
        fai = OUTDIR + "/{sample}/reference/ref.fasta.fai",
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
        contigs = OUTDIR + "/{sample}/denovo_seed/megahit/final.contigs.fa",
    params:
        outdir = OUTDIR + "/{sample}/denovo_seed/megahit",
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
        contigs = OUTDIR + "/{sample}/denovo_seed/megahit/final.contigs.fa",
        ref     = OUTDIR + "/{sample}/reference/ref.fasta",
        fai     = OUTDIR + "/{sample}/reference/ref.fasta.fai",
    output:
        seeded = OUTDIR + "/{sample}/denovo_seed/seeded_ref.fasta",
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
        ref = OUTDIR + "/{sample}/denovo_seed/seeded_ref.fasta",
    output:
        fai = OUTDIR + "/{sample}/denovo_seed/seeded_ref.fasta.fai",
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
        consensus = OUTDIR + "/{sample}/{line}/iterative/final_consensus.fasta",
        bam       = OUTDIR + "/{sample}/{line}/iterative/analysis.bam",
        bai       = OUTDIR + "/{sample}/{line}/iterative/analysis.bam.bai",
        stats     = OUTDIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
    params:
        sample          = "{sample}",
        line            = "{line}",
        aligner         = ALIGNER,
        max_iterations  = MAX_ITERATIONS,
        min_iterations  = MIN_ITERATIONS,
        enforce_inframe = ENFORCE_INFRAME,
        extend_ref      = EXTEND_REFERENCE,
        outdir          = OUTDIR + "/{sample}/{line}/iterative",
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
        bam = OUTDIR + "/{sample}/{line}/iterative/analysis.bam",
    output:
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
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
    """Run Trinity de novo assembly via Docker."""
    input:
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = OUTDIR + "/{sample}/{line}/trinity/assembly/Trinity.fasta",
    params:
        memory = TRINITY_MEMORY,
    threads: THREADS
    shell:
        """
        set -euo pipefail
        READS_DIR=$(realpath $(dirname {input.r1}))
        LINE_DIR=$(dirname "$READS_DIR")
        TRINITY_DIR="$LINE_DIR/trinity"
        TRINITY_OUT="$TRINITY_DIR/trinity_assembly"
        R1=$(realpath {input.r1})
        R2=$(realpath {input.r2})

        docker run --rm -v "$LINE_DIR":"$LINE_DIR" trinityrnaseq/trinityrnaseq rm -rf "$TRINITY_OUT"

        docker run --rm \
            -v "$LINE_DIR":"$LINE_DIR" \
            -w "$TRINITY_DIR" \
            trinityrnaseq/trinityrnaseq Trinity \
            --seqType fq \
            --left "$R1" \
            --right "$R2" \
            --max_memory {params.memory}G \
            --CPU {threads} \
            --output "$TRINITY_OUT"

        docker run --rm -v "$LINE_DIR":"$LINE_DIR" trinityrnaseq/trinityrnaseq chown -R $(id -u):$(id -g) "$TRINITY_OUT" "$TRINITY_OUT.Trinity.fasta"

        mkdir -p "$(dirname {output.fasta})"
        cp "$TRINITY_OUT.Trinity.fasta" {output.fasta}
        """


rule spades:
    """Run SPAdes in rnaviral mode."""
    input:
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = OUTDIR + "/{sample}/{line}/spades/contigs.fasta",
    params:
        outdir = OUTDIR + "/{sample}/{line}/spades",
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
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = OUTDIR + "/{sample}/{line}/megahit/final.contigs.fa",
    params:
        outdir = OUTDIR + "/{sample}/{line}/megahit",
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
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = OUTDIR + "/{sample}/{line}/coronaspades/contigs.fasta",
    params:
        outdir = OUTDIR + "/{sample}/{line}/coronaspades",
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
        """


rule iva:
    """Run IVA (Iterative Virus Assembler)."""
    input:
        r1 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R1.fastq.gz",
        r2 = OUTDIR + "/{sample}/{line}/assembly_reads/mapped_R2.fastq.gz",
    output:
        fasta = OUTDIR + "/{sample}/{line}/iva/contigs.fasta",
    params:
        outdir = OUTDIR + "/{sample}/{line}/iva",
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
        ref = OUTDIR + "/{sample}/{line}/iterative/final_consensus.fasta",
        r1  = get_r1,
        r2  = get_r2,
    output:
        bam = OUTDIR + "/{sample}/{line}/ivar/mapped.bam",
        bai = OUTDIR + "/{sample}/{line}/ivar/mapped.bam.bai",
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
            f"{OUTDIR}/{wc.sample}/{wc.line}/ivar", wc.line
        ),
        ref = OUTDIR + "/{sample}/{line}/iterative/final_consensus.fasta",
    output:
        consensus = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
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

                samtools mpileup -aa -A -d 0 -Q 0 \
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
        ref = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        r1  = get_r1,
        r2  = get_r2,
    output:
        bam = OUTDIR + "/{sample}/{line}/final/mapped.bam",
        bai = OUTDIR + "/{sample}/{line}/final/mapped.bam.bai",
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
        consensus = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
    output:
        polished = OUTDIR + "/{sample}/{line}/{assembler}_polish/{assembler}_guided_polished.fasta",
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
        consensus = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
    output:
        paf = OUTDIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.paf",
        bam = OUTDIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.bam",
        vcf = OUTDIR + "/{sample}/{line}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.vcf.gz",
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
        ref      = OUTDIR + "/{sample}/reference/ref.fasta",
        iter_stats = OUTDIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        ivar     = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished = expand(
            OUTDIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
    output:
        alignment = OUTDIR + "/{sample}/{line}/final_alignment/aligned.fasta",
    params:
        assemblers = ASSEMBLERS,
        outdir     = OUTDIR,
    threads: THREADS
    run:
        import subprocess, os, glob, re, tempfile

        # Discover round directories (created by iterative_mapping)
        iter_base = os.path.join(
            str(params.outdir), wildcards.sample, wildcards.line, "iterative"
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
        alignment = OUTDIR + "/{sample}/{line}/final_alignment/aligned.fasta",
    output:
        report    = OUTDIR + "/{sample}/{line}/homoplasy/homoplasy_report.tsv",
        distances = OUTDIR + "/{sample}/{line}/homoplasy/pairwise_distances.tsv",
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
        ref        = OUTDIR + "/{sample}/reference/ref.fasta",
        iter_stats = OUTDIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        consensus  = OUTDIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished   = expand(
            OUTDIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
        comparison = expand(
            OUTDIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_vs_consensus/{assembler}_vs_consensus.paf",
            assembler=ASSEMBLERS,
        ),
        alignment  = OUTDIR + "/{sample}/{line}/final_alignment/aligned.fasta",
        final_bam  = lambda wc: _analysis_bam(
            f"{OUTDIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        ivar_bam = lambda wc: _analysis_bam(
            f"{OUTDIR}/{wc.sample}/{wc.line}/ivar", wc.line
        ),
    output:
        summary = OUTDIR + "/{sample}/{line}/summary.tsv",
    params:
        check_stops = CHECK_STOPS,
        assemblers  = ASSEMBLERS,
        outdir      = OUTDIR,
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
            str(params.outdir), wildcards.sample, wildcards.line, "iterative"
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
# IQ-TREE COMBINED PHYLOGENETIC ANALYSIS
# ═══════════════════════════════════════════════════════════════

rule iqtree_combined:
    """Build ML tree per segment across all samples + GenBank references using IQ-TREE."""
    input:
        ivar = expand(
            OUTDIR + "/{sample}/{{line}}/ivar/final_consensus.fasta",
            sample=SAMPLES.keys(),
        ),
        polished = expand(
            OUTDIR + "/{sample}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            sample=SAMPLES.keys(),
            assembler=ASSEMBLERS,
        ),
        genbank = GENBANK_DIR + "/{segment}.full_len.cds.clean.uniq.fa",
    output:
        alignment = OUTDIR + "/iqtree/{line}/{segment}/aligned.fasta",
        treefile  = OUTDIR + "/iqtree/{line}/{segment}/iqtree.treefile",
        iqtree    = OUTDIR + "/iqtree/{line}/{segment}/iqtree.iqtree",
        log       = OUTDIR + "/iqtree/{line}/{segment}/iqtree.log",
    params:
        samples    = list(SAMPLES.keys()),
        assemblers = ASSEMBLERS,
        outdir     = OUTDIR,
    threads: THREADS
    run:
        import subprocess, os, tempfile

        def extract_segment(fasta_path, segment, sample):
            """Extract the sequence for a given segment from a multi-segment FASTA."""
            seqs = read_fasta(fasta_path)
            target = f"{segment}_{sample}"
            if target in seqs:
                return seqs[target]
            # Fallback: find key starting with segment prefix
            for name, seq in seqs.items():
                if name.startswith(f"{segment}_"):
                    return seq
            return None

        seg = wildcards.segment
        line = wildcards.line
        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        # ── Collect sequences ──
        collected = {}

        # ivar consensus per sample
        for sample in params.samples:
            ivar_path = os.path.join(
                str(params.outdir), sample, line, "ivar", "final_consensus.fasta"
            )
            if os.path.exists(ivar_path):
                seq = extract_segment(ivar_path, seg, sample)
                if seq and seq.replace("N", ""):
                    collected[f"{sample}_ivar"] = seq

        # Polished consensus per sample × assembler
        for sample in params.samples:
            for assembler in params.assemblers:
                pol_path = os.path.join(
                    str(params.outdir), sample, line,
                    f"{assembler}_polish",
                    f"{assembler}_guided_polished.fasta",
                )
                if os.path.exists(pol_path):
                    seq = extract_segment(pol_path, seg, sample)
                    if seq and seq.replace("N", ""):
                        collected[f"{sample}_{assembler}"] = seq

        # GenBank references
        genbank_seqs = read_fasta(str(input.genbank))
        for name, seq in genbank_seqs.items():
            collected[f"genbank_{name}"] = seq

        if len(collected) < 4:
            # Not enough sequences for meaningful tree; write stub outputs
            with open(str(output.alignment), "w") as f:
                for name, seq in collected.items():
                    f.write(f">{name}\n{seq}\n")
            with open(str(output.treefile), "w") as f:
                f.write("(no_tree);\n")
            with open(str(output.iqtree), "w") as f:
                f.write("# Too few sequences for IQ-TREE analysis\n")
            with open(str(output.log), "w") as f:
                f.write("# Too few sequences for IQ-TREE analysis\n")
            print(f"[WARN] {seg}/{line}: only {len(collected)} sequences, skipping IQ-TREE")
            return

        # ── Write unaligned FASTA ──
        unaligned_path = os.path.join(out_dir, "unaligned.fasta")
        with open(unaligned_path, "w") as f:
            for name, seq in collected.items():
                f.write(f">{name}\n{seq}\n")

        print(f"[INFO] {seg}/{line}: {len(collected)} sequences collected for IQ-TREE")

        # ── MAFFT alignment ──
        result = subprocess.run(
            ["mafft", "--auto", "--thread", str(threads), unaligned_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr}")

        with open(str(output.alignment), "w") as f:
            f.write(result.stdout)

        # ── IQ-TREE ──
        prefix = os.path.join(out_dir, "iqtree")
        cmd = [
            "iqtree",
            "-s", str(output.alignment),
            "-m", "MFP",          # ModelFinder Plus
            "-bb", "1000",        # ultrafast bootstrap
            "-alrt", "1000",      # SH-aLRT branch test
            "-T", str(threads),
            "--prefix", prefix,
            "-redo",              # overwrite previous run
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"IQ-TREE failed: {result.stderr}")

        print(f"[INFO] {seg}/{line}: IQ-TREE completed → {output.treefile}")


rule iqtree_homoplasy_summary:
    """Parse IQ-TREE outputs and compute CI/RI/HI via Fitch parsimony."""
    input:
        iqtree_files = expand(
            OUTDIR + "/iqtree/{line}/{segment}/iqtree.iqtree",
            line=LINES,
            segment=SEGMENTS,
        ),
        tree_files = expand(
            OUTDIR + "/iqtree/{line}/{segment}/iqtree.treefile",
            line=LINES,
            segment=SEGMENTS,
        ),
        aln_files = expand(
            OUTDIR + "/iqtree/{line}/{segment}/aligned.fasta",
            line=LINES,
            segment=SEGMENTS,
        ),
    output:
        summary       = OUTDIR + "/iqtree/iqtree_summary.tsv",
        taxon_report  = OUTDIR + "/iqtree/taxon_homoplasy.tsv",
    run:
        import os, re
        from collections import Counter
        from Bio import Phylo
        from io import StringIO

        # ── Fitch parsimony (bottom-up) ──
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

        # ── Fitch top-down assignment ──
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

        # ── Build parent lookup ──
        def build_parent_map(tree):
            parents = {}
            for clade in tree.find_clades(order="level"):
                for child in clade.clades:
                    parents[id(child)] = clade
            return parents

        # ── Combined analysis: global metrics + per-taxon ──
        def analyse_tree(tree_path, aln_path):
            """Return (global_metrics_dict, per_taxon_stats_dict).

            Per-taxon stats include patristic-distance-aware convergence:
              - For each convergent site, the patristic distance to the nearest
                other taxon sharing the same derived state is recorded.
              - 'mean_conv_dist': mean patristic distance of convergent changes
              - 'near_convergent' / 'distant_convergent': count of convergent
                changes where the nearest sharer is below/above the median
                pairwise patristic distance (i.e. local vs distant homoplasy).
            """
            na_global = {
                "parsimony_score": "NA", "ci": "NA", "ri": "NA", "hi": "NA",
            }
            with open(tree_path) as f:
                tree_str = f.read().strip()
            if tree_str.startswith("(no_tree)") or not tree_str:
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
            if content.startswith("# Too few"):
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

        # ── Build summary ──
        path_map = {}
        for iq_path in input.iqtree_files:
            parts = str(iq_path).split("/")
            segment = parts[-2]
            line = parts[-3]
            base_dir = os.path.dirname(str(iq_path))
            path_map[(line, segment)] = base_dir

        summary_rows = []
        all_taxon_rows = []

        for (line, segment), base_dir in sorted(path_map.items()):
            iq_path = os.path.join(base_dir, "iqtree.iqtree")
            tree_path = os.path.join(base_dir, "iqtree.treefile")
            aln_path = os.path.join(base_dir, "aligned.fasta")

            info = parse_iqtree(iq_path)
            global_metrics, taxon_stats = analyse_tree(tree_path, aln_path)
            info.update(global_metrics)
            summary_rows.append((line, segment, info))

            # ── Per-taxon rows (sample taxa + GenBank aggregate) ──
            genbank_agg = {
                "terminal_changes": 0, "convergent": 0,
                "autapomorphic": 0, "near_convergent": 0,
                "distant_convergent": 0, "conv_dist_sum": 0.0,
                "conv_dist_n": 0, "count": 0,
            }

            for taxon, ts in sorted(taxon_stats.items()):
                is_genbank = taxon.startswith("genbank_")
                if is_genbank:
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

                # Parse sample and method from taxon name (e.g. "B1_ivar")
                parts_t = taxon.rsplit("_", 1)
                if len(parts_t) == 2:
                    sample, method = parts_t
                else:
                    sample, method = taxon, "unknown"

                tc = ts["terminal_changes"]
                all_taxon_rows.append({
                    "line": line, "segment": segment,
                    "taxon": taxon, "sample": sample, "method": method,
                    "terminal_changes": tc,
                    "convergent_changes": ts["convergent"],
                    "near_convergent": ts["near_convergent"],
                    "distant_convergent": ts["distant_convergent"],
                    "autapomorphic_changes": ts["autapomorphic"],
                    "total_variable_sites": ts["total_variable"],
                    "homoplasy_fraction": f"{ts['homoplasy_frac']:.4f}",
                    "mean_convergent_distance": f"{ts['mean_conv_dist']:.4f}",
                })

            # GenBank aggregate row for comparison
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
                    "taxon": f"genbank_mean (n={gc})",
                    "sample": "genbank", "method": "reference",
                    "terminal_changes": round(gtc / gc, 1),
                    "convergent_changes": round(gconv / gc, 1),
                    "near_convergent": round(gnear / gc, 1),
                    "distant_convergent": round(gdist / gc, 1),
                    "autapomorphic_changes": round(gauto / gc, 1),
                    "total_variable_sites": taxon_stats[
                        next(k for k in taxon_stats if k.startswith("genbank_"))
                    ]["total_variable"] if taxon_stats else 0,
                    "homoplasy_fraction": f"{gfrac:.4f}",
                    "mean_convergent_distance": f"{gmean_d:.4f}",
                })

        # ── Write global summary ──
        hdr = [
            "line", "segment", "num_sequences", "alignment_length",
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
            for line, segment, info in summary_rows:
                f.write("\t".join([
                    line, segment,
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

        print(f"IQ-TREE summary written to {output.summary}")
        print(f"Per-taxon homoplasy written to {output.taxon_report}")
