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

