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


