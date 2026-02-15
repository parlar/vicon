
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


