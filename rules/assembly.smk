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


