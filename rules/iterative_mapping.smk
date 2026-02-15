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


