# ═══════════════════════════════════════════════════════════════
# GENBANK CDS VALIDATION
# ═══════════════════════════════════════════════════════════════

rule validate_genbank_refs:
    """Validate GenBank reference sequences as proper CDSes (ATG start, no internal stops, length %3==0)."""
    input:
        fasta = get_genbank_ref,
    output:
        fasta  = COMMON_DIR + "/genbank_validated/{segment}.validated.fa",
        report = COMMON_DIR + "/genbank_validated/{segment}.validation_report.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})
        python3 {workflow.basedir}/scripts/validate_cds.py \
            {input.fasta} {output.fasta} {output.report}
        """


# ═══════════════════════════════════════════════════════════════
# AUTOMATIC REFERENCE SELECTION (optional)
# ═══════════════════════════════════════════════════════════════

rule autoref_map:
    """Map host-filtered reads against all GenBank refs for one segment."""
    input:
        r1 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R1.fastq.gz",
        r2 = SAMPLE_DIR + "/{sample}/host_filtered/clean_R2.fastq.gz",
        ref = get_validated_genbank_ref,
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
        ref = get_validated_genbank_ref,
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
            # Rename contig to start with segment name so rename_ref can match it
            seg = wildcards.segment
            f.write(f">{seg}_{chosen}\n")
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
        fasta = get_validated_genbank_ref,
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
