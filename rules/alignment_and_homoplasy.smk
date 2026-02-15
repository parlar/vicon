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


