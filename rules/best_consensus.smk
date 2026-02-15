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


