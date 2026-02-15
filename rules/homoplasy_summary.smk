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


