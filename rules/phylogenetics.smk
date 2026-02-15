# ═══════════════════════════════════════════════════════════════
# IQ-TREE TWO-PHASE PHYLOGENETIC ANALYSIS
# Phase 1: backbone tree from GenBank references (per segment)
# Phase 2: place each sample's consensus sequences onto backbone
# ═══════════════════════════════════════════════════════════════

rule iqtree_backbone:
    """Phase 1: Build backbone ML tree from GenBank references only (once per segment)."""
    input:
        genbank = COMMON_DIR + "/genbank_subsampled/{segment}.subsampled.fa",
    output:
        alignment = COMMON_DIR + "/iqtree/backbone/{segment}/aligned.fasta",
        treefile  = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.treefile",
        iqtree    = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.iqtree",
        log       = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.log",
    threads: THREADS
    run:
        import subprocess, os

        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        genbank_seqs = read_fasta(str(input.genbank))

        if len(genbank_seqs) < 4:
            with open(str(output.alignment), "w") as f:
                for name, seq in genbank_seqs.items():
                    f.write(f">{name}\n{seq}\n")
            with open(str(output.treefile), "w") as f:
                f.write("(no_tree);\n")
            with open(str(output.iqtree), "w") as f:
                f.write("# Too few sequences for IQ-TREE backbone\n")
            with open(str(output.log), "w") as f:
                f.write("# Too few sequences for IQ-TREE backbone\n")
            print(f"[WARN] {wildcards.segment}: only {len(genbank_seqs)} GenBank refs, skipping backbone")
            return

        # Write unaligned GenBank FASTA
        unaligned_path = os.path.join(out_dir, "unaligned.fasta")
        with open(unaligned_path, "w") as f:
            for name, seq in genbank_seqs.items():
                f.write(f">{name}\n{seq}\n")

        print(f"[INFO] {wildcards.segment}: {len(genbank_seqs)} GenBank refs for backbone")

        # MAFFT alignment
        result = subprocess.run(
            ["mafft", "--auto", "--thread", str(threads), unaligned_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr}")
        with open(str(output.alignment), "w") as f:
            f.write(result.stdout)

        # IQ-TREE full ML search with ModelFinder
        prefix = os.path.join(out_dir, "backbone")
        cmd = [
            "iqtree",
            "-s", str(output.alignment),
            "-m", "MFP",          # ModelFinder Plus
            "-bb", "1000",        # ultrafast bootstrap
            "-alrt", "1000",      # SH-aLRT branch test
            "-T", str(threads),
            "--prefix", prefix,
            "-redo",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"IQ-TREE backbone failed: {result.stderr}")

        print(f"[INFO] {wildcards.segment}: backbone tree completed → {output.treefile}")


rule iqtree_place_sample:
    """Phase 2: Place a sample's consensus sequences onto the fixed backbone tree."""
    input:
        backbone_tree = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.treefile",
        backbone_aln  = COMMON_DIR + "/iqtree/backbone/{segment}/aligned.fasta",
        backbone_iq   = COMMON_DIR + "/iqtree/backbone/{segment}/backbone.iqtree",
        ivar = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        polished = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/{assembler}_polish/{assembler}_guided_polished.fasta",
            assembler=ASSEMBLERS,
        ),
        iter_stats = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
    output:
        alignment  = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/aligned.fasta",
        treefile   = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.treefile",
        iqtree     = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.iqtree",
        log        = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/iqtree.log",
        placements = SAMPLE_DIR + "/{sample}/{line}/iqtree/{segment}/consensus_placements.tsv",
    params:
        assemblers = ASSEMBLERS,
        sample_dir = SAMPLE_DIR,
    threads: THREADS
    run:
        import subprocess, os, glob, re, shutil

        def extract_segment(fasta_path, segment, sample):
            """Extract the sequence for a given segment from a multi-segment FASTA."""
            seqs = read_fasta(fasta_path)
            target = f"{segment}_{sample}"
            if target in seqs:
                return seqs[target]
            for name, seq in seqs.items():
                if f"_{segment}_{sample}" in name or name.startswith(f"{segment}_"):
                    return seq
            return None

        def run_single_placement(taxon_name, seq, backbone_aln, backbone_tree,
                                 model, sub_dir, n_threads):
            """Place one consensus sequence on the backbone tree.

            Returns (iq_path, tree_path, aln_path, log_path) or None on failure.
            """
            os.makedirs(sub_dir, exist_ok=True)

            # Write single-sequence FASTA
            query_fasta = os.path.join(sub_dir, "query.fasta")
            with open(query_fasta, "w") as fh:
                fh.write(f">{taxon_name}\n{seq}\n")

            # MAFFT --add: add this one sequence to the backbone alignment
            res = subprocess.run(
                ["mafft", "--add", query_fasta, "--keeplength",
                 "--thread", str(n_threads), backbone_aln],
                capture_output=True, text=True,
            )
            if res.returncode != 0:
                print(f"[WARN] MAFFT --add failed for {taxon_name}: {res.stderr}")
                return None

            aln_path = os.path.join(sub_dir, "aligned.fasta")
            with open(aln_path, "w") as fh:
                fh.write(res.stdout)

            # IQ-TREE with constraint backbone
            prefix = os.path.join(sub_dir, "iqtree")
            cmd = [
                "iqtree",
                "-s", aln_path,
                "-g", backbone_tree,
                "-m", model,
                "-T", str(n_threads),
                "--prefix", prefix,
                "-redo",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                print(f"[WARN] IQ-TREE failed for {taxon_name}: {res.stderr}")
                return None

            return (
                os.path.join(sub_dir, "iqtree.iqtree"),
                os.path.join(sub_dir, "iqtree.treefile"),
                aln_path,
                os.path.join(sub_dir, "iqtree.log"),
            )

        seg = wildcards.segment
        line = wildcards.line
        sample = wildcards.sample
        out_dir = os.path.dirname(str(output.alignment))
        os.makedirs(out_dir, exist_ok=True)

        # Check backbone is valid
        with open(str(input.backbone_tree)) as f:
            tree_str = f.read().strip()
        if tree_str.startswith("(no_tree)") or not tree_str:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# No backbone tree available\n")
            with open(str(output.placements), "w") as f:
                f.write("# No backbone tree available\n")
            print(f"[WARN] {seg}/{line}/{sample}: no backbone tree, skipping placement")
            return

        # Parse model from backbone IQ-TREE report
        best_model = "MFP"
        with open(str(input.backbone_iq)) as f:
            for bline in f:
                m = re.search(r"Best-fit model according to BIC:\s+(.+)", bline)
                if m:
                    best_model = m.group(1).strip()
                    break

        # ── Collect sample consensus sequences ──
        # Dict of method_label -> (taxon_name, sequence)
        collected = {}

        # Reference used for iterative mapping
        ref_path = os.path.join(
            str(params.sample_dir), sample, line, "iterative", "round_1", "reference.fasta"
        )
        if os.path.exists(ref_path):
            seq = extract_segment(ref_path, seg, sample)
            if seq and seq.replace("N", ""):
                collected["reference"] = (f"{sample}_reference", seq)

        # Iterative mapping rounds
        iter_base = os.path.join(str(params.sample_dir), sample, line, "iterative")
        round_dirs = sorted(
            [d for d in glob.glob(os.path.join(iter_base, "round_*")) if os.path.isdir(d)],
            key=lambda d: int(re.search(r"round_(\d+)", d).group(1)),
        )
        for rd in round_dirs:
            n = re.search(r"round_(\d+)", rd).group(1)
            cons_path = os.path.join(rd, "consensus.fasta")
            if os.path.exists(cons_path):
                seq = extract_segment(cons_path, seg, sample)
                if seq and seq.replace("N", ""):
                    collected[f"round_{n}"] = (f"{sample}_round_{n}", seq)

        # iVar consensus
        if os.path.exists(str(input.ivar)):
            seq = extract_segment(str(input.ivar), seg, sample)
            if seq and seq.replace("N", ""):
                collected["ivar"] = (f"{sample}_ivar", seq)

        # Polished consensus per assembler
        for i, assembler in enumerate(params.assemblers):
            pol_path = str(input.polished[i])
            if os.path.exists(pol_path):
                seq = extract_segment(pol_path, seg, sample)
                if seq and seq.replace("N", ""):
                    collected[assembler] = (f"{sample}_{assembler}", seq)

        if len(collected) == 0:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# No sample sequences for this segment\n")
            with open(str(output.placements), "w") as f:
                f.write("# No sample sequences for this segment\n")
            print(f"[WARN] {seg}/{line}/{sample}: no sample sequences, skipping")
            return

        print(f"[INFO] {seg}/{line}/{sample}: {len(collected)} consensus sequences "
              f"to place individually on backbone")

        # ── Individual placements ──
        placements_dir = os.path.join(out_dir, "placements")
        placement_results = []  # (method, taxon_name, iq_path, tree_path, aln_path)

        for method, (taxon_name, seq) in collected.items():
            sub_dir = os.path.join(placements_dir, method)
            paths = run_single_placement(
                taxon_name, seq,
                str(input.backbone_aln), str(input.backbone_tree),
                best_model, sub_dir, threads,
            )
            if paths:
                iq_path, tree_path, aln_path, log_path = paths
                placement_results.append((method, taxon_name, iq_path, tree_path, aln_path))
                print(f"[INFO] {seg}/{line}/{sample}/{method}: individual placement done")
            else:
                print(f"[WARN] {seg}/{line}/{sample}/{method}: placement failed, skipping")

        # ── Populate top-level outputs for backward compatibility ──
        # Prefer the ivar placement; fall back to last successful placement
        main_method = None
        for pref in ["ivar"]:
            for m, tn, iq, tr, al in placement_results:
                if m == pref:
                    main_method = (m, tn, iq, tr, al)
                    break
            if main_method:
                break
        if not main_method and placement_results:
            main_method = placement_results[-1]

        if main_method:
            _, _, main_iq, main_tree, main_aln = main_method
            main_log = main_iq.replace("iqtree.iqtree", "iqtree.log")
            shutil.copy2(main_aln, str(output.alignment))
            shutil.copy2(main_tree, str(output.treefile))
            shutil.copy2(main_iq, str(output.iqtree))
            if os.path.exists(main_log):
                shutil.copy2(main_log, str(output.log))
            else:
                with open(str(output.log), "w") as f:
                    f.write("# See placements/ subdirectory\n")
        else:
            with open(str(output.alignment), "w") as f:
                f.write("")
            for p in [output.treefile, output.iqtree, output.log]:
                with open(str(p), "w") as f:
                    f.write("# All individual placements failed\n")

        # ── Write consensus_placements.tsv ──
        placement_header = [
            "sample", "line", "segment", "method", "taxon_name",
            "tree_path", "alignment_path", "iqtree_path",
        ]
        with open(str(output.placements), "w") as f:
            f.write("\t".join(placement_header) + "\n")
            for method, taxon_name, iq_path, tree_path, aln_path in placement_results:
                f.write("\t".join([
                    sample, line, seg, method, taxon_name,
                    tree_path, aln_path, iq_path,
                ]) + "\n")

        print(f"[INFO] {seg}/{line}/{sample}: all placements completed "
              f"({len(placement_results)} individual)")


