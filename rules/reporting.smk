rule sample_report:
    """Generate a comprehensive HTML report consolidating all sample results."""
    input:
        summary      = SAMPLE_DIR + "/{sample}/{line}/summary.tsv",
        conv_stats   = SAMPLE_DIR + "/{sample}/{line}/iterative/convergence_stats.tsv",
        final_bam    = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        homoplasy    = SAMPLE_DIR + "/{sample}/{line}/homoplasy/homoplasy_report.tsv",
        pairwise     = SAMPLE_DIR + "/{sample}/{line}/homoplasy/pairwise_distances.tsv",
        ranking      = SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
        site_homo    = SAMPLE_DIR + "/{sample}/{line}/best_consensus/site_homoplasies.tsv",
        best_fasta   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/best_consensus.fasta",
        validation   = SAMPLE_DIR + "/{sample}/{line}/best_consensus/validation.tsv",
        read_support = SAMPLE_DIR + "/{sample}/{line}/best_consensus/read_support.tsv",
        placements   = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/iqtree/{segment}/consensus_placements.tsv",
            segment=SEGMENTS,
        ),
    output:
        report = SAMPLE_DIR + "/{sample}/{line}/report.html",
    params:
        segments = SEGMENTS,
    shell:
        """
        python {workflow.basedir}/scripts/sample_report.py \
            --sample {wildcards.sample} --line {wildcards.line} \
            --summary {input.summary} \
            --convergence {input.conv_stats} \
            --final-bam {input.final_bam} \
            --homoplasy {input.homoplasy} \
            --pairwise {input.pairwise} \
            --ranking {input.ranking} \
            --site-homoplasy {input.site_homo} \
            --best-fasta {input.best_fasta} \
            --validation {input.validation} \
            --read-support {input.read_support} \
            --placements {input.placements} \
            --segments {params.segments} \
            --output {output.report}
        """


rule sample_readme:
    """Generate a README describing the directory structure for a sample/line."""
    input:
        report = SAMPLE_DIR + "/{sample}/{line}/report.html",
    output:
        readme = SAMPLE_DIR + "/{sample}/{line}/README.md",
    run:
        assemblers = ASSEMBLERS
        content = f"""\
# {wildcards.sample} / {wildcards.line}

Output directory for sample **{wildcards.sample}**, analysis line **{wildcards.line}**.

## Directory structure

### Top-level files

| File | Description |
|------|-------------|
| `report.html` | Comprehensive HTML report with coverage plots, convergence graphs, consensus quality, homoplasy analysis, and phylogenetic placement |
| `summary.tsv` | Per-segment coverage statistics (mapped reads, mean depth, completeness) across all pipeline stages |
| `README.md` | This file |

### `iterative/` — Iterative reference mapping

Reads are mapped against the auto-selected reference and the consensus is
refined over multiple rounds until convergence.

| File | Description |
|------|-------------|
| `convergence_stats.tsv` | Per-round mapped reads and genome completeness |
| `final_consensus.fasta` | Consensus sequence after the final iteration |
| `analysis.bam` | BAM from the final iteration (used downstream) |
| `round_N/mapped.bam` | Alignments for iteration N |
| `round_N/reference.fasta` | Reference used in iteration N |
| `round_N/consensus.fasta` | Consensus called after iteration N |
| `round_N/variants.vcf.gz` | Filtered variant calls for iteration N |
| `round_N/variants.raw.vcf.gz` | Unfiltered variant calls |
| `round_N/variants.inframe.vcf.gz` | In-frame filtered variants |

### `ivar/` — iVar consensus calling

Reads are re-mapped to the iterative consensus and iVar generates a new
consensus using a frequency-based approach.

| File | Description |
|------|-------------|
| `mapped.bam` | Alignments to the iterative consensus |
| `reference.fasta` | Reference used (iterative consensus) |
| `final_consensus.fasta` | iVar consensus sequence |

### `final/` — Final mapping

Reads are mapped one last time to the iVar consensus to produce the
definitive alignment used for depth plots and read-support analysis.

| File | Description |
|------|-------------|
| `mapped.bam` | Final alignments |
| `reference.fasta` | iVar consensus used as reference |

### `assembly_reads/` — Mapped read extraction

| File | Description |
|------|-------------|
| `mapped_R1.fastq.gz` | Forward reads that mapped during iterative mapping |
| `mapped_R2.fastq.gz` | Reverse reads that mapped during iterative mapping |

### De novo assemblers

Each assembler has two directories: the raw assembly output and a
`*_polish/` directory where contigs are aligned to the iVar consensus
and polished.

"""
        for asm in assemblers:
            content += f"""\
#### `{asm}/` — {asm} raw assembly output

Raw contigs/scaffolds produced by {asm}. Contains assembler-specific
files (contigs, scaffolds, logs, graphs).

#### `{asm}_polish/` — Assembly-guided polishing

| File | Description |
|------|-------------|
| `{asm}_guided_polished.fasta` | Polished consensus: iVar consensus corrected by assembly evidence |
| `assembly_vs_consensus.bam` | Assembly contigs aligned to the iVar consensus |
| `variants.vcf.gz` | Variants between assembly and iVar consensus |

"""
        content += """\
### `final_alignment/` — Multi-method alignment

| File | Description |
|------|-------------|
| `aligned.fasta` | MAFFT alignment of all consensus methods (iVar, iterative rounds, assembler-polished, reference) for the sample |

### `homoplasy/` — Homoplasy analysis

| File | Description |
|------|-------------|
| `homoplasy_report.tsv` | Delta score (treelikeness) and parsimony metrics per segment |
| `pairwise_distances.tsv` | Pairwise patristic distances between consensus methods |

### `best_consensus/` — Consensus ranking and selection

| File | Description |
|------|-------------|
| `consensus_ranking.tsv` | All consensus methods ranked by homoplastic site count per segment |
| `best_consensus.fasta` | The top-ranked consensus for each segment |
| `hybrid_consensus.fasta` | Hybrid consensus combining best bases from multiple methods |
| `site_homoplasies.tsv` | Per-site homoplasy classification (convergent / autapomorphic) |
| `read_support.tsv` | Allele frequencies from the final BAM at disputed (homoplastic) sites |
| `validation.tsv` | Comparison of best/hybrid consensus vs. original rank-1 method |
| `validation/` | IQ-TREE runs for best and hybrid consensus placement |

### `iqtree/` — Phylogenetic placement

Per-segment subdirectories (`L/`, `M/`, `S/`), each containing:

| File | Description |
|------|-------------|
| `aligned.fasta` | Backbone alignment with GenBank references |
| `iqtree.treefile` | Maximum-likelihood tree |
| `iqtree.iqtree` | IQ-TREE report (model, log-likelihood, tree length) |
| `sample_seqs.fasta` | Sample consensus sequences used as queries |
| `consensus_placements.tsv` | Placement summary for each consensus method |
| `placements/<method>/` | Per-method IQ-TREE placement run |
"""
        with open(str(output.readme), "w") as f:
            f.write(content)

        print(f"[INFO] {wildcards.line}/{wildcards.sample}: README written → {output.readme}")


rule sample_root_readme:
    """Generate a README for the sample root directory."""
    input:
        readmes = expand(
            SAMPLE_DIR + "/{{sample}}/{line}/README.md",
            line=LINES,
        ),
    output:
        readme = SAMPLE_DIR + "/{sample}/README.md",
    run:
        lines_str = ", ".join(f"`{l}`" for l in LINES)
        content = f"""\
# {wildcards.sample}

Output directory for sample **{wildcards.sample}**.

## Directory structure

| Directory | Description |
|-----------|-------------|
| `host_filtered/` | Host-depleted reads (human and bank vole removed via bowtie2) |
| `autoref/` | Automatic reference selection — per-segment best GenBank reference |
| `reference/` | Combined reference FASTA (L + M + S best references concatenated) |
| `dedup/` | Analysis line with PCR duplicate removal (see `dedup/README.md`) |
| `nodedup/` | Analysis line without duplicate removal (see `nodedup/README.md`) |

### `host_filtered/`

| File | Description |
|------|-------------|
| `clean_R1.fastq.gz` | Forward reads after host removal |
| `clean_R2.fastq.gz` | Reverse reads after host removal |
| `human_bowtie2.log` | Bowtie2 log for human genome filtering |
| `bankvole_bowtie2.log` | Bowtie2 log for bank vole genome filtering |

### `autoref/`

Per-segment subdirectories (`L/`, `M/`, `S/`), each containing:

| File | Description |
|------|-------------|
| `best.fasta` | Best-matching GenBank reference for this segment |
| `ref_stats.tsv` | Mapping statistics for all candidate references |
| `mapped.bam` | Alignments used for reference evaluation |

Top-level file:

| File | Description |
|------|-------------|
| `best_refs.fasta` | Combined best references (all segments) |

### `reference/`

| File | Description |
|------|-------------|
| `ref.fasta` | Combined reference FASTA used as starting point for iterative mapping |
| `ref.fasta.fai` | samtools faidx index |

### Analysis lines ({lines_str})

Each analysis line runs the full pipeline independently. The `dedup` line
removes PCR duplicates before consensus calling; the `nodedup` line retains
all reads. See the README.md inside each line directory for a detailed
description of all subdirectories and files.
"""
        with open(str(output.readme), "w") as f:
            f.write(content)

        print(f"[INFO] {wildcards.sample}: root README written → {output.readme}")


rule read_support:
    """Check read-level allele frequencies at sites flagged as homoplastic."""
    input:
        site_homo  = SAMPLE_DIR + "/{sample}/{line}/best_consensus/site_homoplasies.tsv",
        final_bam  = lambda wc: _analysis_bam(
            f"{SAMPLE_DIR}/{wc.sample}/{wc.line}/final", wc.line
        ),
        ivar_cons  = SAMPLE_DIR + "/{sample}/{line}/ivar/final_consensus.fasta",
        placements = expand(
            SAMPLE_DIR + "/{{sample}}/{{line}}/iqtree/{segment}/consensus_placements.tsv",
            segment=SEGMENTS,
        ),
    output:
        support = SAMPLE_DIR + "/{sample}/{line}/best_consensus/read_support.tsv",
    run:
        import subprocess, os, re
        from collections import defaultdict

        sample = wildcards.sample
        line = wildcards.line

        # ── Read site_homoplasies.tsv ──
        sites = []
        with open(str(input.site_homo)) as f:
            header = f.readline()
            if header.startswith("#"):
                # No data
                with open(str(output.support), "w") as out:
                    out.write("# No homoplastic sites to validate\n")
                return
            for row_line in f:
                fields = row_line.strip().split("\t")
                if len(fields) >= 9:
                    sites.append({
                        "sample": fields[0], "line": fields[1],
                        "segment": fields[2], "method": fields[3],
                        "position": int(fields[4]),  # 1-based alignment position
                        "sample_state": fields[5],
                        "ancestral_state": fields[6],
                        "classification": fields[7],
                        "convergent_partners": fields[8],
                    })

        if not sites:
            with open(str(output.support), "w") as out:
                out.write("# No homoplastic sites to validate\n")
            return

        # ── Build alignment position → ungapped position mapping ──
        # For each (segment, method), read the aligned sequence and build
        # aln_pos (1-based) → ungapped_pos (1-based) mapping
        placement_dirs = {}  # (segment, method) -> base directory
        for tsv_path in input.placements:
            with open(str(tsv_path)) as f:
                hdr = f.readline()
                if hdr.startswith("#"):
                    continue
                for row_line in f:
                    fields = row_line.strip().split("\t")
                    if len(fields) >= 8:
                        seg, method, taxon = fields[2], fields[3], fields[4]
                        aln_p = fields[6]
                        placement_dirs[(seg, method)] = {
                            "aln": aln_p, "taxon": taxon,
                        }

        aln_to_ungapped = {}  # (segment, method) -> {aln_pos: ungapped_pos}
        for (seg, method), info in placement_dirs.items():
            aln_path = info["aln"]
            taxon = info["taxon"]
            if not os.path.exists(aln_path):
                continue
            aln_seqs = read_fasta(aln_path)
            if taxon not in aln_seqs:
                continue
            aligned_seq = aln_seqs[taxon]
            mapping = {}
            ungapped = 0
            for i, base in enumerate(aligned_seq):
                if base != "-":
                    ungapped += 1
                    mapping[i + 1] = ungapped  # 1-based
            aln_to_ungapped[(seg, method)] = mapping

        # ── Get contig names from iVar consensus ──
        ivar_seqs = read_fasta(str(input.ivar_cons))
        seg_to_contig = {}
        for contig_name in ivar_seqs:
            for seg in SEGMENTS:
                if seg in contig_name:
                    seg_to_contig[seg] = contig_name
                    break

        # ── Run samtools mpileup per contig ──
        pileup_data = {}  # contig -> {pos: {A:n, C:n, G:n, T:n, del:n, depth:n}}
        for seg, contig in seg_to_contig.items():
            result = subprocess.run(
                ["samtools", "mpileup", "-a", "-A", "-Q", "0",
                 "-r", contig, str(input.final_bam)],
                capture_output=True, text=True,
            )
            contig_pileup = {}
            for pline in result.stdout.strip().split("\n"):
                if not pline:
                    continue
                parts = pline.split("\t")
                if len(parts) < 5:
                    continue
                pos = int(parts[1])
                ref_base = parts[2].upper()
                depth = int(parts[3])
                bases_str = parts[4].upper()

                counts = {"A": 0, "C": 0, "G": 0, "T": 0, "del": 0}
                i = 0
                while i < len(bases_str):
                    c = bases_str[i]
                    if c in (".", ","):
                        counts[ref_base] = counts.get(ref_base, 0) + 1
                    elif c in ("A", "C", "G", "T"):
                        counts[c] += 1
                    elif c == "*" or c == "#":
                        counts["del"] += 1
                    elif c == "^":
                        i += 1  # skip quality char
                    elif c in ("+", "-"):
                        # Indel: read the length digits then skip that many bases
                        j = i + 1
                        while j < len(bases_str) and bases_str[j].isdigit():
                            j += 1
                        indel_len = int(bases_str[i+1:j]) if j > i + 1 else 0
                        i = j + indel_len - 1
                        if c == "-":
                            counts["del"] += 1
                    i += 1

                counts["depth"] = depth
                counts["ref_base"] = ref_base
                contig_pileup[pos] = counts
            pileup_data[seg] = contig_pileup

        # ── Annotate each site with read support ──
        support_rows = []
        for site in sites:
            seg = site["segment"]
            method = site["method"]
            aln_pos = site["position"]

            # Map alignment position to ungapped position
            mapping = aln_to_ungapped.get((seg, method), {})
            ungapped_pos = mapping.get(aln_pos)
            if ungapped_pos is None:
                continue

            # Look up pileup
            contig_pileup = pileup_data.get(seg, {})
            pile = contig_pileup.get(ungapped_pos, {})
            if not pile:
                continue

            depth = pile.get("depth", 0)
            ref_base = pile.get("ref_base", "N")
            cons_base = site["sample_state"].upper()
            cons_count = pile.get(cons_base, 0)
            cons_freq = cons_count / depth if depth > 0 else 0.0

            if cons_freq > 0.8:
                verdict = "supported"
            elif cons_freq >= 0.5:
                verdict = "weak"
            else:
                verdict = "contradicted"

            support_rows.append({
                "sample": sample, "line": line,
                "segment": seg, "method": method,
                "aln_position": aln_pos,
                "contig_position": ungapped_pos,
                "sample_state": site["sample_state"],
                "ancestral_state": site["ancestral_state"],
                "classification": site["classification"],
                "depth": depth,
                "ref_base": ref_base,
                "A_count": pile.get("A", 0),
                "C_count": pile.get("C", 0),
                "G_count": pile.get("G", 0),
                "T_count": pile.get("T", 0),
                "del_count": pile.get("del", 0),
                "consensus_allele_freq": f"{cons_freq:.3f}",
                "read_support_verdict": verdict,
            })

        # Write output
        sup_hdr = [
            "sample", "line", "segment", "method",
            "aln_position", "contig_position",
            "sample_state", "ancestral_state", "classification",
            "depth", "ref_base",
            "A_count", "C_count", "G_count", "T_count", "del_count",
            "consensus_allele_freq", "read_support_verdict",
        ]
        with open(str(output.support), "w") as f:
            f.write("\t".join(sup_hdr) + "\n")
            for r in support_rows:
                f.write("\t".join([str(r[c]) for c in sup_hdr]) + "\n")

        print(f"[INFO] {line}/{sample}: read support for {len(support_rows)} sites → {output.support}")


rule cross_sample_report:
    """Generate a global HTML overview comparing all samples."""
    input:
        summaries = expand(
            SAMPLE_DIR + "/{sample}/{line}/summary.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
        rankings = expand(
            SAMPLE_DIR + "/{sample}/{line}/best_consensus/consensus_ranking.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
        validations = expand(
            SAMPLE_DIR + "/{sample}/{line}/best_consensus/validation.tsv",
            sample=SAMPLES.keys(), line=LINES,
        ),
    output:
        report = COMMON_DIR + "/overview_report.html",
    shell:
        """
        python {workflow.basedir}/scripts/cross_sample_report.py \
            --summaries {input.summaries} \
            --rankings {input.rankings} \
            --validations {input.validations} \
            --output {output.report}
        """
