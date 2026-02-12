# ViCon — Viral Consensus Pipeline

ViCon is a Snakemake pipeline for generating high-quality consensus genomes from segmented RNA viruses (e.g. hantaviruses) using paired-end Illumina data. It combines iterative reference-guided mapping, multiple de novo assemblers, assembly-guided polishing, and phylogenetic analysis with homoplasy metrics. ViCon automatically selects the best consensus method per segment based on phylogenetic placement, constructs a hybrid consensus from non-homoplastic bases, validates the result against the backbone tree, and checks read support at disputed positions. Per-sample and cross-sample HTML reports consolidate all quality metrics.

## Requirements

- **Linux (x86_64)**
- [Pixi](https://pixi.sh/) package manager
- [Apptainer](https://apptainer.org/) (required only for Trinity assembly)

All other dependencies are installed automatically by Pixi.

## Installation

```bash
git clone <repository-url>
cd vicon
pixi install
```

This creates a `.pixi/` environment with all tools (samtools, bcftools, minimap2, strobealign, IQ-TREE, etc.) pinned via `pixi.lock`.

## Input data

ViCon expects the following inputs:

### Paired-end FASTQ files

One pair per sample in a single directory, named `{sample}_R1_*.fastq.gz` and `{sample}_R2_*.fastq.gz`. The sample name is extracted from everything before the first underscore.

### Per-sample reference FASTA

A multi-segment FASTA per sample at `{ref_dir}/{sample}_best_refs.fasta`, containing one sequence per genomic segment. Contig names must start with the corresponding segment name from the `segments` config list (e.g. `L`, `M`, `S` for hantaviruses).

Alternatively (and by default), enable `auto_select_refs: true` and provide reference FASTA files per segment via the `genbank_refs` config option (or place them in `genbank_dir`).

### Host genome indices

Bowtie2 indices for host removal. For hantavirus from bank voles, both human and bank vole indices are used sequentially.

## Configuration

Edit `config.yaml` before running. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fastq_dir` | `../clean_fasta` | Directory with input FASTQ files |
| `ref_dir` | `../sample-specific_refgenomes` | Directory with per-sample reference FASTAs |
| `outdir` | `../results` | Output directory |
| `segments` | `[L, M, S]` | Genomic segments to process (e.g. `[genome]` for single-segment viruses) |
| `genbank_dir` | `../genbank_reference` | Reference database directory (fallback for `genbank_refs`) |
| `genbank_refs` | (auto from `genbank_dir`) | Per-segment reference FASTA paths (`L`, `M`, `S`) |
| `human_bt2_index` | — | Bowtie2 index prefix for human genome |
| `bankvole_bt2_index` | — | Bowtie2 index prefix for bank vole genome |
| `aligner` | `strobealign` | Read aligner (`bwa-mem2` or `strobealign`) |
| `max_iterations` | `10` | Maximum iterative mapping rounds |
| `min_iterations` | `2` | Minimum rounds before convergence check |
| `threads` | `8` | Threads per job |
| `lines` | `[dedup, nodedup]` | Analysis lines to run in parallel |
| `enforce_inframe` | `true` | Filter out-of-frame indels |
| `check_stops` | `true` | Filter variants introducing internal stop codons |
| `assemblers` | `[trinity, spades, coronaspades, megahit, iva]` | De novo assemblers to run |
| `denovo_seed` | `false` | De novo assembly to seed iterative mapping |
| `extend_reference` | `false` | Extend consensus using soft-clipped reads |
| `auto_select_refs` | `true` | Automatic reference selection |

All directory paths in `config.yaml` are relative to the `vicon/` directory.

### Running a subset of samples

Either list samples explicitly in `config.yaml`:

```yaml
samples:
  - B1
  - F1
```

Or override on the command line (see below).

## Running the pipeline

### Full run

```bash
cd vicon
pixi run snakemake --cores 64 --use-singularity
```

Or use the provided wrapper:

```bash
./run_snakemake.sh
```

### Subset run

```bash
./run_subset.sh
```

This runs with 20 cores on FASTQ files from `../subset_fastq` and writes results to `../results_subset`.

### Overriding config on the command line

Any config parameter can be overridden:

```bash
pixi run snakemake --cores 32 \
  --config fastq_dir=../subset_fastq outdir=../results_subset denovo_seed=true
```

### Dry run

Preview what will be executed without running anything:

```bash
pixi run snakemake --cores 1 -n
```

## Pipeline stages

1. **Adapter trimming** — fastp automatic adapter detection and quality trimming of raw reads
2. **Host read removal** — Bowtie2 sequential filtering against human and host genomes
3. **Automatic reference selection** (optional) — Pick best reference per segment from a database
4. **Reference preparation** — Match contigs to configured segment names; index for aligner
5. **De novo reference seeding** (optional) — MEGAHIT assembly scaffolded against the reference to build a sample-specific starting reference
6. **Iterative mapping** — Convergence-based iterative alignment with optional reference elongation using soft-clipped reads
7. **De novo assembly** — Up to five assemblers (Trinity via Apptainer, SPAdes, coronaSPAdes, MEGAHIT, IVA)
8. **iVar consensus** — Per-contig consensus calling from remapped reads
9. **Assembly-guided polishing** — Polish consensus using each de novo assembly via minimap2 + bcftools
10. **Multi-sequence alignment** — MAFFT alignment of all consensus stages per sample
11. **Phylogenetic analysis** — IQ-TREE backbone trees per segment from GenBank references, then individual placement of each consensus method onto the backbone
12. **Homoplasy analysis** — Fitch parsimony CI/RI/HI, delta score, per-taxon convergent change classification
13. **Best consensus selection** — Site-level Fitch analysis ranks consensus methods by homoplastic signal; best method selected per segment; phylogenetically-informed hybrid consensus constructed from non-homoplastic bases across methods
14. **Consensus validation** — Best and hybrid consensus sequences placed back on the backbone tree to verify that homoplasy is reduced compared to individual methods
15. **Read support analysis** — Allele frequencies extracted from the final BAM at positions flagged as homoplastic; sites classified as supported (>80% reads agree), weak (50–80%), or contradicted (<50%)
16. **Reporting** — Per-sample HTML report with coverage/depth plots, consensus quality, ranking, homoplasy details, validation, and read support; cross-sample HTML overview with heatmaps and method frequency charts

Stages 6–15 run in parallel for dedup and nodedup analysis lines.

## Output structure

```
results/
  samples/
    {sample}/
      trimmed/                    # Adapter-trimmed FASTQs and fastp QC reports
      host_filtered/              # Host-removed FASTQ files
      reference/                  # Renamed and indexed reference
      denovo_seed/                # (if enabled) Seeded reference
      {line}/                     # dedup/ or nodedup/
        iterative/                # Iterative mapping rounds and final consensus
          round_1/                # Per-round BAM, VCF, consensus
          round_2/
          ...
          final_consensus.fasta
          analysis.bam
          convergence_stats.tsv
        assembly_reads/           # Reads extracted for de novo assembly
        trinity/                  # Trinity assembly
        spades/                   # SPAdes assembly
        coronaspades/             # coronaSPAdes assembly
        megahit/                  # MEGAHIT assembly
        iva/                      # IVA assembly
        ivar/                     # iVar consensus
        final/                    # Final mapping to iVar consensus
        {assembler}_polish/       # Assembly-guided polished consensus
        final_alignment/          # MAFFT alignment of all stages
        homoplasy/                # Delta score and pairwise distances
        iqtree/{segment}/         # Per-segment IQ-TREE placement results
          consensus_placements.tsv  # Manifest of individual placements
          placements/{method}/    # Per-method placement (tree, alignment, IQ-TREE report)
        best_consensus/           # Best consensus selection results
          consensus_ranking.tsv   # Methods ranked by homoplasy metrics
          site_homoplasies.tsv    # Per-site homoplasy classification
          best_consensus.fasta    # Best method per segment
          hybrid_consensus.fasta  # Phylogenetically-informed hybrid
          validation.tsv          # Validation metrics for best/hybrid
          read_support.tsv        # Allele frequencies at disputed sites
        summary.tsv               # Per-segment QC summary
        report.html               # Per-sample HTML report
  common/
    iqtree/
      backbone/{segment}/         # Backbone trees from GenBank references
      iqtree_summary.tsv          # Global phylogenetic metrics
      individual_placement_homoplasy.tsv  # Per-method homoplasy report
    overview_report.html          # Cross-sample HTML overview
```

## Key output files

| File | Description |
|------|-------------|
| `{sample}/{line}/ivar/final_consensus.fasta` | Primary consensus sequence |
| `{sample}/{line}/{asm}_polish/{asm}_guided_polished.fasta` | Assembly-polished consensus |
| `{sample}/{line}/best_consensus/best_consensus.fasta` | Best consensus per segment (lowest homoplasy) |
| `{sample}/{line}/best_consensus/hybrid_consensus.fasta` | Hybrid consensus (non-homoplastic majority rule) |
| `{sample}/{line}/best_consensus/consensus_ranking.tsv` | All methods ranked by homoplasy metrics |
| `{sample}/{line}/best_consensus/site_homoplasies.tsv` | Per-site homoplasy classification |
| `{sample}/{line}/best_consensus/validation.tsv` | Best/hybrid validation vs original method |
| `{sample}/{line}/best_consensus/read_support.tsv` | Allele frequencies at disputed sites |
| `{sample}/{line}/summary.tsv` | Per-segment coverage, depth, completeness, N counts |
| `{sample}/{line}/iterative/convergence_stats.tsv` | Mapped reads and completeness per iteration |
| `{sample}/{line}/homoplasy/homoplasy_report.tsv` | Delta score and site statistics |
| `{sample}/{line}/report.html` | Per-sample HTML report with plots |
| `common/iqtree/backbone/{segment}/backbone.treefile` | Backbone ML tree from GenBank references |
| `common/iqtree/iqtree_summary.tsv` | Model selection, tree statistics, CI/RI/HI |
| `common/iqtree/individual_placement_homoplasy.tsv` | Per-method homoplasy across all samples |
| `common/overview_report.html` | Cross-sample HTML overview with heatmaps |

## Optional features

### De novo reference seeding

Set `denovo_seed: true` to run a quick MEGAHIT assembly before iterative mapping. The de novo contigs are scaffolded against the reference to produce a sample-specific starting reference. This helps when samples diverge substantially from the available reference.

### Reference elongation

Set `extend_reference: true` to extend the consensus beyond original reference boundaries during iterative mapping. Soft-clipped reads at contig termini are used to build consensus extensions, recovering terminal regions not present in the starting reference.

### Automatic reference selection

Enabled by default (`auto_select_refs: true`). Provide reference FASTA files per segment via the `genbank_refs` config option, or place them in `genbank_dir` using the default naming convention. The pipeline maps reads against all candidates and selects the best reference per segment based on completeness and mapped read count.

## License

[To be specified]
