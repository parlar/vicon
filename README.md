# ViCon — Viral Consensus Pipeline

ViCon is a Snakemake pipeline for generating consensus genomes from segmented RNA viruses (e.g. hantaviruses) using paired-end Illumina data. It combines iterative reference-guided mapping, multiple de novo assemblers, assembly-guided polishing, and phylogenetic analysis with homoplasy metrics.

## Requirements

- **Linux (x86_64)**
- [Pixi](https://pixi.sh/) package manager
- Docker (required only for Trinity assembly)

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

A multi-segment FASTA per sample at `{ref_dir}/{sample}_best_refs.fasta`, containing one sequence per genomic segment. Contigs are automatically labelled L, M, S by size.

Alternatively, enable `auto_select_refs: true` and provide a directory of reference databases (one FASTA per segment: `L.full_len.cds.clean.uniq.fa`, `M.full_len.cds.clean.uniq.fa`, `S.full_len.cds.clean.uniq.fa`).

### Host genome indices

Bowtie2 indices for host removal. For hantavirus from bank voles, both human and bank vole indices are used sequentially.

## Configuration

Edit `config.yaml` before running. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fastq_dir` | `../clean_fasta` | Directory with input FASTQ files |
| `ref_dir` | `../sample-specific_refgenomes` | Directory with per-sample reference FASTAs |
| `outdir` | `../results` | Output directory |
| `genbank_dir` | `../genbank_reference` | Reference database for auto-selection |
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
| `auto_select_refs` | `false` | Automatic reference selection |

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
pixi run snakemake --cores 64
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

1. **Host read removal** — Bowtie2 sequential filtering against human and host genomes
2. **Automatic reference selection** (optional) — Pick best reference per segment from a database
3. **Reference preparation** — Label contigs as L, M, S by size; index for aligner
4. **De novo reference seeding** (optional) — MEGAHIT assembly scaffolded against the reference to build a sample-specific starting reference
5. **Iterative mapping** — Convergence-based iterative alignment with optional reference elongation using soft-clipped reads
6. **De novo assembly** — Up to five assemblers (Trinity, SPAdes, coronaSPAdes, MEGAHIT, IVA)
7. **iVar consensus** — Per-contig consensus calling from remapped reads
8. **Assembly-guided polishing** — Polish consensus using each de novo assembly via minimap2 + bcftools
9. **Multi-sequence alignment** — MAFFT alignment of all consensus stages per sample
10. **Phylogenetic analysis** — IQ-TREE ML trees per segment with ModelFinder, ultrafast bootstrap, and SH-aLRT
11. **Homoplasy analysis** — Fitch parsimony CI/RI/HI, delta score, per-taxon convergent change classification

Stages 5-9 run in parallel for dedup and nodedup analysis lines.

## Output structure

```
results/
  {sample}/
    host_filtered/          # Host-removed FASTQ files
    reference/              # Renamed and indexed reference
    denovo_seed/            # (if enabled) Seeded reference
    {line}/                 # dedup/ or nodedup/
      iterative/            # Iterative mapping rounds and final consensus
        round_1/            # Per-round BAM, VCF, consensus
        round_2/
        ...
        final_consensus.fasta
        analysis.bam
        convergence_stats.tsv
      assembly_reads/       # Reads extracted for de novo assembly
      trinity/              # Trinity assembly
      spades/               # SPAdes assembly
      coronaspades/         # coronaSPAdes assembly
      megahit/              # MEGAHIT assembly
      iva/                  # IVA assembly
      ivar/                 # iVar consensus
      final/                # Final mapping to iVar consensus
      {assembler}_polish/   # Assembly-guided polished consensus
      final_alignment/      # MAFFT alignment of all stages
      homoplasy/            # Delta score and pairwise distances
      summary.tsv           # Per-segment QC summary
  iqtree/
    {line}/{segment}/       # IQ-TREE output per line and segment
    iqtree_summary.tsv      # Global phylogenetic metrics
    taxon_homoplasy.tsv     # Per-taxon homoplasy report
```

## Key output files

| File | Description |
|------|-------------|
| `{sample}/{line}/ivar/final_consensus.fasta` | Primary consensus sequence |
| `{sample}/{line}/{asm}_polish/{asm}_guided_polished.fasta` | Assembly-polished consensus |
| `{sample}/{line}/summary.tsv` | Per-segment coverage, depth, completeness, N counts |
| `{sample}/{line}/iterative/convergence_stats.tsv` | Mapped reads and completeness per iteration |
| `{sample}/{line}/homoplasy/homoplasy_report.tsv` | Delta score and site statistics |
| `iqtree/{line}/{segment}/iqtree.treefile` | ML phylogenetic tree (Newick) |
| `iqtree/iqtree_summary.tsv` | Model selection, tree statistics, CI/RI/HI |
| `iqtree/taxon_homoplasy.tsv` | Per-taxon terminal changes and convergence classification |

## Optional features

### De novo reference seeding

Set `denovo_seed: true` to run a quick MEGAHIT assembly before iterative mapping. The de novo contigs are scaffolded against the reference to produce a sample-specific starting reference. This helps when samples diverge substantially from the available reference.

### Reference elongation

Set `extend_reference: true` to extend the consensus beyond original reference boundaries during iterative mapping. Soft-clipped reads at contig termini are used to build consensus extensions, recovering terminal regions not present in the starting reference.

### Automatic reference selection

Set `auto_select_refs: true` and provide reference FASTA databases in `genbank_dir` (one per segment). The pipeline maps reads against all candidates and selects the best reference per segment based on completeness and mapped read count.

## License

[To be specified]
