# ViCon: an iterative consensus pipeline for segmented viral genomes from Illumina sequencing data

Rebecca Lantto<sup>1</sup>, Pär Larsson<sup>1</sup>, Sobiah Rauf<sup>1</sup>, Anne Tuiskunen<sup>1</sup>

<sup>1</sup> Department of Clinical Microbiology, Umeå University, Umeå, Sweden

Correspondence: Anne Tuiskunen

---

## Abstract

Generating accurate consensus genomes from segmented RNA viruses remains challenging due to high mutation rates, host contamination, uneven coverage, and the absence of closely related reference genomes. Here we present ViCon (Viral Consensus), a Snakemake pipeline that combines iterative reference-guided mapping, multiple de novo assemblers, and assembly-guided polishing to produce high-quality consensus sequences for segmented viruses. The pipeline processes paired-end Illumina reads through host removal, convergence-based iterative mapping, iVar consensus calling, and polishing with up to five independent de novo assemblies (Trinity, SPAdes, coronaSPAdes, MEGAHIT, IVA). To reduce reference bias for divergent samples, ViCon optionally constructs a sample-specific starting reference via de novo assembly prior to iterative mapping, and can extend consensus sequences beyond original reference boundaries using soft-clipped reads. To assess and optimise sequence reliability, ViCon places each consensus method individually onto a backbone phylogeny constructed from GenBank references using IQ-TREE, performs site-level Fitch parsimony to quantify homoplasy, and automatically selects the best consensus per segment based on the fewest convergent changes. A phylogenetically-informed hybrid consensus is constructed by preferring non-homoplastic bases at positions where methods disagree. The selected consensus is validated by re-placement on the backbone tree, and allele frequencies from the sequencing reads are checked at disputed positions to confirm read-level support. Per-sample and cross-sample HTML reports with coverage plots, ranking tables, and heatmaps consolidate all quality metrics. The pipeline runs parallel dedup and nodedup analysis lines, supports automatic reference selection from local reference sequence databases, and enforces reading-frame integrity through out-of-frame indel and internal stop codon filtering. ViCon is implemented in Snakemake with Pixi-managed dependencies and is freely available at [URL].

**Keywords:** viral genomics, consensus genome, iterative mapping, segmented virus, hantavirus, Snakemake, phylogenetics

---

## Introduction

Segmented RNA viruses, including hantaviruses (Bunyavirales), influenza viruses, and rotaviruses, pose significant challenges for whole-genome sequencing and consensus generation. Their segmented genomes require simultaneous assembly of multiple genomic segments (e.g. L, M, and S for hantaviruses), each of which may have different coverage profiles and evolutionary rates (Jonsson et al. 2010). Furthermore, clinical and environmental samples often contain substantial host-derived reads that must be removed prior to analysis, and segments can undergo reassortment, complicating phylogenetic interpretation (Klempa 2018).

Reference-guided mapping is the most common approach for generating viral consensus genomes, but its accuracy depends heavily on the similarity between the reference and the sample (Wymant et al. 2018). For divergent viral lineages, a single round of mapping to a distant reference can introduce reference bias, where the consensus inherits bases from the reference rather than reflecting the true sample sequence. This problem is particularly acute for highly variable RNA viruses (Fedonin et al. 2019). Iterative mapping, in which the consensus from one round serves as the reference for the next, can reduce this bias (Shepard et al. 2016; Wymant et al. 2018) but requires a reliable convergence criterion to avoid overrefinement or premature termination.

De novo assembly provides a reference-free alternative but is complicated by uneven coverage, short read lengths relative to repeat regions, and the difficulty of resolving multiple genomic segments from a single assembly graph (Yang et al. 2012). In practice, no single assembler consistently outperforms all others across diverse viral datasets (Ponten et al. 2024; Gruning et al. 2019), suggesting that integrating results from multiple assemblers may improve consensus accuracy (Wan et al. 2015).

Several existing tools address parts of this problem. IRMA (Iterative Refinement Meta-Assembler) pioneered adaptive iterative mapping for influenza and other segmented viruses (Shepard et al. 2016). Shiver demonstrated the importance of sample-specific reference construction for reducing bias in HIV genome assembly (Wymant et al. 2018). VirGenA combines iterative mapping with de novo reassembly of divergent regions (Fedonin et al. 2019). iVar provides primer-aware consensus calling from mapped reads (Grubaugh et al. 2019). V-pipe offers a modular Snakemake-based pipeline for viral variant calling and consensus generation (Posada-Cespedes et al. 2021). VirAmp integrates multiple assemblers in a Galaxy-based workflow (Wan et al. 2015). nf-core/viralrecon (Patel et al. 2023) is a widely adopted Nextflow pipeline supporting both amplicon and metagenomics protocols with single-pass mapping and optional de novo assembly. SIGNAL (Nasir et al. 2024) and CoVpipe2 (Tsai et al. 2024) provide Snakemake- and Nextflow-based pipelines tailored for SARS-CoV-2. More recently, PVGA introduced an iterative alignment graph approach that unifies assembly and polishing (Song et al. 2025). However, few tools combine iterative mapping, multi-assembler de novo assembly, assembly-guided polishing, and phylogenetic quality assessment in a single reproducible workflow optimised for segmented viruses.

Here we present ViCon (Viral Consensus), a Snakemake pipeline designed for end-to-end processing of segmented viral genomes from paired-end Illumina sequencing data. ViCon integrates host removal, automatic reference selection, optional de novo reference seeding, convergence-based iterative mapping with optional reference elongation, five independent de novo assemblers, iVar consensus calling, assembly-guided polishing, phylogenetic analysis with homoplasy metrics, phylogeny-guided best consensus selection, hybrid consensus construction, validation, read support analysis, and comprehensive HTML reporting. The pipeline is designed for hantavirus surveillance but is applicable to other segmented viruses with minimal configuration changes.

---

## Methods

### Pipeline overview

ViCon is implemented as a Snakemake (Mölder et al. 2021) workflow with dependencies managed by Pixi (prefix.dev 2024), a modern package manager built on the conda ecosystem that provides automatic lockfile-based reproducibility. Packages are sourced from the bioconda (Grüning et al. 2018) and conda-forge channels. The pipeline accepts paired-end Illumina FASTQ files and a per-sample multi-segment reference FASTA (or uses automatic reference selection), and produces consensus sequences, phylogenetic trees, and quality metrics. Figure 1 shows the overall workflow.

The pipeline processes each sample through the following stages: (1) host read removal, (2) optional automatic reference selection, (3) reference preparation and segment labelling, (4) optional de novo reference seeding, (5) convergence-based iterative mapping with optional reference elongation, (6) de novo assembly with multiple assemblers, (7) iVar per-contig consensus calling, (8) assembly-guided polishing, (9) multi-sequence alignment of all pipeline stages, (10) phylogenetic backbone tree construction and individual consensus placement, (11) homoplasy analysis, (12) best consensus selection and hybrid consensus construction, (13) consensus validation by re-placement on the backbone tree, (14) read support analysis at disputed sites, and (15) per-sample and cross-sample HTML reporting. Stages 5–14 are run in parallel for two analysis lines: one with PCR duplicate removal ("dedup") and one without ("nodedup"), enabling direct comparison of deduplication effects on consensus accuracy (Ebbert et al. 2016).

### Host read removal

Paired-end reads are filtered against host reference genomes using Bowtie2 (Langmead and Salzberg 2012) in `--very-sensitive` mode. For hantavirus samples from bank voles (*Myodes glareolus*), reads are sequentially mapped against a human genome index and a bank vole genome index. Unmapped read pairs from both filtering steps are retained for downstream analysis. This sequential approach ensures removal of both human contamination (from sample handling) and host reads, a critical preprocessing step for metagenomic and clinical viral samples (Matranga et al. 2014).

### Automatic reference selection

When a closely related reference is not available a priori, ViCon can automatically select the best reference for each genomic segment. This is particularly important given the substantial genetic diversity among hantavirus lineages, even within a single geographic region (Lundkvist et al. 2008; Johansson et al. 2024). Host-filtered reads are mapped with Minimap2 (Li 2018) against a user-provided database of reference sequences (e.g. curated from GenBank) for each segment (L, M, S). For each candidate reference, ViCon computes the fraction of positions covered at >2x depth (completeness) and the total number of mapped reads. References are ranked by completeness (primary) and mapped read count (secondary), and the best reference per segment is selected. The per-segment references are combined into a single multi-segment FASTA for downstream analysis.

### Reference preparation

Reference contigs are labelled by size as L (largest), M (middle), and S (smallest), following the standard naming convention for hantavirus genome segments (Jonsson et al. 2010). The reference is indexed for the selected aligner (BWA-MEM2 or strobealign) and samtools.

### De novo reference seeding (optional)

When the available reference genome is distantly related to the sample, a single reference can introduce substantial mapping bias even with iterative refinement. ViCon provides an optional de novo seeding step, inspired by the "assemble first" approach of shiver (Wymant et al. 2018) and viral-ngs (Park et al. 2015), which constructs a sample-specific starting reference before iterative mapping begins.

When enabled (`denovo_seed: true`), host-filtered reads are assembled with MEGAHIT (Li et al. 2015), chosen for its speed and low memory footprint. The resulting contigs are aligned to the original reference using Minimap2 in assembly-to-reference mode (`-ax asm5`). Variants identified in this alignment are called with bcftools and applied to the reference, producing a "seeded" reference that incorporates sample-specific variation while preserving the coordinate system of the original reference. This seeded reference then serves as the starting point for iterative mapping, typically accelerating convergence and improving consensus accuracy for divergent samples.

### Iterative mapping

ViCon performs convergence-based iterative mapping to progressively improve the consensus sequence and reduce reference bias. This approach is conceptually similar to the iterative strategies employed by IRMA (Shepard et al. 2016), shiver (Wymant et al. 2018), and VirGenA (Fedonin et al. 2019), but uses a dual-metric convergence criterion tailored for segmented viral genomes. In each iteration:

1. Reads are aligned to the current reference using BWA-MEM2 (Vasimuddin et al. 2019) or strobealign (Sahlin 2022), sorted, and assigned read groups.
2. For the dedup line, PCR duplicates are removed using Picard MarkDuplicates (Broad Institute 2019).
3. Mapping statistics are computed: total mapped reads and genome completeness (fraction of reference positions with depth >=2x).
4. Variants are called using bcftools mpileup and bcftools call (Danecek et al. 2021) with haploid ploidy, normalised, and used to generate a new consensus.
5. Optionally, out-of-frame indels are filtered (indels whose length is not a multiple of 3), and variants that would introduce internal stop codons are removed using a custom Python filter.

The iteration continues until both mapped reads and completeness stop improving (i.e. neither metric increases compared to the previous round), subject to a configurable minimum (default: 2) and maximum (default: 10) number of iterations. This dual-metric criterion prevents premature convergence when only one metric has plateaued while the other continues to improve. Convergence statistics (mapped reads and completeness per round) are recorded for quality assessment.

**Reference elongation (optional).** When enabled (`extend_reference: true`), ViCon extends the consensus beyond the original reference boundaries after each iteration round, inspired by IRMA's reference elongation strategy (Shepard et al. 2016). Reads mapped near contig boundaries frequently contain soft-clipped bases representing genuine viral sequence not captured by the reference. A custom Python script (`extend_consensus.py`) identifies reads with soft clips of ≥20 bp near the 5' and 3' ends of each contig, builds a majority-rule consensus extension (requiring ≥3 supporting reads and ≥60% base agreement), and appends the extension to the contig. This enables recovery of terminal regions that may be missing from the starting reference, progressively extending coverage across iterations.

### Reading frame integrity filters

ViCon includes two filters to preserve open reading frames in the consensus, ensuring that the resulting sequences are biologically plausible for downstream functional analyses:

**Out-of-frame indel filter.** After variant calling, bcftools filter removes indels whose net length change is not a multiple of 3, preventing frameshift mutations that would disrupt the coding sequence.

**Internal stop codon filter.** A custom Python script (`filter_stop_variants.py`) evaluates each remaining variant in the context of the full contig sequence. If applying a variant would introduce a premature stop codon (TAA, TAG, or TGA) that was not present in the reference, the variant is excluded. This filter assumes that each contig represents a single coding sequence starting at position 1, which is the standard architecture for hantavirus segments (Jonsson et al. 2010).

### De novo assembly

To provide reference-independent sequence information for polishing, ViCon runs up to five de novo assemblers on reads extracted from the final iterative-mapping BAM. The use of multiple assemblers follows the rationale that no single assembler consistently performs best across all viral datasets (Ponten et al. 2024; Sutton et al. 2019), and their complementary algorithmic strategies can capture different aspects of the viral genome:

- **Trinity** (Grabherr et al. 2011): RNA-Seq assembler designed for transcript reconstruction, run via Docker container
- **SPAdes** (Bankevich et al. 2012): in `--rnaviral` mode, leveraging RNA viral-specific heuristics
- **coronaSPAdes** (Meleshko et al. 2022): optimised for single-stranded RNA viral genomes
- **MEGAHIT** (Li et al. 2015): memory-efficient de Bruijn graph assembler suitable for large datasets
- **IVA** (Hunt et al. 2015): iterative virus assembler specifically designed for viral genome assembly from short reads

Each assembler is configured independently with user-specified memory limits and thread counts.

### iVar consensus calling

The original host-filtered reads are remapped to the final iterative-mapping consensus using the selected aligner. Per-contig consensus sequences are generated using iVar (Grubaugh et al. 2019) with a minimum depth threshold of 2x. Positions below this threshold are called as N. This step produces the primary consensus output for each analysis line.

### Assembly-guided polishing

Each de novo assembly is aligned against the iVar consensus using Minimap2 in assembly-to-reference mode (`-ax asm5`). Variants identified in this alignment are called with bcftools and applied to the consensus, producing an assembly-guided polished sequence for each assembler. This approach is conceptually analogous to long-read polishing strategies (Zimin and Salzberg 2020) but uses short-read de novo assemblies as the correction source. This step can correct errors in the mapping-based consensus by incorporating information from reference-free assembly, particularly in regions of low coverage or high divergence from the original reference.

### Assembly-to-consensus comparison

Each de novo assembly is also compared against the iVar consensus to generate PAF alignment files and VCF variant files. These outputs allow users to assess concordance between mapping-based and assembly-based approaches and identify regions of disagreement, providing an additional quality control layer.

### Multi-sequence alignment

For each sample and analysis line, all consensus sequences from every pipeline stage are aligned with MAFFT (Katoh and Standley 2013) in `--auto` mode. The alignment includes: the original reference, each iterative-mapping round consensus, the iVar consensus, and all assembly-polished sequences. This per-segment alignment enables visual inspection and quantitative comparison of how the consensus evolves through the pipeline.

### Phylogenetic analysis

ViCon constructs maximum-likelihood phylogenetic trees per segment using IQ-TREE (Minh et al. 2020) in a two-phase approach. Building segment-specific trees is important for detecting reassortment, which is well documented in hantaviruses (Klempa 2018; Kim et al. 2016).

**Backbone tree construction.** For each segment, GenBank reference sequences (optionally subsampled to a configurable maximum) are aligned with MAFFT (Katoh and Standley 2013) and used to construct a backbone maximum-likelihood tree with IQ-TREE using ModelFinder Plus (MFP; Kalyaanamoorthy et al. 2017) for automatic model selection, 1000 ultrafast bootstrap replicates (Hoang et al. 2018), and 1000 SH-aLRT replicates. This backbone tree provides a stable phylogenetic framework for subsequent sample placement.

**Individual consensus placement.** Each consensus method generated by the pipeline (reference, iterative rounds, iVar, assembly-polished variants) is placed individually onto the backbone tree. For each consensus, MAFFT `--add --keeplength` appends the sequence to the backbone alignment without altering existing columns, and IQ-TREE is run with the backbone tree as a topological constraint (`-g`), optimising only the placement of the new taxon. This produces a separate tree and alignment for each consensus method per segment per sample, enabling direct comparison of phylogenetic behaviour across methods. A manifest file (`consensus_placements.tsv`) records the method origin, taxon name, and output paths for each placement.

### Homoplasy analysis

ViCon computes multiple metrics to assess phylogenetic signal quality and potential artifacts:

**Per-sample metrics.** For each sample's multi-stage alignment, ViCon computes: variable sites, parsimony-informative sites, pairwise Hamming distances between all consensus stages, and the delta score (Holland et al. 2002) as a measure of treelikeness. The delta score quantifies the degree to which pairwise distance data conform to a tree-like pattern using quartet analysis: values near 0 indicate perfectly tree-like signal, while values approaching 1 indicate substantial phylogenetic conflict, potentially arising from recombination, reassortment, or systematic methodological artifacts.

**Global phylogenetic metrics.** From the IQ-TREE phylogeny, Fitch parsimony (Fitch 1971) is used to compute the consistency index (CI; Kluge and Farris 1969), retention index (RI; Farris 1989), and homoplasy index (HI = 1 - CI). These metrics quantify the extent of homoplasy in the dataset: CI measures the ratio of minimum possible character changes to observed changes, RI accounts for the maximum possible homoplasy, and HI directly quantifies the proportion of observed change attributable to homoplasy. Low CI values (below 0.5) suggest substantial homoplasy (Crispell et al. 2019), which may arise from convergent evolution, recombination, or sequencing artifacts.

**Per-taxon homoplasy.** For each terminal taxon (sample consensus), ViCon identifies terminal branch changes via Fitch top-down state assignment and classifies them as:
- **Autapomorphic**: unique to the taxon (no other taxon shares the derived state)
- **Convergent**: shared with at least one other taxon independently (homoplasy)

Convergent changes are further classified by the patristic distance to the nearest taxon sharing the same derived state, relative to the median pairwise distance in the tree:
- **Near-convergent**: nearest sharer is below the median distance (potentially due to local recombination or shared methodological bias)
- **Distant-convergent**: nearest sharer is above the median distance (consistent with true convergent evolution or systematic error)

GenBank reference sequences are aggregated and reported as a mean for comparison against sample-derived sequences.

**Individual placement homoplasy.** In addition to the global analysis, each individual consensus placement (from the per-method backbone placement step) is analysed separately with per-site Fitch parsimony. This produces a comprehensive table (`individual_placement_homoplasy.tsv`) reporting tree statistics (log-likelihood, tree length, parsimony score), CI/RI/HI, and per-taxon change classification for every combination of sample, segment, and consensus method.

### Best consensus selection

ViCon automatically selects the best consensus method per segment based on phylogenetic quality metrics, rather than relying on a single predetermined method.

**Site-level Fitch analysis.** For each segment and each consensus method, site-level Fitch parsimony is performed on the individual placement tree. At each variable alignment column, the bottom-up pass identifies whether a change occurred on the sample's terminal branch, and the top-down state assignment classifies each change as convergent (homoplastic) or autapomorphic. The results are recorded in `site_homoplasies.tsv` with the alignment position, sample state, ancestral state, and classification for every flagged site.

**Method ranking.** Consensus methods are ranked per segment by a composite criterion: fewest convergent (homoplastic) sites (primary), fewest autapomorphic sites (secondary), and fewest total changes (tertiary). This ranking prioritises phylogenetic consistency — the method whose consensus introduces the least homoplasy when placed on the backbone tree. The full ranking is recorded in `consensus_ranking.tsv`.

**Best consensus.** The rank-1 method per segment is selected as the best consensus. The corresponding sequences are extracted from the individual placement results and written to `best_consensus.fasta` as a multi-segment FASTA.

**Hybrid consensus.** A phylogenetically-informed hybrid consensus is constructed by examining each alignment position across all consensus methods. At positions where all methods agree, the shared base is used. At positions where methods disagree, ViCon preferentially selects the base that is not associated with homoplasy in any method. If multiple non-homoplastic options exist, majority rule is applied. If all options are homoplastic, the base from the rank-1 method is used as a fallback. This approach can resolve individual homoplastic positions that persist in even the best single-method consensus, producing a hybrid sequence with potentially lower overall homoplasy.

### Consensus validation

To verify that the best and hybrid consensus sequences are phylogenetically sound, ViCon places them back onto the backbone tree using the same MAFFT + IQ-TREE constraint approach used for individual placements. Site-level Fitch analysis is then repeated on the resulting trees. The validation output (`validation.tsv`) reports homoplastic sites, autapomorphic sites, total changes, and total variable sites for each of: the best consensus, the hybrid consensus, and the original rank-1 method. This enables direct comparison — a well-constructed best or hybrid consensus should show equal or fewer homoplastic sites than the original method.

### Read support analysis

At positions flagged as homoplastic by the site-level Fitch analysis, ViCon checks whether the sequencing reads support the consensus call. Alignment column positions from `site_homoplasies.tsv` are mapped to ungapped contig coordinates by counting non-gap characters in the sample's aligned sequence. For each flagged position, allele counts (A, C, G, T, deletions) are extracted from the final BAM file using samtools mpileup (Danecek et al. 2021). The consensus allele frequency is computed as the fraction of reads supporting the called base, and each site receives a verdict:
- **Supported**: >80% of reads agree with the consensus call
- **Weak**: 50–80% of reads agree
- **Contradicted**: <50% of reads agree (reads favour a different base)

Contradicted sites may indicate sequencing errors, low-frequency variants, or systematic biases in the consensus method. The results are recorded in `read_support.tsv`.

### Summary reporting

For each sample and analysis line, ViCon produces a comprehensive TSV summary including: per-segment reference length, mapped read counts, mean depth, and genome completeness at each mapping stage (iterative rounds, iVar, final); N counts and N fractions for each consensus stage; and internal stop codon flags.

### HTML reporting

**Per-sample reports.** For each sample and analysis line, ViCon generates an interactive HTML report (`report.html`) consolidating all results into a single document with embedded plots. The report includes: iterative mapping convergence plots (mapped reads and completeness per round), per-position sequencing depth plots for each segment, coverage summary tables, consensus quality metrics (N counts, internal stop codons), the full consensus method ranking with the best method highlighted, site-level homoplasy details, phylogenetic placement summaries, consensus validation comparisons with improvement/degradation badges, and read support analysis with verdict summaries. Depth plots are generated using samtools depth and rendered as area plots with mean depth annotations.

**Cross-sample overview.** A global HTML report (`overview_report.html`) compares results across all samples with: coverage completeness heatmaps (samples × segments), mean depth heatmaps, a bar chart of best-method frequency showing which consensus methods are selected most often, homoplasy heatmaps showing the number of homoplastic sites per sample and segment for the best method, and a validation summary table. These visualisations enable rapid identification of problematic samples or segments across large surveillance datasets.

### Implementation and configuration

ViCon is implemented as a single Snakemake workflow with a YAML configuration file. Key configurable parameters include:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `aligner` | strobealign | Read aligner (bwa-mem2 or strobealign) |
| `max_iterations` | 10 | Maximum iterative mapping rounds |
| `min_iterations` | 2 | Minimum rounds before convergence check |
| `lines` | [dedup, nodedup] | Analysis lines to run |
| `enforce_inframe` | true | Filter out-of-frame indels |
| `check_stops` | true | Check for internal stop codons |
| `assemblers` | [trinity, spades, coronaspades, megahit, iva] | De novo assemblers |
| `denovo_seed` | false | De novo assembly to seed iterative mapping |
| `extend_reference` | false | Extend consensus using soft-clipped reads |
| `auto_select_refs` | false | Automatic reference selection |
| `threads` | 8 | Threads per job |

Samples are auto-detected from paired-end FASTQ files in the input directory, or a subset can be specified explicitly. Dependencies are managed by Pixi with packages from bioconda (Grüning et al. 2018) and conda-forge, ensuring reproducibility through lockfile-based dependency resolution.

---

## Results

### Application to hantavirus surveillance

ViCon was developed for and applied to Puumala hantavirus (PUUV) surveillance samples collected from bank voles (*Myodes glareolus*) in Sweden. PUUV causes nephropathia epidemica, the most common hantavirus disease in Europe, and exhibits substantial genetic diversity with multiple co-circulating lineages in Scandinavia (Lundkvist et al. 2008; Johansson et al. 2024). The pipeline successfully processed samples through all stages, generating consensus sequences for the L (~6.5 kb), M (~3.6 kb), and S (~1.8 kb) genomic segments.

The iterative mapping stage typically converged within 2–4 rounds, with mapped read counts and genome completeness stabilising rapidly. The convergence criterion (requiring simultaneous improvement in both metrics) prevented unnecessary iterations while ensuring adequate refinement of the consensus.

Running five de novo assemblers in parallel provided multiple independent assemblies for polishing. While individual assemblers occasionally failed on low-coverage samples or produced fragmented assemblies, the multi-assembler approach ensured that at least one assembly was available for polishing in most cases. Assembly-guided polishing corrected a small number of positions in the mapping-based consensus, primarily in regions of low coverage or high divergence from the initial reference.

The parallel dedup/nodedup analysis lines revealed that PCR duplicate removal had a measurable but generally modest effect on consensus accuracy for samples with adequate sequencing depth, consistent with previous findings on the effect of duplicate removal in deep sequencing data (Ebbert et al. 2016). For low-coverage segments, nodedup consensus sequences sometimes achieved higher completeness at the cost of potential PCR duplicate-driven errors.

### Phylogenetic quality assessment

IQ-TREE analysis with ModelFinder selected appropriate substitution models for each segment, with GTR+G variants commonly selected for the more variable L and M segments. The homoplasy metrics (CI, RI, HI) provided a quantitative assessment of phylogenetic signal quality. Per-taxon homoplasy analysis distinguished between autapomorphic changes (unique to individual samples, potentially representing true variation or sequencing errors) and convergent changes (shared across the tree, potentially indicating systematic biases or genuine convergent evolution). The near/distant convergent classification based on patristic distance provided additional resolution for interpreting the source of homoplastic signal.

Segment-specific phylogenetic trees enabled detection of potential reassortment events, a phenomenon well documented in hantaviruses including PUUV (Klempa 2018; Razzauti et al. 2009). Incongruent tree topologies between segments can indicate reassortment, which ViCon facilitates through its per-segment tree construction approach.

### Consensus method comparison and selection

Individual phylogenetic placement of each consensus method onto the backbone tree revealed variation in homoplasy levels across methods. Iterative mapping consensus sequences from early rounds typically showed higher convergent change counts, reflecting residual reference bias, while later rounds, iVar consensus, and assembly-polished sequences generally showed lower homoplasy. The automatic ranking system consistently selected the method with the fewest convergent changes, which varied across segments and samples — no single method dominated universally, confirming the value of evaluating multiple approaches.

The hybrid consensus construction resolved individual homoplastic positions at sites where at least one method produced a non-homoplastic base. Validation by re-placement on the backbone tree confirmed that the best and hybrid consensus sequences maintained equal or reduced homoplasy compared to the original best single-method consensus. Read support analysis at disputed sites provided an independent line of evidence: the majority of flagged positions showed strong read support (>80% consensus allele frequency), while a small fraction of contradicted sites (<50% support) highlighted positions where the consensus call was uncertain and may warrant manual review.

### Reporting

The per-sample HTML reports provided a consolidated view of all pipeline metrics in a single document, replacing the need to inspect multiple TSV files. Coverage depth plots allowed rapid identification of low-coverage regions, while the consensus ranking tables and validation badges enabled immediate assessment of consensus quality. The cross-sample overview report proved particularly useful for surveillance datasets with many samples, enabling rapid identification of samples with unusually high homoplasy, low coverage, or anomalous best-method selection patterns through colour-coded heatmaps.

---

## Discussion

### Comparison to existing pipelines

ViCon occupies a distinct niche among viral genome assembly pipelines by combining iterative mapping, multi-assembler de novo assembly, assembly-guided polishing, and integrated phylogenetic quality assessment in a single workflow optimised for segmented viruses (Table 1).

IRMA (Shepard et al. 2016) is the most directly comparable tool, as it was also designed for segmented viruses (primarily influenza) and employs iterative refinement. However, IRMA focuses on read gathering optimisation and variant phasing rather than multi-assembler polishing and homoplasy assessment. Shiver (Wymant et al. 2018) pioneered sample-specific reference construction to reduce mapping bias for HIV but does not include de novo assembly polishing or phylogenetic quality metrics. VirGenA (Fedonin et al. 2019) combines iterative mapping with de novo reassembly of divergent regions but operates on individual genomes rather than multi-sample surveillance datasets. V-pipe (Posada-Cespedes et al. 2021) provides comprehensive variant calling and haplotype reconstruction but uses a single-pass mapping strategy without iterative refinement. VirAmp (Wan et al. 2015) integrates multiple assemblers in a Galaxy interface but lacks the iterative mapping and phylogenetic assessment components. nf-core/viralrecon (Patel et al. 2023) is perhaps the most widely used viral consensus pipeline, supporting both amplicon and metagenomics protocols with single-pass mapping and multiple de novo assemblers, but does not include iterative refinement, assembly-guided polishing, or phylogenetic quality metrics. SIGNAL (Nasir et al. 2024) and CoVpipe2 (Tsai et al. 2024) are optimised for amplicon-based SARS-CoV-2 sequencing and do not address the specific challenges of segmented viral genomes.

ViCon's unique contribution is the integration of these approaches into a cohesive pipeline with a closed-loop quality assurance system: iterative mapping reduces reference bias, multiple assemblers provide complementary de novo information, assembly-guided polishing combines both approaches, and phylogenetic analysis with per-taxon homoplasy metrics provides built-in quality assessment. Critically, ViCon goes beyond quality reporting to actively use phylogenetic information for consensus selection — automatically choosing the best method per segment, constructing a hybrid consensus from non-homoplastic bases, validating the result by re-placement on the backbone tree, and confirming read-level support at disputed positions. The integrated HTML reporting with coverage plots, ranking tables, and cross-sample heatmaps further distinguishes ViCon from existing pipelines that produce only tabular output files.

### Reference bias mitigation

The iterative mapping approach progressively adapts the reference to the sample, reducing reference bias compared to single-pass mapping. Wymant et al. (2018) demonstrated that mapping to a sample-specific reference recovered a median of 205 additional bases compared to mapping to the closest available reference for HIV genomes. ViCon extends this principle through convergence-based iteration, where the consensus from each round serves as the reference for the next, and the dual-metric convergence criterion provides an objective stopping point that balances accuracy against computational cost.

Two optional features further reduce reference bias. De novo reference seeding constructs a sample-specific starting reference before iterative mapping, following the "assemble first" strategy of shiver (Wymant et al. 2018) and viral-ngs (Park et al. 2015). By incorporating sample-derived variants into the starting reference, the first mapping round already achieves higher sensitivity and specificity, accelerating convergence. Reference elongation, inspired by IRMA (Shepard et al. 2016), extends the consensus beyond the original reference boundaries using soft-clipped reads at contig termini. This is particularly valuable when the reference is shorter than the sample genome, a situation that arises frequently when reference genomes represent partial sequences or when the sample lineage has terminal insertions.

### Multi-assembler consensus

By running five independent de novo assemblers and using their outputs for polishing, ViCon leverages the complementary strengths of different assembly algorithms. Ponten et al. (2024) demonstrated that assembler choice significantly impacts viral genome quality, and that different assemblers perform optimally under different conditions of coverage depth, genome complexity, and read quality. This is particularly valuable for segmented viruses, where different segments may favour different assemblers depending on coverage depth and complexity. The VirAmp pipeline (Wan et al. 2015) similarly recognised the value of combining multiple assembly strategies, integrating both de novo and reference-guided approaches.

### Reading frame integrity

The dual filtering strategy (out-of-frame indels and internal stop codons) preserves open reading frames in the consensus, which is critical for downstream analyses of protein-coding genes. This is especially important for hantaviruses, where each segment encodes a single major open reading frame: L encodes the RNA-dependent RNA polymerase, M encodes the glycoprotein precursor, and S encodes the nucleocapsid protein (Jonsson et al. 2010).

### Integrated quality assessment

The combination of per-stage summary statistics, multi-sequence alignments, phylogenetic trees, and homoplasy metrics provides a comprehensive view of consensus quality. The per-taxon homoplasy analysis can identify samples with unusually high homoplastic signal, which may indicate sequencing artifacts, contamination, or recombination. The consistency and retention indices (Kluge and Farris 1969; Farris 1989) are classical measures of phylogenetic data quality, while the delta score (Holland et al. 2002) provides an alignment-based measure of treelikeness that can detect non-tree-like signal arising from recombination or reassortment — processes of particular relevance for segmented viruses (Klempa 2018). Crispell et al. (2019) demonstrated the utility of identifying homoplasies on phylogenies for downstream interpretation, an approach ViCon extends with per-taxon distance-aware convergence classification.

### Phylogeny-guided consensus selection

Rather than relying on a single consensus method, ViCon's phylogeny-guided selection evaluates all available consensus sequences in a common phylogenetic framework and selects the one that introduces the least homoplasy. This approach is motivated by the observation that different consensus methods may perform differently across samples and segments depending on coverage depth, reference divergence, and assembly quality (Ponten et al. 2024). By ranking methods based on convergent changes — which represent phylogenetically inconsistent signal — the pipeline objectively identifies the most trustworthy consensus for each segment.

The hybrid consensus construction extends this principle to the site level. At positions where methods disagree, phylogenetic information is used to prefer bases that are not associated with homoplasy in any method. This targeted correction avoids the wholesale replacement of one consensus with another and instead resolves only the specific positions where phylogenetic evidence suggests an error. The validation step provides a closed-loop check: by re-placing the best and hybrid consensus on the backbone tree, ViCon confirms that the selection and hybridisation process did not inadvertently introduce new homoplastic signal.

The read support analysis adds an independent, non-phylogenetic line of evidence. Sequencing read allele frequencies at disputed positions can distinguish between genuine variants (strong read support) and potential errors or low-frequency contaminants (weak or contradicted support). This multi-layered approach — phylogenetic ranking, hybrid construction, validation, and read support — provides a robust framework for consensus quality assessment that goes beyond the capabilities of any single metric.

### Reproducibility

The use of Snakemake for workflow management (Mölder et al. 2021) and Pixi for dependency management (prefix.dev 2024) ensures that the pipeline can be reproduced across different computing environments. Pixi provides lockfile-based dependency resolution using the conda ecosystem, including bioconda (Grüning et al. 2018), ensuring exact reproducibility of the software environment. The YAML configuration file makes it straightforward to adapt the pipeline to different viral targets.

### Limitations

ViCon assumes that each genomic segment contains a single open reading frame starting at position 1, which is valid for hantavirus segments but may require adaptation for viruses with overlapping reading frames or non-standard gene structures. The pipeline currently supports paired-end Illumina data only; adaptation for long-read technologies (Oxford Nanopore, PacBio) would require modifications to the mapping and assembly stages. The Trinity assembler is run via Docker, which adds a container dependency that may not be available in all computing environments. The phylogeny-guided consensus selection relies on a representative set of GenBank reference sequences for backbone tree construction; performance may be reduced for novel or poorly sampled viral lineages where the backbone tree does not adequately represent the phylogenetic neighbourhood of the sample. The hybrid consensus construction assumes that disagreements between methods at individual sites can be resolved independently, which may not hold for positions involved in epistatic interactions.

### Future directions

Planned extensions include: support for long-read sequencing data (Oxford Nanopore, PacBio), integration of primer trimming for amplicon-based protocols, and automated recombination and reassortment detection through formal phylogenetic incongruence testing between segment trees. The read support analysis could be extended to incorporate strand bias and mapping quality metrics for more refined error detection.

---

## Availability

ViCon is implemented in Python and Bash, orchestrated by Snakemake, with dependencies managed by Pixi. The source code is available at [URL] under an open-source licence.

---

## References

Bankevich A, Nurk S, Antipov D, et al. (2012) SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 19(5):455–477.

Broad Institute (2019) Picard toolkit. https://broadinstitute.github.io/picard/

Crispell J, Balaz D, Gordon SV (2019) HomoplasyFinder: a simple tool to identify homoplasies on a phylogeny. *Microbial Genomics* 5(1):e000245.

Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10(2):giab008.

Ebbert MTW, Wadsworth ME, Staley LA, et al. (2016) Evaluating the necessity of PCR duplicate removal from next-generation sequencing data and a comparison of approaches. *BMC Bioinformatics* 17(Suppl 7):239.

Farris JS (1989) The retention index and the rescaled consistency index. *Cladistics* 5(4):417–419.

Fedonin GG, Fantin YS, Favorov AV, Shipulin GA, Neverov AD (2019) VirGenA: a reference-based assembler for variable viral genomes. *Briefings in Bioinformatics* 20(1):15–27.

Fitch WM (1971) Toward defining the course of evolution: minimum change for a specific tree topology. *Systematic Zoology* 20(4):406–416.

Grabherr MG, Haas BJ, Yassour M, et al. (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nature Biotechnology* 29(7):644–652.

Grubaugh ND, Gangavarapu K, Quick J, et al. (2019) An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. *Genome Biology* 20(1):8.

Grüning B, Dale R, Sjödin A, et al. (2018) Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nature Methods* 15(7):475–476.

Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018) UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35(2):518–522.

Holland BR, Huber KT, Dress A, Moulton V (2002) Delta plots: a tool for analyzing phylogenetic distance data. *Molecular Biology and Evolution* 19(12):2051–2059.

Hunt M, Gall A, Ong SH, et al. (2015) IVA: accurate de novo assembly of RNA virus genomes. *Bioinformatics* 31(14):2374–2376.

Johansson P, Lagerqvist N, Lindgren P-E, Lundkvist Å (2024) Nephropathia epidemica caused by Puumala virus in bank voles, Scania, southern Sweden. *Emerging Infectious Diseases* 30(4):795–798.

Jonsson CB, Figueiredo LTM, Vapalahti O (2010) A global perspective on hantavirus ecology, epidemiology, and disease. *Clinical Microbiology Reviews* 23(2):412–441.

Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14(6):587–589.

Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4):772–780.

Kim W-K, No JS, Lee S-H, et al. (2016) Genetic diversity and reassortment of Hantaan virus tripartite RNA genomes in nature, the Republic of Korea. *PLOS Neglected Tropical Diseases* 10(6):e0004650.

Klempa B (2018) Reassortment events in the evolution of hantaviruses. *Virus Genes* 54(5):638–646.

Kluge AG, Farris JS (1969) Quantitative phyletics and the evolution of anurans. *Systematic Zoology* 18(1):1–32.

Langmead B, Salzberg SL (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods* 9(4):357–359.

Li D, Liu CM, Luo R, Sadakane K, Lam TW (2015) MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31(10):1674–1676.

Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094–3100.

Lundkvist Å, Hukic M, Hörling J, et al. (2008) Puumala hantavirus genetic variability in an endemic region (Northern Sweden). *Infection, Genetics and Evolution* 8(3):286–296.

Matranga CB, Andersen KG, Winnicki S, et al. (2014) Enhanced methods for unbiased deep sequencing of Lassa and Ebola RNA viruses from clinical and biological samples. *Genome Biology* 15(11):519.

Meleshko D, Hajirasouliha I, Bankevich A (2022) coronaSPAdes: from biosynthetic gene clusters to RNA viral assemblies. *Bioinformatics* 38(1):1–8.

Minh BQ, Schmidt HA, Chernomor O, et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37(5):1530–1534.

Mölder F, Jablonski KP, Letcher B, et al. (2021) Sustainable data analysis with Snakemake. *F1000Research* 10:33.

Nasir JA, Kozak RA, Aftanas P, et al. (2024) SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL), a Snakemake workflow for rapid and bulk analysis of Illumina sequencing of SARS-CoV-2 genomes. *NAR Genomics and Bioinformatics* 6(4):lqae176.

Patel H, Varona S, Monzón S, et al. (2023) nf-core/viralrecon: nf-core/viralrecon v2.6.0 — Rhodium Raccoon. *Zenodo*. https://doi.org/10.5281/zenodo.3901628

Park DJ, Dudas G, Wohl S, et al. (2015) Ebola virus epidemiology, transmission, and evolution during seven months in Sierra Leone. *Cell* 161(7):1516–1526.

Ponten TS, Smaling R, Geurts van Kessel CH, Koenig J, Baaijens JA (2024) Comparative evaluation of open-source bioinformatics pipelines for full-length viral genome assembly. *Viruses* 16(12):1824.

Posada-Cespedes S, Seez D, Aquino Y, et al. (2021) V-pipe: a computational pipeline for assessing viral genetic diversity from high-throughput data. *Bioinformatics* 37(11):1673–1680.

prefix.dev (2024) Pixi: package management made easy. https://pixi.sh/

Razzauti M, Plyusnina A, Henttonen H, Plyusnin A (2009) Analysis of Puumala hantavirus in a bank vole population in northern Finland: evidence for co-circulation of two genetic lineages and frequent reassortment between strains. *Journal of General Virology* 90(8):1923–1931.

Sahlin K (2022) Effective sequence similarity detection with strobemers. *Genome Research* 32(7):1381–1394.

Shepard SS, Meno S, Bahl J, et al. (2016) Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler. *BMC Genomics* 17:708.

Song Z, Cai D, Sun Y, Wang L (2025) PVGA: a precise viral genome assembler using an iterative alignment graph. *GigaScience* 14:giaf063.

Sutton TD, Clooney AG, Ryan FJ, Ross RP, Hill C (2019) Choice of assembly software has a critical impact on virome characterisation. *Microbiome* 7:12.

Tsai Y-C, Hölzer M, Tautenhahn N, et al. (2024) Lessons learned: overcoming common challenges in reconstructing the SARS-CoV-2 genome from short-read sequencing data via CoVpipe2. *F1000Research* 12:1091.

Vasimuddin M, Misra S, Li H, Aluru S (2019) Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In *IEEE International Parallel and Distributed Processing Symposium (IPDPS)*, pp. 314–324.

Wan Y, Renner DW, Albert I, Bhatt AS (2015) VirAmp: a galaxy-based viral genome assembly pipeline. *GigaScience* 4:19.

Wymant C, Blanquart F, Golubchik T, et al. (2018) Easy and accurate reconstruction of whole HIV genomes from short-read sequence data with shiver. *Virus Evolution* 4(1):vey007.

Yang X, Charlebois P, Gnerre S, et al. (2012) De novo assembly of highly diverse viral populations. *BMC Genomics* 13:475.

Zimin AV, Salzberg SL (2020) The genome polishing tool POLCA makes fast and accurate corrections in genome assemblies. *PLOS Computational Biology* 16(6):e1007981.
