# Integrated Mouse Methylome Analysis

## Overview

This project analyzes whole-genome bisulfite sequencing (WGBS) data from developing mouse kidney and liver at E14.5 and E15.5.

The workflow covers:

- WGBS quality control and preprocessing
- Bisulfite-aware alignment and methylation calling
- Genome-wide methylation analysis
- Differential methylation between kidney and liver
- Methylation-state segmentation
- ChromHMM state integration
- ATAC-seq and histone-mark integration
- RNA-seq integration
- Gene-level interpretation of selected loci

## Dataset

- **Samples:** 8 single-end WGBS samples
- **Tissues:** Kidney, Liver
- **Developmental stages:** E14.5, E15.5
- **Replicates:** 2 biological replicates per tissue/timepoint
- **Reference:** mm10 / GRCm38 reduced reference
- **Main comparison:** Kidney vs Liver
- **Genomic region:** chr18–chr19

## WGBS Processing

WGBS reads were processed using the nf-core/methylseq workflow.

Main steps:

1. FastQC
2. Adapter trimming with Trim Galore!
3. Bisulfite-aware alignment with Bismark/Bowtie2
4. Deduplication
5. Methylation extraction
6. Bismark reporting
7. Qualimap
8. Preseq
9. MultiQC

## Quality Control

All eight samples passed the FastQC per-base sequence quality check.

All eight samples showed the expected failure in per-base sequence content due to bisulfite conversion.

Alignment efficiency was consistently high across samples, ranging from approximately 97–99%.

Global CpG methylation was higher in kidney (~77–79%) than liver (~60–62%).

## Methylome Analysis

RnBeads was used for:

- Methylation quality control
- Coverage filtering
- CpG annotation
- Exploratory analysis
- Differential methylation

Filtering included:

- Coverage threshold of 10
- SNP-overlapping site removal
- Low-coverage masking
- High-coverage outlier filtering

The final dataset contained approximately 1.16 million retained CpG sites.

## Exploratory Analysis

The methylation profiles were explored using:

- PCA
- PCA variance/scree plots
- Methylation density plots
- CpG-context methylation distributions
- Heatmap of the 1,000 most variable CpGs
- Sample-level boxplots
- CpG and region-level violin plots

The main pattern was a clear separation of samples by tissue, with kidney and liver forming distinct methylation groups.

## Custom Region Analysis

A custom genomic annotation containing 6,145 regions was incorporated into RnBeads.

After removal of regions containing missing values, 6,113 regions were used for exploratory analysis.

The 100 most variable regions were analyzed and showed strong tissue-specific clustering.

## Differential Methylation

Differential methylation was assessed primarily between kidney and liver.

The analysis included:

- Tiling regions
- Genes
- Promoters
- CpG islands
- Custom regions

The strongest differential methylation was observed in CpG islands and promoter-associated regions.

The top custom DMRs were inspected using IGV and compared between kidney and liver.

## Methylation-State Segmentation

MethylSeekR-based segmentation was performed on merged kidney and liver methylomes.

The resulting methylation states were:

- **UMR:** Unmethylated Region
- **LMR:** Low Methylated Region
- **HMR:** Highly Methylated Region
- **PMD:** Partially Methylated Domain

Segmentation tracks were exported as BED files and inspected together with methylation signal in IGV.

## Multi-Omic Integration

DNA methylation was integrated with:

- ChromHMM 15-state chromatin annotations
- ATAC-seq accessibility
- Histone modifications
- RNA-seq expression

General patterns included:

- Low methylation at promoters
- High accessibility at active promoters
- LMRs at enhancer-associated regions
- High methylation across many transcribed gene bodies
- H3K4me3/H3K27ac/H3K9ac at active regulatory regions
- H3K36me3 across actively transcribed gene bodies
- H3K27me3/H3K9me3 in repressive chromatin

## TSS-Centered Analysis

Signals were summarized around transcription start sites using:

- 500 bp upstream of the TSS
- 1,500 bp downstream of the TSS

The analysis was strand-aware so that the same biological window was applied correctly to both plus- and minus-strand genes.

deepTools was used for signal summarization and visualization.

## DMR Integration

Differentially methylated regions were intersected with:

- Tissue-specific ChromHMM heterochromatin states
- Differentially accessible ATAC-seq peaks

Observed overlaps included:

- Kidney heterochromatin: 257 DMRs
- Liver heterochromatin: 676 DMRs
- Kidney accessible peaks: 105 DMRs
- Liver accessible peaks: 499 DMRs

## Gene-Level Examples

### Abcc2

Abcc2 was used as a liver example showing coordinated differences involving:

- Differential expression
- DNA methylation
- ATAC-seq accessibility
- ChromHMM states
- Histone modifications

### Kcnn2

Kcnn2 was used as a kidney example showing:

- Differential expression
- Differential accessibility
- Promoter-associated chromatin states
- Low promoter methylation
- Active promoter histone marks

## Main Findings

- Tissue identity was the strongest source of methylation variation.
- Kidney showed higher global CpG methylation than liver.
- Promoters were generally hypomethylated.
- Gene bodies were generally highly methylated.
- Methylation states showed consistent relationships with ChromHMM states.
- Active regulatory regions were associated with accessibility and active histone marks.
- DMRs overlapped tissue-specific chromatin and accessibility features.
- Integrating methylation with other epigenomic layers provided a broader view of tissue-specific regulation.

## Main Tools

- nf-core/methylseq
- FastQC
- Trim Galore!
- Bismark
- Bowtie2
- Qualimap
- Preseq
- MultiQC
- RnBeads
- MethylSeekR
- IGV
- ChromHMM
- deepTools
- R
- GenomicRanges
- ggplot2

## Reproducibility

The analysis was performed in the COMPEPIWS/de.NBI Cloud environment using Conda, Nextflow/nf-core, R and command-line genomic tools.
The reduced mm10 reference containing chromosomes 18 and 19 was used throughout the WGBS analysis.

## Author

Sakshi Parate

M.Sc. Bioinformatics, Saarland University
