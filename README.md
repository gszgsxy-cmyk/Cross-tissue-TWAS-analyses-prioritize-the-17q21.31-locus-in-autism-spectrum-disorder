# ASD TWAS Analysis Code

This repository contains the analysis scripts and selected reproducibility resources associated with the following publication:

Li X, Zhang Y, Li R, et al. Cross-tissue transcriptome-wide association analyses prioritize the 17q21.31 locus in autism spectrum disorder. *Progress in Neuro-Psychopharmacology & Biological Psychiatry*. [https://doi.org/10.1016/j.pnpbp.2026.112014](https://doi.org/10.1016/j.pnpbp.2026.112014).

## Study overview

This study integrated autism spectrum disorder (ASD) genome-wide association study (GWAS) summary statistics with genetically predicted gene expression across 29 tissues using Joint-Tissue Imputation (JTI) and FUSION. Follow-up analyses included summary-data-based Mendelian randomization (SMR) using brain and non-brain expression and methylation QTL resources, regional FUSION conditional/joint analysis, BrainMeta–ASD colocalization, and comparisons with previously reported genetic and molecular evidence for ASD.

After Benjamini–Hochberg correction across all valid gene–tissue tests for each method, JTI identified 264 significant associations representing 64 genes, and FUSION identified 561 associations representing 128 genes. The two methods shared 40 genes, with 152 genes in their union. Regional conditional/joint analysis retained the Brain Cerebellar Hemisphere *MAPT-AS1* model in both the core and extended primary analyses. BrainMeta–ASD colocalization evaluated 23 probes; none met the prespecified criteria for strong colocalization.

## Repository contents

- `scripts/` — Analysis and figure-preparation scripts for TWAS, SMR, colocalization, regional analyses, and external-evidence comparisons.
- `data/` — Selected processed results, plotting inputs, annotations, and provenance records.
- `docs/` — Analysis maps and reproducibility documentation.
- `reproducibility/` — The regional FUSION reproducibility archive and its supporting records.
- `third_party/` — Upstream FUSION and MR-JTI/PrediXcan-R code and associated license files.

The [figure and table mapping](Analysis_Figure_Table_Mapping.tsv) links the article's figures and tables to the corresponding repository materials. The [data guide](data/README.md) describes the processed results included in the repository.

## Data and reproducibility

Detailed tabulated results are provided in the article's supplementary materials. This repository provides selected processed summaries, plotting inputs, and analysis scripts.

Original GWAS summary statistics, expression-prediction weights, QTL datasets, and LD reference genotypes are obtained from their respective providers and are not redistributed in this repository. Access and use are subject to the conditions set by the original providers.

For module-specific input requirements and available execution and validation information, see the [analysis coverage and validation guide](docs/reproducibility/analysis_coverage_and_validation.md).

The [FUSION regional reproducibility archive](reproducibility/FUSION_GRCh38_Reproducibility.zip) provides commands, coordinate checks, and diagnostic records for the regional conditional/joint analyses.

## Citation and licensing

Please cite the article above and the original methods and data resources used in any work based on this repository.

Third-party components retain their respective upstream licenses; see the corresponding license files in `third_party/`. No license is specified for the repository as a whole.
