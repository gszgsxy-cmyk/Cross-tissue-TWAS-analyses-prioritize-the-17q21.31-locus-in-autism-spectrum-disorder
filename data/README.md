# Data and supplementary tables

This directory contains selected table exports and plotting resources. The
article's supplementary workbook contains the complete tabulations.
[source_exports.json](source_exports.json) identifies the export snapshot,
source worksheets, header rows and record counts. Blank values and source
field definitions are retained.

| Article material | Repository resource | Contents |
| --- | --- | --- |
| S5/S6 | [Test counts](final_summaries/S5_Test_counts.tsv), [gene union](final_summaries/S6_Integrated.tsv) | Saved primary counts and the 152-gene union with method flags. |
| S9A/S9B | [FUSION archive](../reproducibility/FUSION_GRCh38_Reproducibility.zip) | Regional QC and diagnostic records; full feature tables are in the article workbook. |
| S9C Part A | [Regional summaries](final_summaries/S9C_Sensitivity.tsv) | 39 nonblank source records, including four analyses and explanatory/diagnostic rows. |
| S9C Part B | Article workbook; coordinate subset in the FUSION archive | Complete 2,281 model-member records and 32 fields are in the workbook. |
| S9D | Article workbook | Five-transcript cross-method summary; no separate export is bundled. |
| S10A/S10B | [Coloc posteriors](final_summaries/S10A_Coloc_primary.tsv), [workflow](../scripts/06_coloc/) | Saved results for 23 probes; full settings and QC are in the workbook. |
| S11A–I | [16-tissue summary](final_summaries/S11D_Tissue_summary.tsv), [SMR source](../scripts/04_smr_execution/), [mQTL source](../scripts/07_mqtl/) | Selected totals and source code; complete resource-specific results are in the workbook. |
| S12 | [External-evidence source](../scripts/08_external_evidence/) | Historical SFARI/Fu processing; external source inputs are required. |
| S13A/S13B | [Source flags](external_evidence/S13A_Molecular_convergence.tsv), [evidence details](external_evidence/S13B_Evidence_details.tsv) | Stored assignments for 152 genes. S13C enrichment results are in the workbook. |
| S14 Part A | [Huang source trace](external_evidence/Huang_S14_Source_Trace.tsv) | 22 original associations for the 10 overlapping genes. |
| S14 Part B | [Huang strict set](external_evidence/Huang_S14_Strict5.tsv) | Five strict combined-P associations and union-overlap flags. |
| Figures 2/3 A/B and D | [A/B grids and orders](../scripts/10_final_figures/circular_heatmaps/data/), [D grids](../scripts/10_final_figures/cross_tissue_D/data/) | Saved Z values, significance flags, missing cells and plotting orders. |

## Source definitions

S9C Part B contains core_353 (353), extended_358 (358), core_559 (559) and
extended_1011 (1,011) records, in source order. The original
`Supplementary_Data_1_FUSION_17q21.31_model_membership.tsv` is retained in the
authors' project records and is not bundled here.
The [regional notes](../docs/reproducibility/FUSION_regional_analysis_notes.md)
explain the archive's earlier companion-file labels.

For Huang et al. (2021), [original S1](https://doi.org/10.1371/journal.pgen.1009309)
contains 18 genes represented by 31 expression/splicing associations.
`Original_S1_row` refers to the original Excel row; `Discovery_P`,
`Replication_P` and `Combined_meta_P` come from columns R, U and X.
Only combined P from X defines the separate strict P < 4e-7 comparison.
The primary broad overlap remains 10 genes. S13 source counts use Huang,
Matoba and Gandal; SFARI/Fu are genetic benchmarks and Walker supplies
developmental context.

## Access

GWAS summaries, prediction weights, QTL datasets and LD genotypes are not
redistributed here; obtain them from the providers under their access terms
(article Table S1). For study materials covered by the article's
reasonable-request provision, contact the corresponding author. Selected
exports in this repository do not imply unrestricted access to all source data.
