# Data resources and supplementary-table mapping

These TSVs are direct value exports from the study workbook at the export snapshot (SHA256
`a92841ea4e8c7b91ca1e0b0db2b4681f291e4d2f16d9a26c2308de7267889bac`).
[`source_exports.json`](source_exports.json) records worksheets, original header rows and record counts.
No association testing, BH adjustment, gene matching or external-data analysis was run.
Blank source values remain blank; booleans and full available floating-point precision
are preserved. These excerpts do not replace the article's complete workbook.
The subsequent S9C Part A/B organization retains the exported content; the
export-time checksum remains provenance for these unchanged TSVs.

## Article tables and repository resources

| Article material | Repository resource | Coverage |
| --- | --- | --- |
| S5/S6 | [Test counts](final_summaries/S5_Test_counts.tsv), [integrated genes](final_summaries/S6_Integrated.tsv) | Saved tests/associations and 152-gene union, with method flags. |
| S9A/S9B | [FUSION archive](../reproducibility/FUSION_GRCh38_Reproducibility.zip) | Regional QC, retained-output comparisons and diagnostic records. Complete S9A feature details and S9B summaries are in the article workbook. |
| S9C Part A | [S9C_Sensitivity.tsv](final_summaries/S9C_Sensitivity.tsv) | Pre-existing value export: 39 nonblank records, including four main analysis rows, diagnostic/annotation rows and explanatory material. It is not a new clean 39-analysis dataset. |
| S9C Part B | Article supplementary workbook; archived coordinate subset in the FUSION ZIP | Complete 2,281 records × 32 fields are in the article workbook. The archive's model-coordinate QC contains selected coordinate fields, not the complete member table. |
| S9D | Article supplementary workbook | Descriptive five-transcript cross-method summary. No separate S9D data export is bundled. |
| S10A/S10B | [S10A posterior results](final_summaries/S10A_Coloc_primary.tsv), [coloc workflow](../scripts/06_coloc/) | Saved 23-probe posteriors; complete settings/QC are in article S10B and require full probe-level inputs for execution. |
| S11A–I | [S11D tissue summary](final_summaries/S11D_Tissue_summary.tsv), [mQTL source](../scripts/07_mqtl/), [SMR source](../scripts/04_smr_execution/) | Saved 16-tissue summary plus execution/post-processing code; complete SMR, HEIDI, overlap and attrition tables are in the article workbook. |
| S12 | [Benchmark integration](../scripts/08_external_evidence/) | Historical scripts and their input filenames map to the final genetic-benchmark tables; external source inputs are obtained separately. |
| S13A/S13B | [Convergence flags](external_evidence/S13A_Molecular_convergence.tsv), [evidence details](external_evidence/S13B_Evidence_details.tsv) | All 152 genes and stored source assignments. S13C enrichment results are in the article workbook. |
| S14 Part A | [Huang source trace](external_evidence/Huang_S14_Source_Trace.tsv) | 22 original associations for 10 overlapping genes. |
| S14 Part B | [Huang strict five](external_evidence/Huang_S14_Strict5.tsv) | Five strict combined-P associations and union-overlap flags. |
| Figures 2/3 | [Circular grids/orders](../scripts/10_final_figures/circular_heatmaps/data/), [D grids](../scripts/10_final_figures/cross_tissue_D/data/) | Existing Z values, missing cells, significance flags and orders for plotting. |

## S9C model-member provenance

The current title is **Supplementary Table S9C. Regional FUSION conditional/joint
analysis and gene–tissue model membership**. Part A contains sensitivity
summaries, generated boundaries and existing diagnostic comparisons. Part B
contains all 2,281 records in original order: core_353 (353), extended_358 (358),
core_559 (559) and extended_1011 (1,011).

The source TSV, `Supplementary_Data_1_FUSION_17q21.31_model_membership.tsv`, is
retained in the authors' local project records for provenance and is not bundled
in this repository. Its SHA256 is
`a2d782898100e701850691964c1e5a7c93ad6ecfc2d91159d2ce3ad6d52577b9`.
The current article supplies the full member records in **S9C Part B**. The
older companion-file labels inside the unchanged archive describe its
preparation snapshot; see the [regional notes](../docs/reproducibility/FUSION_regional_analysis_notes.md).

## External molecular evidence

Huang sources: Huang et al. (2021), original S1 Table,
https://doi.org/10.1371/journal.pgen.1009309 . Part A of final S14 contains
22 associations for the 10 overlapping genes. Part B contains five strict associations.
Discovery_P, Replication_P and Combined_meta_P correspond to original Excel R, U and X.
Only X determines strict combined P < 4e-7. Both expression and splicing are included.
Original_S1_row preserves the original Excel row, not the row in these exported TSVs.
Missing splicing IDs/statistics are not inferred. Stable IDs use the project's mapping.
The primary Huang broad set (18 genes/31 associations) and 10-gene overlap are retained;
the strict five-gene comparison is separate. SFARI/Fu are genetic benchmarks and Walker
is developmental context, outside the final Huang/Matoba/Gandal three-source count.
