# Final article and historical file identifiers

The final article uses S12 for SFARI/Fu genetic benchmarks, S13 for external molecular
evidence, and S14 for Huang source trace and strict-threshold comparison. Historical
inputs named S13A/B_SFARI or S14A/B_Fu2022 retain their original filenames; these are
not the final numbered tables. Their original names are retained because the
historical scripts reference them.

The non-brain scripts named Figure4C in scripts/05_nonbrain_rebuild correspond to
the non-brain forest plot now reported as Supplementary Figure S3. Final Figure 4C
is Hatton EUR whole-blood mQTL-SMR, with its preparation in scripts/09_figures.
Fourteen source gene symbols in the non-brain records represent 13 stable IDs;
KIZ and PLK1S1 share an ID. Preserve all 59 associations.

The long FUSION processing record in scripts/03_twas_summary_figures starts by
inspecting the historical all29.allchr file, removes old Whole Blood entries,
and subsequently loads the corrected GTExv8WB.updated.GENCODEv26 file. It is a
sequential historical record, not a one-step current entry point. The corrected
Whole Blood panel has 8,648 rows, all named GTExv8.ALL.Whole_Blood; the old
GTEx.Whole_Blood panel is absent from the locked final association file.

Final input basenames and SHA256 identifiers are in data/final_input_identifiers.tsv.
Full association files, prediction weights, GWAS, LD and QTL inputs are not bundled.
The available TSVs are saved-result excerpts consumed by the read-only verifier.
Historical analysis scripts retain their original input requirements.

Figures 2/3 A-B use brain/non-brain-specific FDR; D uses all-29-tissue FDR.
The shared plotted Z range is -6.1 to 6.1. Inputs and factor/order files are preserved.
Current module and figure coverage is in Analysis_Figure_Table_Mapping.tsv.
The author confirmed on 2026-11-22 that the revised S2 layout was used in the final Supplementary Figures submitted to the journal. The historical directory and script names are retained.

The original FUSION reproducibility ZIP is preserved byte for byte, including
historical Major1/2 filenames and companion-workbook identifiers. Those identifiers
describe the preparation snapshot. The final article supplies the current workbook
with the 2,281-row, 32-field member records integrated into **S9C Part B**.
S9C Part A retains the summaries, generated boundaries and diagnostic comparisons.
The full source TSV is retained in the authors' local project records for
provenance and is not bundled in this repository; the ZIP itself is unchanged.
The [regional notes](reproducibility/FUSION_regional_analysis_notes.md) explain
how its historical companion labels map to the current article organization.

Module inputs and current versus archived validation evidence are described in
the [coverage guide](reproducibility/analysis_coverage_and_validation.md).
