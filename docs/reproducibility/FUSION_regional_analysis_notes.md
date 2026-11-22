# Regional FUSION analysis: coordinates, selection and diagnostic evidence

This note connects the archived FUSION evidence to Supplementary Table S9C.
The evidence paths below refer to members of the unchanged
[FUSION reproducibility ZIP](../../reproducibility/FUSION_GRCh38_Reproducibility.zip).
They contain the recorded commands, analyses and diagnostic checks.

## Coordinates and input roles

Regional post-processing uses positions from the prepared LD-reference `.bim`
and prediction-model P0/P1 fields. GWAS summaries supply SNP identifiers,
association statistics and alleles for reference matching. The LD positions
and GTEx v8 model positions used here are on **GRCh38**, with GENCODE v26
gene annotation.

The core interval is chr17:45,414,716–46,836,265; the extended interval is
chr17:44,414,716–47,836,265. Regional mQTL analyses separately use the
**GRCh37/hg19** core, chr17:43,492,082–44,913,631. Those QTL coordinates do
not define the positions used in the FUSION LD reference.

The plotting variant uses the hg38 gene-track reference; the change is recorded
in `reproducibility/source/GRCh38_plot.diff`. Archived coordinate checks include:

| Evidence file inside the ZIP | Recorded result |
| --- | --- |
| `reproducibility/evidence/independent_GRCh38_coordinate_comparison.tsv` | All 32,169 chr17 LD variants have status MATCH_POSITION_ALLELES against the independently obtained GRCh38 reference. |
| `reproducibility/evidence/GWAS_reference_matching.tsv` | 13,780 same-orientation and 12,587 swapped-orientation variants: 26,367 matches; 5,802 variants missing from the GWAS input. |
| `reproducibility/evidence/model_coordinate_annotation_QC.tsv` | All 2,281 records satisfy P1 = P0 + 1 and match the strand-aware GENCODE v26 transcription start site. |
| `reproducibility/evidence/snp_allele_coordinate_QC.tsv` | Summary of the coordinate and allele checks. |

The reference identifier and compressed checksum are recorded in the archive's
`README_EN.md`. The archive identifies the reference used for these recorded checks. Unavailable GWAS statistics remain missing.

## Candidate sets and generated groups

Primary candidates overlap the predefined interval and meet marginal FUSION
all-tissue BH-FDR < 0.05. All-testable sensitivity sets use the available
overlapping models irrespective of marginal FDR.

| Analysis ID | Models | Genes | Predefined inclusion interval, GRCh38 | FUSION-generated group, GRCh38 |
| --- | ---: | ---: | --- | --- |
| core_353 | 353 | 28 | chr17:45,414,716–46,836,265 | chr17:44,834,540–47,428,979 |
| extended_358 | 358 | 33 | chr17:44,414,716–47,836,265 | chr17:44,236,113–47,606,009 |
| core_559 | 559 | 32 | chr17:45,414,716–46,836,265 | chr17:44,834,540–47,428,979 |
| extended_1011 | 1,011 | 82 | chr17:44,414,716–47,836,265 | chr17:43,906,472–48,430,072 |

Each analysis yields one merged group. FUSION derives boundaries from
LD-reference SNP positions, a 100-kb extension and merging of overlapping
intervals; models are assigned using P0/P1. Generated group boundaries
therefore differ from the interval used to select input models.

`reproducibility/evidence/regional_groups.tsv` records the independent boundary
reconstruction; historical logs provide group counts. `Archived_regional_checks.tsv`,
`four_analysis_QC.tsv` and the default-run `run.regional_boundaries.tsv` files
provide complementary stored checks.

## Global and locus-specific selection

The globally retained feature in each analysis is MAPT-AS1 / ENSG00000264589.2 /
GTExv8.ALL.Brain_Cerebellar_Hemisphere. Global and locus-specific selection
are distinct outputs.

| Default diagnostic run | Global retained | Locus-specific selected | Locus-specific negative variance evaluations | Computational SNP NaN Z | Missing-input SNP NA Z |
| --- | ---: | ---: | ---: | ---: | ---: |
| core_353 | 1 | 1 | 0 | 0 | 20 |
| extended_358 | 1 | 1 | 0 | 0 | 38 |
| core_559 | 1 | 10 | 260 | 502 | 20 |
| extended_1011 | 1 | 19 | 532 | 1,151 | 60 |

The all-testable locus-specific iterations include negative conditional
variance evaluations and extreme conditional Z values (maximum absolute
values about 8,253 and 119,159). Subsequent SNP conditioning produces negative
variances and computational NaN. The 10/19 selected models describe unstable
selections, not reliable independent signals. The 260/532 counts are
variance-evaluation records across iterations, not counts of distinct models.
The global guard excludes six/seven candidates in these diagnostic runs;
the final global set still contains one model.

The primary 353/358-model diagnostic selections show none of the corresponding
negative-variance or SNP-NaN instability. Computational NaN and missing-input
NA are separately counted and are not replaced by zero.

Detailed stage explanations are in the archive's
`Major2_numerical_stability_notes_EN.md`. Supporting evidence includes:

- `reproducibility/evidence/Major2_threshold_stability_summary.tsv`;
- `reproducibility/evidence/selection_stage_diagnostics.tsv`;
- `reproducibility/evidence/SNP_nonfinite_classification.tsv`;
- `reproducibility/evidence/regional_retained_models.tsv`.

## Default and common Z thresholds

| Input models | Recorded default absolute Z threshold |
| --- | ---: |
| 353 | 3.80528140728608 |
| 358 | 3.80876112046778 |
| 559 | 3.91757172039658 |
| 1,011 | 4.05818311062166 |

The default corresponds to two-sided P < 0.05/M for M input models. The
explicit common threshold was 4.058183. In the eight default/fixed diagnostic
comparisons, global retention, locus-specific selected sets and selection
trajectories were unchanged. The common threshold did not resolve the
all-testable numerical instability. The saved threshold summary records these comparisons.

## Historical records and diagnostic runs

The `diagnostic_20261016` directories contain diagnostic executions generated
on 16 October 2026 with R 3.6.0. They are distinct from original analysis-run
logs. `Original_vs_diagnostic_outputs.tsv` records comparisons with existing
historical outputs, including byte-identical default global outputs.

`reproducibility/Source_code_and_versions_EN.md` identifies the supplied local
and instrumented sources by checksum. A unique historical FUSION Git commit
and complete independent historical stdout/stderr were not recovered; the
current copies are identified as recovered local sources. Instrumentation
changes are recorded separately. R 4.4.1 in this archive applies to the
independent coordinate helper, not the diagnostic FUSION runs.

`reproducibility/Commands_and_parameters_EN.md` and
`Diagnostic_parameters_20261016.tsv` record options and input identifiers.
Re-execution requires the external GWAS, LD and prediction-weight resources.

## Current Supplementary Table S9C organization

**Supplementary Table S9C. Regional FUSION conditional/joint analysis and
gene–tissue model membership** combines:

- **Part A:** four regional summaries, generated boundaries, the threshold
  comparison and existing locus-specific diagnostic membership/notes.
- **Part B:** all 2,281 model-member records with the 32 original fields in
  source order.

The full source TSV is retained in the authors' local project records, not bundled
in this repository. Its SHA256 is
`a2d782898100e701850691964c1e5a7c93ad6ecfc2d91159d2ce3ad6d52577b9`.
The complete records are now in the article workbook's S9C Part B. The
repository's [S9C TSV](../../data/final_summaries/S9C_Sensitivity.tsv) preserves
the pre-existing Part A export. The ZIP's model-coordinate QC is a coordinate
subset, not the full 32-field member table.

The ZIP's original README, column dictionary and companion filenames capture
the earlier attachment organization. Their references to a separate member-data
attachment identify that archival snapshot; the current location is S9C Part B.
The unchanged archive is 810,699 bytes, with SHA256
`f52792140f945eb3d124c367afabc166de5194f58fe3d864ad79de509106293b`.
