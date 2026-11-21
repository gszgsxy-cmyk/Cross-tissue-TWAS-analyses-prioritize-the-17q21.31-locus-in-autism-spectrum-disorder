# Analysis coverage, inputs and validation

This guide distinguishes execution evidence, saved-result checks and external
input requirements. Complete tabulations are in the article supplementary
workbook; repository data provide selected results and plotting resources.

## Coverage by module

| Module | Available material | Validation and additional inputs |
| --- | --- | --- |
| [third_party](../../third_party/) | FUSION source, MR-JTI/PrediXcan-R files and upstream licenses. | Upstream implementations. The study's complete original JTI launch orchestration is not supplied; GWAS, weights and upstream dependencies are external. |
| [02_fusion](../../scripts/02_fusion/) | Regional input preparation, path adjustment and GRCh38 processing/plotting variants. | Source coverage and archived execution evidence. These scripts combine calculations and drawing; full GWAS, LD and weights are external. |
| [03_twas_summary_figures](../../scripts/03_twas_summary_figures/) | JTI/FUSION result processing and earlier heatmap/summary source. | Historical sequential records requiring association outputs and annotations. The FUSION record includes the transition to corrected GTEx v8 Whole Blood results. |
| [04_smr_execution](../../scripts/04_smr_execution/) | SMR shell execution, QTL extraction and non-brain completion commands. | Available execution stages, not one launcher for all resources. Requires SMR, GWAS, QTL BESD and LD inputs. |
| [05_nonbrain_rebuild](../../scripts/05_nonbrain_rebuild/) | Six R scripts and old-input mapping for 16-tissue aggregation, per-tissue FDR and forest plots. | Historical source chain requiring original SMR outputs and annotations. S11D supplies saved totals; the repository checker does not execute the R chain. |
| [06_coloc](../../scripts/06_coloc/) | Extraction, batch preparation, full-cis analysis and detailed result verification; saved 23-probe S10A summary. | R script enforces R 4.4.1/coloc 5.2.3. Full probe-level cis inputs/intermediates are external; saved posteriors and flags are checked. |
| [07_mqtl](../../scripts/07_mqtl/) | Hatton EUR processing/enrichment and Hannon fetal-brain post-processing. | Requires existing SMR outputs and annotations. Scripts can calculate summary FDR/enrichment when run; this update only inspects them. Complete results are in S11E–I. |
| [08_external_evidence](../../scripts/08_external_evidence/) | Benchmark/integration scripts, S13A/B flags, S14 trace/strict records and S4/S5 counts. | Historical processing needs original external resources. Current checks count saved flags and verify source fields without rerunning integration or enrichment. |
| [09_figures](../../scripts/09_figures/) | Hatton forest, Figure 4 assembly and overlap-panel source. | Uses historical TSVs and panels; complete final Figure 4 source-panel assembly is not bundled. |
| [Circular heatmaps](../../scripts/10_final_figures/circular_heatmaps/README.md) | Standalone R plotter, four grids and eight explicit order files. | Archived successful R execution and drawing-stream comparison on 11 November 2026, recorded in VALIDATION.md. Current input/grid/star checks; no repeated rendering. |
| [D panels](../../scripts/10_final_figures/cross_tissue_D/README.md) | Larger-axis script and complete 21×29 / 32×29 grids. | Original plot blocks retained with adapted input/output routing. Grid/flag checks pass; no current R syntax/execution test. |
| [Figure 5](../../scripts/10_final_figures/figure5_boundaries/README.md) | Annotation script and glyph-compatible vector base. | Check-only validates the base, coordinate calibration and in-memory MB-to-Mb unit correction. The corrected drawing was rendered with Poppler and matched the approved TIFF pixel for pixel. Drawing annotates the existing vector result rather than independently recreating it from FUSION inputs. Optional PDFium raster antialiasing can differ from historical Poppler. |
| [S4/S5](../../scripts/10_final_figures/supplementary_S4_S5/figures.py) | Existing drawing functions and saved three-source counts. | Check-only validates counts. Archived plot text was previously compared with the submitted supplementary PDF; no figures were regenerated here. |
| [S2 candidate](../../scripts/10_final_figures/supplementary_S2_candidate/README.md) | Five-CpG data, archived vector panels, assembly script and historical source. | Inputs/check-only tested. Adoption in the final article attachment awaits author confirmation. Rendering requires the specified Arial Bold font. |

Figures 2/3 C panels, their complete final composite assembly, and complete
portable final S1/S3 workflows are not supplied as standalone modules. Related
historical sources remain mapped in the
[figure/table guide](../../Analysis_Figure_Table_Mapping.tsv).

## Validation levels

| Level | Meaning |
| --- | --- |
| Executed checks | Integrity verifier and three Python check-only entry points read bundled inputs and saved results; no association calculation or figure rendering. |
| Syntax inspection | Python sources compiled without executing analysis code or writing bytecode. Shell/R sources inventoried and inspected; no current R or shell execution validation claimed. |
| Archived execution | Circular-heatmap R validation and FUSION diagnostic records document earlier runs, separately from current checks. |
| External-input workflows | Association, SMR, coloc and integration stages require original resources and/or historical intermediates. |
| Saved results/plotting inputs | Existing values and significance flags support result checks and visualization, without replacing raw analysis inputs. |

The current checking environment uses Python 3.12.14. Native Rscript is
unavailable and WSL access is restricted, so no new R syntax, execution or
rendering validation is reported. Documentation changes refresh the file
manifest but do not modify analysis code.

## Reproducible entry points

From the repository root:

```sh
python tools/verify_collection.py --report ../verification.json
python scripts/10_final_figures/figure5_boundaries/annotate_boundaries.py --check-only
python scripts/10_final_figures/supplementary_S4_S5/figures.py --check-only
python scripts/10_final_figures/supplementary_S2_candidate/assemble_candidate.py --check-only
```

The main checker uses the Python standard library. Figure checks use packages
in [requirements-plotting.txt](../../requirements-plotting.txt). The manifest
covers repository files excluding itself, Git metadata, Python bytecode and
generated `output` directories. The verifier writes only its requested report.

Checks cover saved TWAS totals, 23-probe coloc results, the 16-tissue summary,
regional model counts, Huang trace, three-source flags, and A/B versus D FDR
definitions. Coverage notes are reported separately. The verifier does not
currently check every field in S9C Part B: the full member table is supplied
in the article workbook. The archived model-coordinate QC is a subset of its
fields, as explained in the [regional notes](FUSION_regional_analysis_notes.md).

## Statistical definitions and external resources

- Primary TWAS uses all-29-tissue BH-FDR within each method; exploratory A/B
  panels use brain/non-brain group-specific FDR. D panels use all-tissue FDR.
  Plotting consumes stored flags and Z values without recalculating FDR.
- SMR correction sets follow the resource-specific study definitions;
  non-brain eQTL SMR uses per-tissue FDR. Unavailable HEIDI stays unavailable;
  HEIDI P > 0.05 means no detected heterogeneity.
- Regional primary and all-testable sensitivity sets use different inclusion
  rules. Global and locus-specific selections are distinguished.
- S13 counts Huang, Matoba and Gandal. SFARI/Fu and Walker have separate
  benchmark/context fields. S14's strict Huang set is a sensitivity comparison.

External inputs include ASD GWAS summaries, GTEx v8 weights, European LD
genotypes, BrainMeta and brain mQTL summaries, Hatton EUR and Hannon fetal-brain
resources, 16 non-brain GTEx QTL panels and the external gene-set studies.
Article S1 lists sources and versions; the archive records key FUSION input
identifiers. Provider access and licensing terms apply.

The [data guide](../../data/README.md) distinguishes article tables, repository
excerpts and plotting grids. [paths.example.json](../../config/paths.example.json)
and [configure_paths.py](../../tools/configure_paths.py) configure a separate
copy of historical scripts; the helper does not obtain external datasets.
