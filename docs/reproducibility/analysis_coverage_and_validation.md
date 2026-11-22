# Analysis inputs, entry points and coverage

Use the table below to locate a module and its inputs. Complete result tables
are in the article workbook; selected exports are indexed in the
[data guide](../../data/README.md).

| Module | Available code or entry point | Input requirements and coverage |
| --- | --- | --- |
| [JTI/FUSION implementations](../../third_party/) | Selected upstream code and licenses | GWAS, weights and dependencies are external. The complete original JTI launch orchestration is not supplied. |
| [Regional FUSION](../../scripts/02_fusion/) | `prepare_allmodel_inputs.R`, path repair and GRCh38 post-processing variants | Requires association results, LD and weight files. Scripts combine computation and plotting. Commands and saved diagnostics are in the FUSION ZIP. |
| [TWAS processing](../../scripts/03_twas_summary_figures/) | Historical JTI/FUSION processing and heatmap records | Requires association outputs and annotations; the sequential FUSION record includes the corrected GTEx v8 Whole Blood transition. |
| [SMR execution](../../scripts/04_smr_execution/) | Shell commands and non-brain completion scripts | Requires SMR, GWAS, QTL BESD and LD inputs; these are resource-specific stages, not one launcher for every resource. |
| [Non-brain SMR](../../scripts/05_nonbrain_rebuild/) | `07_scripts/00_config.R` through the aggregation and plotting scripts | Requires original SMR outputs and annotations. Saved 16-tissue totals are in S11D. |
| [BrainMeta coloc](../../scripts/06_coloc/) | `extract_all.sh`, `prepare_batch.py`, `run_final_23probe_fullcis_coloc.R`, `verify_final.R` | Full cis inputs/intermediates are external. The analysis script requires R 4.4.1/coloc 5.2.3; the saved summary has 23 probes. |
| [mQTL follow-up](../../scripts/07_mqtl/) | Hatton and Hannon processing/enrichment scripts | Requires existing SMR outputs and annotations; complete tables are S11E–I. |
| [External evidence](../../scripts/08_external_evidence/) | `prepare.py`, `analyze.py` | Original source files are external. Bundled S13 flags and S14 records support saved-result checks. |
| [Figure 4 source](../../scripts/09_figures/) | Hatton forest, overlap and panel assembly records | Historical TSVs and source panels are not all bundled; complete end-to-end reconstruction is not provided. |
| [Figures 2/3 A–B](../../scripts/10_final_figures/circular_heatmaps/README.md) | R plotter with four grids and explicit order files | Bundled plotting inputs; prior R execution is recorded in VALIDATION.md. |
| [Figures 2D/3D](../../scripts/10_final_figures/cross_tissue_D/README.md) | Larger-axis R plotter and 21×29 / 32×29 grids | Bundled plotting inputs. |
| [Figure 5](../../scripts/10_final_figures/figure5_boundaries/README.md) | Python annotation of the supplied vector base | Coordinate/unit checks and optional PDFium rasterization; this does not reconstruct FUSION statistics. |
| [S4/S5](../../scripts/10_final_figures/supplementary_S4_S5/figures.py) | Plotting functions and saved three-source counts | Bundled inputs; `--check-only` checks existing counts. |
| [S2 layout](#supplementary-figure-s2) | Existing vector panels and Python assembly | Author confirmed use in the final supplementary submission on 2026-11-22. Requires an Arial Bold TTF file for rendering. |

The complete final Figures 2/3 composite assembly and portable final S1/S3
workflows are not supplied. Related historical source remains in the
[figure/table map](../../Analysis_Figure_Table_Mapping.tsv).
Figures 2C/3C were manually assembled by the authors; no standalone plotting
code is provided for these panels.

## Supplementary Figure S2

The assembly requires `pypdf`, `reportlab` and an Arial Bold TTF file for
the headings. Inputs in
`scripts/10_final_figures/supplementary_S2_candidate/data/` are
`S2A_Hannon_forest.pdf`, `S2B_Hannon_Hatton_Venn.pdf` and `S2A_plot_data.tsv`.

From the repository root:

```sh
python scripts/10_final_figures/supplementary_S2_candidate/assemble_candidate.py --check-only
python scripts/10_final_figures/supplementary_S2_candidate/assemble_candidate.py --bold-font PATH_TO_ARIAL_BOLD_TTF --output-dir output/S2
```

The second command writes `output/S2/Supplementary_Figure_S2_Layout_v1.pdf`.
An existing PDF is not overwritten. The entry point produces PDF only.

## Bundled-input checks

Save generated plot files under an `output/` directory or outside the repository.

From the repository root:

```sh
python tools/verify_collection.py --report ../verification.json
python scripts/10_final_figures/figure5_boundaries/annotate_boundaries.py --check-only
python scripts/10_final_figures/supplementary_S4_S5/figures.py --check-only
```

The main verifier uses the standard library. Figure checks use packages in
[requirements-plotting.txt](../../requirements-plotting.txt). These commands
read bundled inputs and saved results; they do not run association analyses
or draw figures. Syntax checks and earlier execution records are separate
from a new end-to-end analysis run.

## Paths and statistical definitions

Copy [paths.example.json](../../config/paths.example.json) to
`config/paths.local.json`, set the roots for your files, then use the existing
helper to make a separate runnable copy:

```sh
python tools/configure_paths.py --config config/paths.local.json --output ../configured-run
```

Placeholders retain historical subdirectory names. Check those expected inputs
against your own resources before running a module. The configured copy has
different code bytes; the original collection's manifest verifies the published
source snapshot. The helper does not obtain external data.

- Primary TWAS uses all-29-tissue BH-FDR within each method. Exploratory A/B
  panels use brain/non-brain group-specific FDR; D stars use all-tissue FDR.
- Non-brain eQTL SMR uses within-tissue FDR. Other SMR resources follow their
  resource-specific definitions; unavailable HEIDI remains unavailable.
- Regional primary and all-testable sets use different inclusion rules.
  [Regional notes](FUSION_regional_analysis_notes.md) distinguish global/locus
  selection, genome builds and numerical stability.
- S13 counts Huang, Matoba and Gandal; SFARI/Fu and Walker serve separate
  benchmark/context roles. The strict Huang comparison is separate in S14.

Sources and versions are listed in article S1. Provider access terms apply;
study materials covered by the article's reasonable-request provision should
be requested from the corresponding author.
