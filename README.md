# ASD TWAS Analysis Code

Analysis code, processed result summaries and figure resources for:

**Cross-tissue transcriptome-wide association analyses prioritize the 17q21.31 locus in autism spectrum disorder.**
Xiaoyang Li, Yuan Zhang, Ruofei Li, et al.
*Progress in Neuropsychopharmacology & Biological Psychiatry*, article 112014.
[DOI: 10.1016/j.pnpbp.2026.112014](https://doi.org/10.1016/j.pnpbp.2026.112014).

The study integrates ASD GWAS summary statistics with genetically predicted
expression across 29 GTEx v8 tissues, followed by molecular-QTL analyses,
regional conditional/joint analysis, colocalization and external-evidence
comparisons. This repository provides analysis scripts, selected tabulated
results, plotting inputs and a regional FUSION reproducibility archive.

## Study overview and main findings

JTI and FUSION were applied to 13 brain and 16 non-brain tissues.

| Method | Valid gene–tissue tests | Significant associations | Significant genes |
| --- | ---: | ---: | ---: |
| JTI | 296,168 | 264 | 64 |
| FUSION | 167,057 | 561 | 128 |

The primary analysis uses Benjamini–Hochberg FDR < 0.05 across all valid tests
in the 29 tissues, separately for each method. The methods share 40 significant
genes; their union contains 152 genes. FUSION uses **GTEx v8 Whole Blood
prediction weights** and GENCODE v26 annotations.

Brain and non-brain analyses use separate within-group BH corrections as an
exploratory analysis. JTI identifies 40 genes in each group, with 21 shared;
FUSION identifies 76 and 85, with 32 shared. Figures 2A/B and 3A/B use these
group-specific significance flags; Figures 2D and 3D use all-tissue FDR.

Follow-up analyses include brain eQTL-SMR, brain mQTL-SMR, Hatton EUR
whole-blood mQTL-SMR, Hannon fetal-brain mQTL-SMR and eQTL-SMR across 16
non-brain tissues. They identify 22 brain eQTL probes, 21 brain mQTL probes
and 59 non-brain probe–tissue associations. BrainMeta–ASD colocalization
examines 23 probes using PP.H4 ≥ 0.80 and PP.H4 > 3 × PP.H3 as the criterion
for strong colocalization; none meets both conditions.

Regional FUSION analyses use the GRCh38 core interval
chr17:45,414,716–46,836,265 and an interval extended by 1 Mb on each side.

| Regional analysis | Models | Genes |
| --- | ---: | ---: |
| Core primary: all-tissue FDR < 0.05 | 353 | 28 |
| Extended primary: all-tissue FDR < 0.05 | 358 | 33 |
| Core all-testable sensitivity | 559 | 32 |
| Extended all-testable sensitivity | 1,011 | 82 |

The MAPT-AS1 Cerebellar Hemisphere model is retained in the primary and
extended analyses. All-testable sensitivity analyses provide global-retention
evidence; their locus-specific numerical diagnostics are documented in the
[regional analysis notes](docs/reproducibility/FUSION_regional_analysis_notes.md).

Among the 152 union genes, 82 have support from at least one of the Huang,
Matoba or Gandal studies, and 15 from at least two; MAPT and XRN2 have support
from all three. These are complementary evidence sources rather than
independent replication. SFARI/Fu comparisons provide genetic benchmarks,
and Walker findings provide developmental-brain context.

## Analysis workflow and repository structure

| Analysis or resource | Repository location | Inputs and available outputs |
| --- | --- | --- |
| JTI/FUSION association methods and TWAS summaries | [MR-JTI](third_party/MR-JTI/UPSTREAM_README.md), [FUSION](third_party/FUSION/README.md), [TWAS processing](scripts/03_twas_summary_figures/) | Upstream implementations and study processing records; external GWAS/weights and final association files feed the summaries. Bundled [saved summaries](data/final_summaries/) contain test counts and the 152-gene set. |
| Regional FUSION conditional/joint analysis | [scripts/02_fusion](scripts/02_fusion/) | Input preparation and GRCh38 post-processing/plotting variants; commands, coordinate QC and diagnostic records in the [archive](reproducibility/FUSION_GRCh38_Reproducibility.zip). |
| SMR execution for brain and non-brain QTL resources | [scripts/04_smr_execution](scripts/04_smr_execution/) | SMR shell commands, including non-brain completion; requires SMR, GWAS, QTL BESD files and LD genotypes. |
| Non-brain SMR aggregation and forest plots | [scripts/05_nonbrain_rebuild](scripts/05_nonbrain_rebuild/) | Existing SMR outputs and annotations → 16-tissue summaries and forest-plot data; [saved tissue summary](data/final_summaries/S11D_Tissue_summary.tsv) included. |
| BrainMeta–ASD colocalization | [scripts/06_coloc](scripts/06_coloc/) | Extraction, preparation, 23-probe full-cis analysis and result checks; [posterior summaries](data/final_summaries/S10A_Coloc_primary.tsv) included. Full cis inputs are external. |
| Hatton EUR and Hannon fetal-brain mQTL follow-up | [scripts/07_mqtl](scripts/07_mqtl/) | Post-processing of existing SMR outputs, annotation, enrichment and cross-resource comparisons; complete tabulations in S11E–S11I. |
| External genetic and molecular evidence | [scripts/08_external_evidence](scripts/08_external_evidence/), [data/external_evidence](data/external_evidence/) | SFARI/Fu and molecular-evidence integration code; bundled S13 flags, Huang S14 source records and saved figure counts. Source databases are obtained separately. |
| Main Figure 4 preparation | [scripts/09_figures](scripts/09_figures/) | Hatton forest and panel assembly records using existing results and source panels. |
| Figures 2/3 A–B and D panels | [circular heatmaps](scripts/10_final_figures/circular_heatmaps/README.md), [D panels](scripts/10_final_figures/cross_tissue_D/README.md) | Bundled plotting grids and explicit orders → circular A/B panels, a shared legend and larger-axis D panels. |
| Figure 5 and supplementary figure resources | [Figure 5](scripts/10_final_figures/figure5_boundaries/README.md), [S4/S5 script](scripts/10_final_figures/supplementary_S4_S5/figures.py), [S2 layout candidate](scripts/10_final_figures/supplementary_S2_candidate/README.md) | Vector-base boundary annotation and MB-to-Mb axis-unit correction, saved-count plots and an archived-panel S2 assembly module, with module-specific inputs. |
| Data, provenance and detailed mappings | [data](data/README.md), [docs](docs/final_analysis_map.md), [figure/table mapping](Analysis_Figure_Table_Mapping.tsv) | Export provenance, input identifiers and links between article numbering and historical filenames. |

Module coverage and external-input requirements are detailed in
[analysis coverage and validation](docs/reproducibility/analysis_coverage_and_validation.md).

## Data availability and supplementary materials

The article's supplementary workbook contains the complete tabulated results.
Repository [data excerpts](data/README.md) and plotting grids provide selected
processed results for checking and figure preparation. Original GWAS, prediction
weights, LD genotypes and QTL resources are obtained from their providers under
the applicable access conditions.

**Supplementary Table S9C. Regional FUSION conditional/joint analysis and
gene–tissue model membership** has two parts in one worksheet:

- **Part A:** four regional-analysis groups, sensitivity summaries,
  FUSION-generated boundaries and existing diagnostic comparisons.
- **Part B:** all **2,281 model-membership records and 32 fields**, retaining
  the original source order. The four groups contain 353, 358, 559 and 1,011 records.

The complete model-membership records are in **S9C Part B** of the article
workbook. The original source TSV is retained in the authors' local project
records for provenance; the full TSV is not bundled in this repository. The repository's
`data/final_summaries/S9C_Sensitivity.tsv` is an export of the existing Part A
content, including summary and explanatory/diagnostic rows; it is not the full
32-field Part B table.

The [data guide](data/README.md) maps S9A–S9D, S10, S11, S13 and S14 to the
available files. For Huang et al. (2021), the primary comparison uses the 18
genes represented by 31 expression/splicing associations in original S1.
S14 Part A traces the 22 associations for the 10 overlapping genes. S14 Part B
reports the separate combined-meta-analysis P < 4.0 × 10⁻⁷ comparison:
five Huang genes, of which MAPT and KIZ overlap the union. Its three-source
sensitivity counts are 80/10, with MAPT supported by all three sources.

## Software and usage

| Component | Requirements and recorded versions |
| --- | --- |
| Repository integrity checker | Python 3.9+; standard library only |
| Python figure modules | Packages in [requirements-plotting.txt](requirements-plotting.txt); module READMEs identify required and optional dependencies |
| Circular heatmaps | R 4.4.1, ggplot2 3.5.1 and circlize 0.4.16 in the archived plotting validation |
| Cross-tissue D panels | R with ggplot2, grid and scales |
| BrainMeta colocalization | R 4.4.1 and coloc 5.2.3, enforced by the analysis script; data.table for input handling |
| SMR analyses | SMR 1.3.1, with the relevant QTL files and LD reference |
| Regional FUSION diagnostics | R 3.6.0 and dependency records in the reproducibility archive |

From the **repository root**, the following commands check bundled inputs and
saved results without running association analyses or drawing figures:

```sh
python tools/verify_collection.py --report ../verification.json
python scripts/10_final_figures/figure5_boundaries/annotate_boundaries.py --check-only
python scripts/10_final_figures/supplementary_S4_S5/figures.py --check-only
python scripts/10_final_figures/supplementary_S2_candidate/assemble_candidate.py --check-only
```

The main checker validates the file manifest, archive integrity, saved result
counts and plotting grids. Its report distinguishes failures from coverage
notes. Rendering commands and package installation instructions are in the
linked plotting-module READMEs.

Historical analysis scripts use configurable input roots. The
[path configuration example](config/paths.example.json) and
[copy configuration helper](tools/configure_paths.py) support a separate
path-configured copy. Prepare the module's external inputs before execution.

## Validation and reproducibility resources

| Resource | Validation evidence |
| --- | --- |
| File collection, saved summaries and source trace | Integrity checker executed against the bundled files; counts, identifiers and stored flags checked |
| Python figure modules | Check-only entry points executed; Python scripts syntax-checked |
| Circular heatmaps | Archived successful R execution and drawing-stream comparison; current input integrity checks |
| R analysis scripts and D-panel plotting | Dependencies and source coverage listed in the [validation guide](docs/reproducibility/analysis_coverage_and_validation.md) |

The [FUSION reproducibility archive](reproducibility/FUSION_GRCh38_Reproducibility.zip)
contains **58 files, 810,699 bytes**, with an internal manifest covering the
other 57 files. Start with `README_EN.md` inside the ZIP for its commands,
code versions, coordinate checks and diagnostic records.

SHA256:
`f52792140f945eb3d124c367afabc166de5194f58fe3d864ad79de509106293b`.

The [regional technical notes](docs/reproducibility/FUSION_regional_analysis_notes.md)
explain genome builds, global/locus-specific selection, archived diagnostic
runs and the current S9C Part A/B organization. The
[coverage guide](docs/reproducibility/analysis_coverage_and_validation.md)
identifies the inputs and execution evidence available for each module.

## Citation and licenses

Please cite Li et al. (2026), **Cross-tissue transcriptome-wide association analyses
prioritize the 17q21.31 locus in autism spectrum disorder**, *Progress in
Neuropsychopharmacology & Biological Psychiatry*, article 112014,
DOI: 10.1016/j.pnpbp.2026.112014, and the methods/resources used in your analysis.

Upstream FUSION and its derived scripts retain the
[GPLv3 license](third_party/FUSION/LICENSE). MR-JTI/PrediXcan-R retains the
[MIT license](third_party/MR-JTI/LICENSE). Upstream descriptions are preserved
under `third_party`. No collection-wide license has been assigned.
