# ASD TWAS circular heatmaps

Reproducible plotting code for the JTI/FUSION brain and non-brain panels.
This package draws the approved 270-degree heatmaps and a separate shared
legend. It performs no TWAS, FUSION model selection, or FDR recalculation.

## Run

Install R packages once:

```r
install.packages(c("ggplot2", "circlize"))
```

From the package directory:

```sh
Rscript scripts/plot_circular_heatmaps.R
```

The script finds `data/` relative to its own location, so it can also be
called from another working directory. Optional paths:

```sh
Rscript scripts/plot_circular_heatmaps.R --data-dir data --output-dir output
```

Existing outputs are protected by default. Use `--overwrite` explicitly to
replace them. Outputs are four vector panel PDFs, a shared legend PDF,
`Plot_QC.tsv`, and runtime `sessionInfo.txt`. Generated output is git-ignored.
No RData workspace, Python, external plotting scripts, or local model files
are required. Tested package versions are recorded in VALIDATION.md.

## Locked scientific definitions

| Panel | Method | Tissue group | Genes | Tissues | Cells | Stars | Missing Z |
|---|---|---|---:|---:|---:|---:|---:|
| Figure2A | JTI | Brain | 40 | 13 | 520 | 152 | 238 |
| Figure2B | JTI | Non-brain | 40 | 16 | 640 | 115 | 312 |
| Figure3A | FUSION | Brain | 76 | 13 | 988 | 263 | 468 |
| Figure3B | FUSION | Non-brain | 85 | 16 | 1360 | 295 | 755 |

Use the provided Z values and existing STAR/STAR_COLOR fields. Stars denote
brain-group or non-brain-group BH-FDR < 0.05, not all-29-tissue FDR.
FDR_ALL is retained for provenance but does not determine these stars.
No new genes are selected and missing values remain NA.

Genes follow the original rectangular top-to-bottom order, clockwise from
3 o'clock to 12 o'clock; tissues follow the original left-to-right order,
outside to inside. Explicit order files replace the private R workspace.
No clustering is performed.

Shared color limits: -6.1 to 6.1, blue #3980b8, white midpoint, red #ef3b3c,
original ggplot2 Lab interpolation, missing grey92. The shared legend ticks
are -6, -3, 0, 3, 6. There is no text in the center or individual color key.
The shared missing swatch represents these A/B panels; other panels should
use the same missing-color setting if that swatch is shared in an assembly.

## Data provenance

The four plot-data TSVs were checked row by row against the locked association
results and submitted supplementary tables. The underlying result basenames
are `ASD.JTI.29tissues.updated.all_withFDR.tsv` and
`ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat`.
The FUSION Whole Blood results use GTEx v8 Whole Blood weights.
The TSVs are unchanged copies of the verified plotting inputs.

The package includes only the 11-column panel data and explicit gene/tissue
orders. It contains no personal paths, model-weight files, saved workspaces,
private logs, manuscript files, credentials, or repository configuration.
`data_manifest.tsv` records relative filenames and SHA256 checksums.

## circlize citation

Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and enhances
circular visualization in R. Bioinformatics. 2014;30(19):2811-2812.
DOI: 10.1093/bioinformatics/btu393.
