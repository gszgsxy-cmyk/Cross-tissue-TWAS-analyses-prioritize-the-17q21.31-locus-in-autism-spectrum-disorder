# Figures 2A/B and 3A/B: circular heatmaps

Draw four 270-degree JTI/FUSION heatmaps and a separate shared legend from
the bundled data and gene/tissue order files.

From this directory, with R packages `ggplot2` and `circlize` installed:

```sh
Rscript scripts/plot_circular_heatmaps.R
Rscript scripts/plot_circular_heatmaps.R --data-dir data --output-dir output
```

Default inputs and outputs are located relative to the script, so the first
command can also be invoked using an absolute script path. Existing outputs
are protected; `--overwrite` permits replacement. Outputs are four panel
PDFs, `Shared_TWAS_Z_Legend.pdf`, `Plot_QC.tsv` and `sessionInfo.txt`.

| Panel | Method | Group | Genes | Tissues |
| --- | --- | --- | ---: | ---: |
| Figure2A | JTI | Brain | 40 | 13 |
| Figure2B | JTI | Non-brain | 40 | 16 |
| Figure3A | FUSION | Brain | 76 | 13 |
| Figure3B | FUSION | Non-brain | 85 | 16 |

Stars use the stored brain/non-brain group-specific BH-FDR < 0.05 flags.
`FDR_ALL` does not determine A/B stars. Genes follow the original rectangular
order clockwise from 3 o'clock; tissues run from outer to inner rings.
There is no clustering. Z limits are -6.1 to 6.1 with blue `#3980b8`, white
and red `#ef3b3c`; missing cells use `grey92`. The center and panels have no
individual legend. The shared legend ticks are -6, -3, 0, 3 and 6.

Inputs originate from `ASD.JTI.29tissues.updated.all_withFDR.tsv` and
`ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat`, including
the GTEx v8 Whole Blood panel. The bundled inputs suffice for plotting;
checksums are listed in [data_manifest.tsv](data_manifest.tsv).

Gu Z, et al. circlize implements and enhances circular visualization in R.
Bioinformatics. 2014;30:2811–2812. DOI: 10.1093/bioinformatics/btu393.
