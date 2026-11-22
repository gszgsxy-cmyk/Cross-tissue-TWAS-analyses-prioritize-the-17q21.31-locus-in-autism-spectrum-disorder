# Figures 2D and 3D

Draw the cross-tissue heatmaps with x-axis labels at 12 pt and y-axis labels
at 13 pt. From this directory, with `ggplot2`, `grid` and `scales` available:

```sh
Rscript plot_D_larger_axis_only.R
Rscript plot_D_larger_axis_only.R output/custom
```

The default output directory is `output/`. Inputs in `data/` are
`Figure_2D_JTI_shared21_alltissueFDR_plotdata.tsv` and
`Figure_3D_FUSION_shared32_alltissueFDR_plotdata.tsv`. These contain 21 genes
across 29 tissues for JTI and 32 genes across 29 tissues for FUSION.

Outputs use stems `Figure_2D_JTI_LargerAxis_v1` and
`Figure_3D_FUSION_LargerAxis_v1`: vector PDF, 600-DPI TIFF and `_preview.png`.
`Axis_Font_Only_QC.tsv` records the input counts and plot comparisons. Existing
figure files are not overwritten.

Stars use stored `FDR_ALL < 0.05` across all 29 tissues, distinct from A/B
group-specific FDR. Existing Z values, ordering, -6.1 to 6.1 color limits
and `grey85` missing cells are retained.
