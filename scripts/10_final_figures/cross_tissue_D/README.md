# Figures 2D and 3D

Run `Rscript plot_D_larger_axis_only.R [output_directory]` with ggplot2, grid and scales.
Inputs are the byte-preserved final plotting TSVs in data/. The original larger-axis
plot blocks are unchanged; only input/output routing is adapted. Outputs are protected
against overwrite. Grids: 21 x 29 and 32 x 29. Stars use FDR_ALL < 0.05, unlike A/B.
Z values, missing values, factor order, fill range (-6.1, 6.1), and aesthetics are retained.
This preparation checks stored inputs and syntax where available; no new figure was drawn.
