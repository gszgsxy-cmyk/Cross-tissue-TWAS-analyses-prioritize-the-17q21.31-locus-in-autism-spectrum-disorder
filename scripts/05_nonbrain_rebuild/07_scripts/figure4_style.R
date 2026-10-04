# Extracted from original Figure4A .Rhistory, theme and style only.
# Status label renamed to requested terminology; numerical results not sourced.
tm <- forest_theme(
base_family = "serif",
base_size = 9.5,
ci_pch = 15,
ci_col = "black",
ci_fill = "#3C5488",
ci_lty = 1,
ci_lwd = 1.3,
ci_Theight = 0.15,
refline_lwd = 1,
refline_lty = "dashed",
refline_col = "grey40",
footnote_cex = 0.72,
footnote_fontface = "italic",
footnote_col = "grey30",
core = list(
bg_params = list(fill = "#FFFFFF", col = NA),
fg_params = list(col = "black")
),
colhead = list(
fg_params = list(col = "black", fontface = 2)
)
)
# Color-blind friendly colors used only for the HEIDI status text.
col_nohet <- "#0072B2"
col_het   <- "#D55E00"
col_na    <- "#7F7F7F"
style_forest <- function(p, display_dat, status_vec, gene_col = NULL){
# Header bottom line
p <- add_border(
p,
part = "header",
where = "bottom",
gp = gpar(lwd = 1.2, col = "black")
)
# Center all columns except the first label column
if(ncol(display_dat) >= 2){
p <- edit_plot(
p,
col = 2:ncol(display_dat),
part = "body",
which = "text",
hjust = unit(0.5, "npc"),
x = unit(0.5, "npc")
)
}
# Gene symbols in italics for the eQTL panel
if(!is.null(gene_col)){
p <- edit_plot(
p,
col = gene_col,
part = "body",
which = "text",
gp = gpar(fontface = "italic")
)
}
# HEIDI-status column: explicit text + color
status_col <- ncol(display_dat)
idx_nohet <- which(status_vec == "No heterogeneity detected")
idx_het   <- which(status_vec == "Heterogeneity detected")
idx_na    <- which(status_vec == "HEIDI unavailable")
if(length(idx_nohet) > 0){
p <- edit_plot(
p, row = idx_nohet, col = status_col,
part = "body", which = "text",
gp = gpar(col = col_nohet, fontface = "bold")
)
}
if(length(idx_het) > 0){
p <- edit_plot(
p, row = idx_het, col = status_col,
part = "body", which = "text",
gp = gpar(col = col_het, fontface = "bold")
)
}
if(length(idx_na) > 0){
p <- edit_plot(
p, row = idx_na, col = status_col,
part = "body", which = "text",
gp = gpar(col = col_na, fontface = "bold")
)
}
p
}
