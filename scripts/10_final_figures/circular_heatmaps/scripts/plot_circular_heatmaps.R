# ASD TWAS circular heatmaps: plotting only, no association/FDR calculation.
# Relative paths; no saved R workspace or private source directory required.
options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)
}
for (pkg in c("ggplot2", "circlize")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg, ". Install it before running this script.")
  }
}
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(circlize))
argv <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
project_root <- dirname(dirname(script_file))
input_dir <- file.path(project_root, "data")
out <- file.path(project_root, "output")
overwrite <- FALSE
iarg <- 1L
while (iarg <= length(argv)) {
  opt <- argv[iarg]
  if (opt == "--overwrite") { overwrite <- TRUE; iarg <- iarg + 1L; next }
  if (!(opt %in% c("--data-dir", "--output-dir")) || iarg == length(argv)) {
    stop("Usage: Rscript scripts/plot_circular_heatmaps.R [--data-dir PATH] [--output-dir PATH] [--overwrite]")
  }
  if (opt == "--data-dir") input_dir <- argv[iarg + 1L]
  if (opt == "--output-dir") out <- argv[iarg + 1L]
  iarg <- iarg + 2L
}
panels <- data.frame(code = c("Figure2A", "Figure2B", "Figure3A", "Figure3B"),
  method = c("JTI", "JTI", "FUSION", "FUSION"),
  group = c("Brain", "Nonbrain", "Brain", "Nonbrain"),
  genes = c(40L,40L,76L,85L), tissues = c(13L,16L,13L,16L),
  expected_stars = c(152L,115L,263L,295L))
outputs <- file.path(out, c(paste0(panels$code, "_Circular_Heatmap.pdf"),
                           "Shared_TWAS_Z_Legend.pdf", "Plot_QC.tsv", "sessionInfo.txt"))
if (!overwrite && any(file.exists(outputs))) stop("Output files exist. Select another --output-dir or use --overwrite.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scale <- scale_fill_gradient2(low = "#3980b8", mid = "white", high = "#ef3b3c",
  midpoint = 0, limits = c(-6.1,6.1), na.value = "grey92")
qc <- list()
for (idx in seq_len(nrow(panels))) {
  cfg <- panels[idx,]; code <- cfg$code; method <- cfg$method; group <- cfg$group
  d <- read.delim(file.path(input_dir,paste0(code,"_",method,"_",group,"_FINAL_plotdata.tsv")),
                 check.names = FALSE, na.strings = "NA")
  go <- read.delim(file.path(input_dir,paste0(code,"_gene_order.tsv")),check.names=FALSE)
  to <- read.delim(file.path(input_dir,paste0(code,"_tissue_order.tsv")),check.names=FALSE)
  stopifnot(identical(as.integer(go$index),seq_len(nrow(go))),
            identical(as.integer(to$ring),seq_len(nrow(to))))
  genes <- go$gene; tissues <- to$tissue; labels <- to$label
  ng <- length(genes); nt <- length(tissues)
  key <- paste(d$DISPLAY_GENE,d$TISSUE,sep="|")
  stopifnot(ng == cfg$genes, nt == cfg$tissues, nrow(d) == ng*nt,
            !anyDuplicated(key), !anyDuplicated(genes), !anyDuplicated(tissues),
            setequal(genes,d$DISPLAY_GENE), setequal(tissues,d$TISSUE))
  fdrcol <- if(group == "Brain") "FDR_BRAIN" else "FDR_NONBRAIN"
  stopifnot(fdrcol %in% names(d), !anyNA(d$STAR),
    identical(d$STAR,ifelse(!is.na(d[[fdrcol]]) & d[[fdrcol]] < .05,"*","")),
    sum(d$STAR=="*") == cfg$expected_stars,
    all(abs(d$TWAS_Z[!is.na(d$TWAS_Z)]) <= 6.1))
  d$gene_clockwise_index <- match(d$DISPLAY_GENE,genes)
  d$tissue_outer_to_inner_index <- match(d$TISSUE,tissues)
  stopifnot(!anyNA(d$gene_clockwise_index),!anyNA(d$tissue_outer_to_inner_index))
  colors <- scale$map(d$TWAS_Z)
  stopifnot(identical(colors=="grey92",is.na(d$TWAS_Z)))
 # 270-degree arc, 90-degree opening in the upper-right quadrant.
 # Preserve the original gene order clockwise, from 3 o'clock to 12 o'clock.
 pdf(file.path(out,paste0(code,"_Circular_Heatmap.pdf")),width=180/25.4,height=180/25.4,
     useDingbats=FALSE,family="Helvetica")
 par(mar=c(0,0,0,0),oma=c(0,0,0,0))
 circos.clear()
 circos.par(start.degree=0,clock.wise=TRUE,gap.after=90,
   cell.padding=c(0,0,0,0),track.margin=c(0,0),
   canvas.xlim=c(-1.47,1.47),canvas.ylim=c(-1.47,1.47),
   points.overflow.warning=FALSE)
 circos.initialize(factors="heat",xlim=c(0,ng))
 circos.trackPlotRegion(ylim=c(0,nt),track.height=.64,bg.border=NA,panel.fun=function(x,y){
  i<-d$gene_clockwise_index;j<-nt-d$tissue_outer_to_inner_index
  circos.rect(i-1,j+.025,i,j+.975,col=colors,border="white",lwd=.14)
  star<-which(d$STAR=="*")
  circos.text(i[star]-.5,j[star]+.5,rep("*",length(star)),col=d$STAR_COLOR[star],
    cex=.38,font=2,facing="inside",niceFacing=TRUE)
  circos.text(seq_len(ng)-.5,rep(nt+.40,ng),genes,
    facing="clockwise",niceFacing=TRUE,adj=c(0,.5),cex=.59)
 })
 # The upper cut ends at 90 degrees. Direct labels align to their own rings.
 # Full source tissue labels are used; no numbering lookup is needed.
 radial<-sapply(seq_len(nt),function(j)circlize(ng,nt-j+.5)[1,2])
 text(.035,1.075,"Tissues (outer to inner)",adj=c(0,0),cex=.56,font=2)
 text(rep(.035,nt),radial,labels,adj=c(0,.5),cex=.53)
 text(-1.40,1.42,paste0(method," \u2014 ",if(group=="Brain")"Brain" else "Non-brain"),
   adj=c(0,.5),cex=.95,font=2)
 text(-1.40,1.32,"* Group-specific BH-FDR < 0.05",adj=c(0,.5),cex=.60,col="grey35")
 actual_arc<-get.cell.meta.data("cell.start.degree")-get.cell.meta.data("cell.end.degree")
 stopifnot(abs(actual_arc-270)<1e-9)
 circos.clear();dev.off()
  qc[[code]] <- data.frame(panel=code,genes=ng,tissues=nt,cells=nrow(d),
    missing_Z=sum(is.na(d$TWAS_Z)),stars=sum(d$STAR=="*"),
    color_min=-6.1,color_max=6.1,arc_degrees=270,gap_degrees=90)
}
write.table(do.call(rbind,qc),file.path(out,"Plot_QC.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
# Standalone common legend for later assembly. No per-panel color key.
pdf(file.path(out,"Shared_TWAS_Z_Legend.pdf"),width=45/25.4,height=100/25.4,
    family="Helvetica",useDingbats=FALSE)
par(mar=c(0,0,0,0))
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
text(.06,.95,"TWAS Z-score",adj=0,cex=.82,font=2)
z <- seq(-6.1,6.1,length.out=521)
for(k in seq_len(520))rect(.32,.19+(k-1)*.00125,.48,min(.84,.19+k*.00125+.0005),
   col=scale$map(z[k]),border=NA)
for(tick in c(-6,-3,0,3,6)) {
 yy <- .19+(tick+6.1)/12.2*.65
 segments(.48,yy,.515,yy,col="grey30",lwd=.5)
 text(.55,yy,as.character(tick),adj=0,cex=.75)
}
rect(.32,.085,.40,.12,col="grey92",border=NA)
text(.46,.1025,"Missing",adj=0,cex=.62)
dev.off()
cat("Saved four panels and a separate shared legend. No statistical analysis was run.\n")
