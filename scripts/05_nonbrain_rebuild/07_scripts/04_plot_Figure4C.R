source("00_config.R")
library(forestploter)
library(grid)
library(gridExtra)
source("figure4_style.R")
d <- fread(file.path(root,"05_Figure4C/Figure4C_GTEx_nonbrain_SMR_forest_data.tsv"))
stopifnot(nrow(d)>0, all(d$FDR_recalc<.05), !anyDuplicated(d$Record_ID))
fmt <- function(x) ifelse(is.na(x),"NA",ifelse(x<.0001,formatC(x,format="e",digits=2),formatC(x,format="f",digits=4)))
tab <- data.frame(Gene=d$Gene,Tissue=gsub("_"," ",d$Tissue),Probe=d$Probe,
                  `Top SNP`=d$topSNP,`Effect (95% CI)`=sprintf("%.4f (%.4f to %.4f)",d$b_SMR,d$CI_lower,d$CI_upper),
                  ` `=paste(rep(" ",22),collapse=" "),FDR=fmt(d$FDR_recalc),
                  `HEIDI P value`=fmt(d$p_HEIDI),`SNPs in HEIDI`=ifelse(is.na(d$nsnp_HEIDI),"NA",as.character(d$nsnp_HEIDI)),
                  `HEIDI assessment`=d$HEIDI_assessment,check.names=FALSE)
a <- max(abs(c(d$CI_lower,d$CI_upper)))
ticks <- pretty(c(-a,a),n=5); lim <- max(abs(ticks))
p <- forest(data=tab,est=d$b_SMR,lower=d$CI_lower,upper=d$CI_upper,ci_column=6,
            sizes=.55,ref_line=0,xlim=c(-lim,lim),ticks_at=ticks,
            xlab="SMR effect estimate (95% CI)",theme=tm)
p <- style_forest(p,tab,d$HEIDI_assessment,gene_col=1)
# Use actual grob dimensions to avoid clipping the additional Tissue column.
pdf(NULL,width=22,height=20)
w <- convertWidth(sum(p$widths),"in",valueOnly=TRUE)+.5
h <- convertHeight(sum(p$heights),"in",valueOnly=TRUE)+1.1
dev.off()
figure <- arrangeGrob(textGrob("C",x=unit(.01,"npc"),hjust=0,
                               gp=gpar(fontfamily="serif",fontface="bold",fontsize=14)),
                      p,ncol=1,heights=unit.c(unit(.4,"in"),unit(h-.4,"in")))
stem <- file.path(root,"05_Figure4C/Figure4C_GTEx_nonbrain_SMR_forest")
pdf(paste0(stem,".pdf"),width=w,height=h,useDingbats=FALSE)
grid.newpage();grid.draw(figure);dev.off()
png(paste0(stem,".png"),width=round(w*300),height=round(h*300),res=300)
grid.newpage();grid.draw(figure);dev.off()
writeLines(sprintf("Rows: %s\nWidth: %.3f in\nHeight: %.3f in\nColors: #0072B2 / #D55E00 / #7F7F7F\n",nrow(d),w,h),file.path(root,"01_QC/figure_dimensions.txt"))
capture.output(sessionInfo(),file=file.path(root,"01_QC/plot_R_sessionInfo.txt"))
cat("Figure saved:",nrow(d),"rows,",w,"x",h,"inches\n")
