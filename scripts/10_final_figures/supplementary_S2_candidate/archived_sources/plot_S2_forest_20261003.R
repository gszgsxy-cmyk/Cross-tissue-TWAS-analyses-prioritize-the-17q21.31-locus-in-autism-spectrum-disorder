if (.Platform$OS.type=='windows') Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')
library(data.table)
library(forestploter)
library(grid)
library(gridExtra)
root <- '__WORKSPACE_ROOT__'
out <- file.path(root,'outputs/Supplementary_Figures_FINAL_WORKING_20261003/S2_panels')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'outputs/2026-09-19_GTEx_nonbrain_SMR_rebuild/07_scripts/figure4_style.R'))
d <- fread(file.path(root,'outputs/fetalbrain_20260922/taierbrain_ASD_sensitivity/06_Hannon_FetalBrain_FDRsig_annotated.tsv'))
stopifnot(nrow(d)==5,all(d$FDR_FetalBrain<.05),!anyDuplicated(d$probeID),all(d$In_17q21.31_core),all(is.na(d$p_HEIDI)),all(is.na(d$nsnp_HEIDI)))
d <- d[order(hg19_BP)]
d[,`:=`(CI_lower=b_SMR-1.96*se_SMR,CI_upper=b_SMR+1.96*se_SMR,HEIDI_category='HEIDI unavailable')]
tab <- data.frame(CpG=d$probeID,`Genomic position`=paste0('chr',d$hg19_CHR,':',d$hg19_BP),`Top SNP`=d$topSNP,`Effect (95% CI)`=sprintf('%.4f (%.4f to %.4f)',d$b_SMR,d$CI_lower,d$CI_upper),` `=paste(rep(' ',22),collapse=' '),FDR=formatC(d$FDR_FetalBrain,format='f',digits=6),`HEIDI P value`='NA',`SNPs in HEIDI`='NA',`HEIDI assessment`=d$HEIDI_category,check.names=FALSE)
a <- max(abs(c(d$CI_lower,d$CI_upper)))
ticks <- pretty(c(-a,a),n=5);lim <- max(abs(ticks))
p <- forest(tab,est=d$b_SMR,lower=d$CI_lower,upper=d$CI_upper,ci_column=5,sizes=.55,ref_line=0,xlim=c(-lim,lim),ticks_at=ticks,xlab='SMR effect estimate (95% CI)',theme=tm)
p <- style_forest(p,tab,d$HEIDI_category)
pdf(NULL,width=22,height=20)
w <- convertWidth(sum(p$widths),'in',valueOnly=TRUE)+.5
h <- convertHeight(sum(p$heights),'in',valueOnly=TRUE)+1.15
dev.off()
figure <- arrangeGrob(textGrob('A   Hannon fetal-brain mQTL-SMR',x=unit(.01,'npc'),hjust=0,gp=gpar(fontfamily='serif',fontface='bold',fontsize=14)),p,ncol=1,heights=unit.c(unit(.4,'in'),unit(h-.4,'in')))
pdf(file.path(out,'S2A_Hannon_forest.pdf'),width=w,height=h,useDingbats=FALSE)
grid.newpage();grid.draw(figure);dev.off()
png(file.path(out,'S2A_Hannon_forest.png'),width=round(w*120),height=round(h*120),res=120)
grid.newpage();grid.draw(figure);dev.off()
fwrite(d[,.(probeID,hg19_CHR,hg19_BP,topSNP,b_SMR,se_SMR,CI_lower,CI_upper,FDR_FetalBrain,p_HEIDI,nsnp_HEIDI,HEIDI_category,Standardized_gene,In_17q21.31_core)],file.path(out,'S2A_plot_data.tsv'),sep='\t',na='NA')
cat('Forest plot saved: 5 CpGs; original forest theme retained.\n')
