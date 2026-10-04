source("00_config.R")
x <- rbindlist(lapply(tissues,function(t) fread(file.path(root,"02_per_tissue_merged",paste0(t,"_allchr_SMR_withFDR.tsv")))) ,fill=TRUE)
stopifnot(uniqueN(x$Tissue)==16, !anyDuplicated(x$Record_ID))
sig <- x[!is.na(FDR_recalc) & FDR_recalc < .05]
setorder(sig, ProbeChr, Probe_bp, probeID, Tissue)
summary <- x[, .(N_total_SMR_tests=.N, N_valid_p_SMR=sum(p_SMR_valid),
    N_FDR_lt_0.05=sum(FDR_recalc<.05,na.rm=TRUE),
    N_HEIDI_no_heterogeneity=sum(FDR_recalc<.05 & HEIDI_assessment=="No heterogeneity detected",na.rm=TRUE),
    N_HEIDI_heterogeneity=sum(FDR_recalc<.05 & HEIDI_assessment=="Heterogeneity detected",na.rm=TRUE),
    N_HEIDI_unavailable=sum(FDR_recalc<.05 & HEIDI_assessment=="HEIDI unavailable",na.rm=TRUE)),by=Tissue]
stopifnot(sum(summary$N_FDR_lt_0.05)==nrow(sig))
write_tsv(x,"03_all16_combined/GTEx_nonbrain_16tissues_SMR_all_results_withFDR.tsv")
write_tsv(summary,"03_all16_combined/GTEx_nonbrain_16tissues_summary.tsv")
write_tsv(sig,"04_FDR_significant/GTEx_nonbrain_16tissues_SMR_FDRsig.tsv")
genes <- sig[,.(N_associations=.N,N_tissues=uniqueN(Tissue),Tissues=paste(sort(unique(Tissue)),collapse="; "),Probes=paste(unique(probeID),collapse="; ")),by=Gene]
write_tsv(genes,"04_FDR_significant/significant_genes_by_tissue.tsv")
print(as.data.frame(summary),row.names=FALSE)
cat("Significant:",nrow(sig),"; genes:",uniqueN(sig$Gene),"; tissues:",uniqueN(sig$Tissue),"\n")
