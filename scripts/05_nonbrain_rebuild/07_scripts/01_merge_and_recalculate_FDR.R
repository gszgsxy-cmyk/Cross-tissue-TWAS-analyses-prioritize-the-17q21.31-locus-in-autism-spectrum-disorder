# Rebuild only derived output files. Original CSV/SMR inputs remain read-only.
source("00_config.R")
audit <- list(); duplicates <- list(); invalids <- list()
required <- c("probeID","ProbeChr","Gene","Probe_bp","topSNP","A1","A2","Freq","b_SMR","se_SMR","p_SMR","p_HEIDI","nsnp_HEIDI")
for (t in tissues) {
  if (t %in% new_tissues) {
    x <- rbindlist(lapply(1:22, function(ch) {
      f <- file.path(new_root,t,paste0("GTEX_",t,"_chr",ch,".smr"))
      stopifnot(file.exists(f))
      z <- fread(f,na.strings=c("NA","NaN",""))
      stopifnot(all(required %in% names(z)),all(z$ProbeChr==ch))
      z[, `:=`(Tissue=t,Chromosome_source=paste0("chr",ch),Chromosome_source_basis="SMR chromosome filename",Source_file=f,Source_row=seq_len(.N))]
      z
    }))
  } else {
    f <- old_map[Tissue==t,Path];x <- fread(f,na.strings=c("NA","NaN",""))
    stopifnot(all(required %in% names(x)),setequal(unique(x$ProbeChr),1:22))
    if ("qtl_name" %in% names(x)) stopifnot(all(x$qtl_name==paste0("eQTL_GTEx_",t)))
    x[, `:=`(Tissue=t,Chromosome_source=paste0("chr",ProbeChr),Chromosome_source_basis="ProbeChr in historical merged CSV",Source_file=f,Source_row=seq_len(.N))]
  }
  x[, p_SMR_valid := is.finite(p_SMR) & p_SMR>=0 & p_SMR<=1]
  x[, FDR_recalc := NA_real_]
  x[p_SMR_valid==TRUE,FDR_recalc:=p.adjust(p_SMR,method="BH")]
  stopifnot(all(is.na(x$p_HEIDI)|(is.finite(x$p_HEIDI)&x$p_HEIDI>=0&x$p_HEIDI<=1)))
  x[, HEIDI_assessment:=fifelse(is.na(p_HEIDI),"HEIDI unavailable",fifelse(p_HEIDI>.05,"No heterogeneity detected","Heterogeneity detected"))]
  x[, Record_ID:=paste(Tissue,probeID,Source_row,sep="|")]
  dup <- x[duplicated(probeID)|duplicated(probeID,fromLast=TRUE)]
  duplicates[[t]]<-dup;invalids[[t]]<-x[p_SMR_valid==FALSE]
  fd <- if("fdr_dan" %in% names(x)) max(abs(x$fdr_dan-x$FDR_recalc),na.rm=TRUE) else NA_real_
  audit[[t]]<-data.table(Tissue=t,N_rows=nrow(x),N_chr=uniqueN(x$ProbeChr),N_valid_p=sum(x$p_SMR_valid),N_zero_p=sum(x$p_SMR==0,na.rm=TRUE),N_duplicate_probe_rows=nrow(dup),N_duplicate_gene_tissue_rows=sum(duplicated(x$Gene)),N_missing_Gene=sum(is.na(x$Gene)|x$Gene==""),N_missing_HEIDI=sum(is.na(x$p_HEIDI)),N_missing_nsnp_HEIDI=sum(is.na(x$nsnp_HEIDI)),Max_difference_historical_FDR=fd)
  write_tsv(x,paste0("02_per_tissue_merged/",t,"_allchr_SMR_withFDR.tsv"))
}
write_tsv(rbindlist(audit),"01_QC/per_tissue_QC.tsv")
write_tsv(rbindlist(duplicates,fill=TRUE),"01_QC/duplicate_probe_records.tsv")
write_tsv(rbindlist(invalids,fill=TRUE),"01_QC/invalid_p_SMR_records.tsv")
capture.output(sessionInfo(),file=file.path(root,"01_QC/R_sessionInfo.txt"))
