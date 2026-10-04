source("00_config.R")
x <- fread(file.path(root,"04_FDR_significant/GTEx_nonbrain_16tissues_SMR_FDRsig.tsv"))
stopifnot(all(is.finite(x$b_SMR)),all(is.finite(x$se_SMR)),all(x$se_SMR>0))
x[, `:=`(CI_lower=b_SMR-1.96*se_SMR,CI_upper=b_SMR+1.96*se_SMR,Probe=probeID)]
x[, Plot_row := seq_len(.N)]
setcolorder(x,c("Plot_row","Gene","Tissue","Probe","ProbeChr","Probe_bp","topSNP","b_SMR","se_SMR","CI_lower","CI_upper","p_SMR","FDR_recalc","p_HEIDI","nsnp_HEIDI","HEIDI_assessment"))
write_tsv(x,"05_Figure4C/Figure4C_GTEx_nonbrain_SMR_forest_data.tsv")
