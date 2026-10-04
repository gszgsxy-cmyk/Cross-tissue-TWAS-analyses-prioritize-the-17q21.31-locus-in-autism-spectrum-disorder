library(data.table)
o<-'outputs/BrainMeta_fullcis_coloc_final_20260924'
inv<-fread(file.path(o,'01_inventory/01_final_23probe_inventory.tsv'))
res<-fread(file.path(o,'04_summary/01_final_23probe_coloc_results.tsv'))
qc<-fread(file.path(o,'04_summary/02_final_23probe_harmonization_QC.tsv'))
checks<-list()
stopifnot(nrow(res)==23,uniqueN(res$Probe)==23,setequal(inv$Probe,res$Probe))
for(k in seq_len(nrow(res))){
 r<-res[k];if(r$Status!='COMPLETED')next
 p<-file.path(o,'03_probe_results',paste0(r$Gene,'__',r$Probe))
 h<-fread(file.path(p,'full_cis_shared_variants.tsv'));d<-readRDS(file.path(p,'exact_coloc_inputs.rds'));f<-readRDS(file.path(p,'coloc_result.rds'));sn<-fread(file.path(p,'coloc_SNP_results.tsv'))
 stopifnot(identical(d$ASD$snp,d$BrainMeta$snp),identical(h$SNP,d$ASD$snp),!anyDuplicated(h$SNP),length(d$ASD$snp)==r$N_final_shared)
 stopifnot(!'MAF'%in%names(d$ASD),!'LD'%in%names(d$ASD),!'sdY'%in%names(d$BrainMeta),d$BrainMeta$N==2865,d$ASD$N==46350,d$ASD$s==18381/46350)
 stopifnot(max(abs(d$ASD$beta-h$Beta_gwas))<1e-12,max(abs(d$BrainMeta$beta-h$Beta_eqtl))<1e-12,max(abs(d$BrainMeta$MAF-pmin(h$Freq_eqtl,1-h$Freq_eqtl)))<1e-12)
 stopifnot(max(abs(d$ASD$varbeta-h$SE_gwas^2))<1e-12,max(abs(d$BrainMeta$varbeta-h$SE_eqtl^2))<1e-12,identical(as.numeric(d$ASD$position),as.numeric(h$BP)))
 stopifnot(nrow(sn)==nrow(h),setequal(sn$SNP,h$SNP),all(sn$position==h$BP[match(sn$SNP,h$SNP)]),abs(sum(sn$SNP.PP.H4)-1)<1e-10)
 pp<-as.numeric(f$summary[paste0('PP.H',0:4,'.abf')]);delta<-max(abs(pp-as.numeric(r[,paste0('PP_H',0:4),with=FALSE])))
 stopifnot(delta<1e-12,abs(sum(pp)-1)<1e-10,r$Study_defined_strong_coloc==(pp[5]>=.8&&pp[5]>3*pp[4]))
 q<-qc[Probe==r$Probe]
 stopifnot(q$N_shared_rsID==q$N_final_shared+q$N_mismatch_removed,q$N_final_shared==sum(q[,.(N_direct,N_swapped,N_complement_direct,N_complement_swapped,N_palindromic_retained)]))
 checks[[r$Probe]]<-data.table(Gene=r$Gene,Probe=r$Probe,Exact_inputs_PASS=TRUE,Variant_identity_PASS=TRUE,QC_accounting_PASS=TRUE,Posterior_PASS=TRUE,max_posterior_serialization_delta=delta)
}
fwrite(rbindlist(checks),file.path(o,'05_logs/verification_checks.tsv'),sep='\t')
cat('Verified saved input/output identities and posterior accounting for',length(checks),'completed probes; no new inference calls.\n')
