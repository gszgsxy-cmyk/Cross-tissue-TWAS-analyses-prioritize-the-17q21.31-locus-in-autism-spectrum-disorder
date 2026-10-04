Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
options(stringsAsFactors=FALSE,warn=1)
library(data.table);library(coloc)
stopifnot(as.character(getRversion())=='4.4.1',as.character(packageVersion('coloc'))=='5.2.3')
out<-'outputs/BrainMeta_fullcis_coloc_final_20260924'
stopifnot(!file.exists(file.path(out,'04_summary/01_final_23probe_coloc_results.tsv')))
wt<-function(d,p)fwrite(d,file.path(out,p),sep='\t',na='NA')
lg<-function(...) {cat(format(Sys.time()),..., '\n');flush.console()}
start<-Sys.time();lg('BATCH START')
inv<-fread(file.path(out,'01_inventory/01_final_23probe_inventory.tsv'))
stopifnot(nrow(inv)==23,uniqueN(inv$Probe)==23,all(inv$BrainMeta_source_status%in%c('A','PASS')))
qs<-fread(file.path(out,'05_logs/query_status.tsv'))
stopifnot(nrow(qs)==23,setequal(qs$Probe,inv$Probe))
capture.output(sessionInfo(),file=file.path(out,'05_logs/R_sessionInfo.txt'))
validallele<-function(a)!is.na(a)&a%in%c('A','C','G','T')
validstat<-function(d)is.finite(d$b)&is.finite(d$se)&d$se>0&is.finite(d$se^2)&d$se^2>0
ispal<-function(a,b)paste0(a,b)%in%c('AT','TA','CG','GC')
comp<-function(a)chartr('ACGT','TGCA',a)
near<-function(a,b) (is.na(a)&is.na(b)) | (is.finite(a)&is.finite(b)&abs(a-b)<=1e-12+1e-6*pmax(abs(a),abs(b)))
direction<-function(a,b,x,y)fifelse(a==x&b==y,1L,fifelse(a==y&b==x,-1L,fifelse(!ispal(a,b)&comp(a)==x&comp(b)==y,1L,fifelse(!ispal(a,b)&comp(a)==y&comp(b)==x,-1L,NA_integer_))))
clean<-function(d,src) {
 z<-copy(d);ex<-list()
 record<-function(a,r){if(nrow(a))ex[[length(ex)+1L]]<<-data.table(source=src,row=a$row,SNP=a$SNP,reason=r)}
 # Resolve conflicting identity across all records before selecting valid statistics.
 du<-z[,.(N=.N),by=SNP][N>1,SNP]
 bad<-character()
 for(s in du){a<-z[SNP==s];pairs<-paste(pmin(a$A1,a$A2),pmax(a$A1,a$A2));if(uniqueN(pairs)!=1 || (src=='BrainMeta' && uniqueN(paste(a$CHR,a$BP))!=1))bad<-c(bad,s)}
 record(z[SNP%in%bad],'ambiguous_duplicate_identity');z<-z[!SNP%in%bad]
 reason<-rep(NA_character_,nrow(z));tag<-function(b,l){b[is.na(b)]<-TRUE;reason[is.na(reason)&b]<<-l}
 tag(is.na(z$SNP)|z$SNP%in%c('','.')|grepl('[[:space:]]',z$SNP),'invalid_identifier')
 if(src=='BrainMeta')tag(!is.finite(z$BP)|z$BP<=0|z$BP!=floor(z$BP)|z$CHR!=17,'invalid_position')
 tag(!validallele(z$A1)|!validallele(z$A2)|z$A1==z$A2,'not_biallelic_ACGT_SNP')
 tag(!validstat(z),'invalid_beta_SE_variance')
 if(src=='BrainMeta')tag(!is.finite(z$f)|z$f<=0|z$f>=1,'invalid_BrainMeta_frequency')
 for(r in unique(na.omit(reason)))record(z[which(reason==r)],r)
 z<-z[is.na(reason)]
 dup<-duplicated(z$SNP)|duplicated(z$SNP,fromLast=TRUE);keep<-which(!dup)
 for(ix in split(which(dup),z$SNP[dup])){
  a<-z[ix];d<-ifelse(a$A1==a$A1[1]&a$A2==a$A2[1],1,-1)
  okay<-all(near(a$b*d,a$b[1]))&&all(near(a$se,a$se[1]))
  if(src=='BrainMeta')okay<-okay&&all(near(ifelse(d==-1,1-a$f,a$f),a$f[1]))
  if(okay){keep<-c(keep,ix[1]);record(a[-1],'equivalent_duplicate_collapsed')}else record(a,'ambiguous_duplicate_statistics')
 }
 list(data=z[sort(keep)],excluded=if(length(ex))rbindlist(ex) else data.table(source=character(),row=integer(),SNP=character(),reason=character()))
}
g0<-fread('__SMR_ROOT__/asd-gwas/gwas_chr17.txt')
stopifnot(all(c('SNP','A1','A2','b','se','p','freq','n')%in%names(g0)),all(g0$n==46350))
g<-g0[,.(row=.I,SNP,CHR=17,BP=NA_real_,A1=toupper(A1),A2=toupper(A2),b=as.numeric(b),se=as.numeric(se),f=as.numeric(freq),p=as.numeric(p))]
gc<-clean(g,'ASD');gw<-gc$data
stopifnot(nrow(gw)>0,!anyDuplicated(gw$SNP))
wt(gc$excluded,'05_logs/ASD_excluded_records.tsv')
wt(data.table(N_ASD_raw=nrow(g),N_ASD_valid=nrow(gw),ASD_frequency_source='external 1000G EUR reference; not observed cohort EAF',ASD_frequency_used_in_coloc=FALSE),'05_logs/ASD_global_QC.tsv')
results<-list();qclist<-list();sourceqc<-list();fails<-list();warnings<-list();calls<-0L
for(k in seq_len(nrow(inv))) {
 gene<-inv$Gene[k];probe<-inv$Probe[k];pb<-inv$Probe_BP[k]
 lg('PROBE START',gene,probe)
 dir<-file.path('03_probe_results',paste0(gene,'__',probe));dir.create(file.path(out,dir),recursive=TRUE,showWarnings=FALSE)
 warnrows<-list();stage<-'query';t0<-Sys.time()
 capturewarn<-function(expr){withCallingHandlers(expr,warning=function(w){warnrows[[length(warnrows)+1L]]<<-data.table(Gene=gene,Probe=probe,stage=stage,warning=conditionMessage(w));lg('WARNING',probe,stage,conditionMessage(w));invokeRestart('muffleWarning')})}
 tryCatch({
  stopifnot(qs[Probe==probe,Exit_code]==0)
  e0<-fread(file.path(out,'02_BrainMeta_queries',paste0('BrainMeta_',probe,'.txt')))
  stopifnot(nrow(e0)>0,all(e0$Probe==probe),all(e0$Probe_bp==pb))
  # Completeness counts are provenance verification, never a filtering target.
  if(nrow(e0)!=inv$Expected_stored_variants[k])stop('Fresh extraction count conflicts with source completeness inventory')
  stage<-'basic_QC'
  e<-e0[,.(row=.I,SNP,CHR=as.numeric(Chr),BP=as.numeric(BP),A1=toupper(A1),A2=toupper(A2),b=as.numeric(b),se=as.numeric(SE),f=as.numeric(Freq),p=as.numeric(p))]
  ec<-clean(e,'BrainMeta');eq<-ec$data
  wt(ec$excluded,file.path(dir,'excluded_records.tsv'))
  pv<-e$p[is.finite(e$p)&e$p>=0&e$p<=1]
  sq<-data.table(Gene=gene,Probe=probe,N_BrainMeta_raw=nrow(e),N_BrainMeta_valid=nrow(eq),Min_BP=min(e$BP),Max_BP=max(e$BP),Span_bp=max(e$BP)-min(e$BP)+1,Min_P=min(pv),Max_P=max(pv),N_P_gt_0.05=sum(pv>.05),N_P_gt_0.5=sum(pv>.5),N_P_gt_0.9=sum(pv>.9),N_P_zero=sum(pv==0),N_invalid_P=nrow(e)-length(pv))
  sourceqc[[probe]]<-sq
  stage<-'matching'
  h<-merge(gw[,.(SNP,A1,A2,b,se,f,p)],eq[,.(SNP,CHR,BP,A1,A2,b,se,f,p)],by='SNP',suffixes=c('_ASD','_BrainMeta'))
  h[,palindromic:=ispal(A1_ASD,A2_ASD)&ispal(A1_BrainMeta,A2_BrainMeta)]
  h[,d:=direction(A1_BrainMeta,A2_BrainMeta,A1_ASD,A2_ASD)]
  h[,alignment:=fifelse(is.na(d),'mismatch',fifelse(palindromic,'palindromic',fifelse(A1_ASD==A1_BrainMeta&A2_ASD==A2_BrainMeta,'direct',fifelse(A1_ASD==A2_BrainMeta&A2_ASD==A1_BrainMeta,'swapped',fifelse(d==1,'complement_direct','complement_swapped')))))]
  nshared<-nrow(h);ct<-table(factor(h$alignment,levels=c('direct','swapped','complement_direct','complement_swapped','palindromic','mismatch')))
  wt(h[alignment=='mismatch'],file.path(dir,'allele_mismatches.tsv'))
  h<-h[alignment!='mismatch'];setorder(h,BP,SNP)
  stopifnot(nrow(h)>0,!anyDuplicated(h$SNP))
  h[,`:=`(Beta_gwas=b_ASD,SE_gwas=se_ASD,Beta_eqtl=b_BrainMeta,SE_eqtl=se_BrainMeta,Freq_eqtl=f_BrainMeta,Beta_eqtl_label_aligned=b_BrainMeta*d,MAF_BrainMeta=pmin(f_BrainMeta,1-f_BrainMeta),ASD_position_source='BrainMeta position for matched rsID',ASD_frequency_source='external 1000G EUR reference; not observed cohort EAF')]
  h[,freq_diff_QC_only:=abs(f_ASD-ifelse(d==-1,1-f_BrainMeta,f_BrainMeta))]
  gaps<-diff(sort(unique(h$BP)));q<-function(pr)if(length(gaps))unname(quantile(gaps,pr))else NA_real_
  amb_eq<-uniqueN(ec$excluded[grepl('^ambiguous_duplicate',reason),SNP]);amb_g<-uniqueN(gc$excluded[grepl('^ambiguous_duplicate',reason)&SNP%in%e$SNP,SNP])
  qc<-data.table(Gene=gene,Probe=probe,N_BrainMeta_raw=nrow(e),N_BrainMeta_valid=nrow(eq),N_shared_rsID=nshared,N_shared_rsID_before_allele_QC=nshared,N_direct=as.integer(ct[1]),N_swapped=as.integer(ct[2]),N_complement_direct=as.integer(ct[3]),N_complement_swapped=as.integer(ct[4]),N_palindromic_retained=sum(h$palindromic),N_mismatch_removed=as.integer(ct[6]),N_allele_mismatch_removed=as.integer(ct[6]),N_ambiguous_duplicate_removed=amb_eq+amb_g,N_BrainMeta_ambiguous_duplicate_IDs=amb_eq,N_ASD_ambiguous_duplicate_IDs_in_probe=amb_g,N_allele_compatible=nrow(h),N_final_shared=nrow(h),Min_BP_final=min(h$BP),Max_BP_final=max(h$BP),Median_gap_bp=q(.5),P95_gap_bp=q(.95),P99_gap_bp=q(.99),Max_gap_bp=if(length(gaps))max(gaps)else NA_real_,N_frequency_difference_gt_0.2_retained=sum(h$freq_diff_QC_only>.2,na.rm=TRUE),N_P_zero_retained=sum(h$p_BrainMeta==0,na.rm=TRUE),Status='INPUT_QC_PASS')
  qclist[[probe]]<-qc;wt(qc,file.path(dir,'QC.tsv'));wt(h,file.path(dir,'full_cis_shared_variants.tsv'))
  d1<-list(beta=h$Beta_gwas,varbeta=h$SE_gwas^2,snp=h$SNP,position=h$BP,type='cc',N=46350,s=18381/46350)
  d2<-list(beta=h$Beta_eqtl,varbeta=h$SE_eqtl^2,snp=h$SNP,position=h$BP,type='quant',N=2865,MAF=h$MAF_BrainMeta)
  saveRDS(list(ASD=d1,BrainMeta=d2),file.path(out,dir,'exact_coloc_inputs.rds'))
  stage<-'check_dataset_ASD';capturewarn(check_dataset(d1))
  stage<-'check_dataset_BrainMeta';capturewarn(check_dataset(d2))
  stage<-'coloc.abf';calls<-calls+1L
  fit<-capturewarn(coloc.abf(dataset1=d1,dataset2=d2,p1=1e-4,p2=1e-4,p12=1e-5))
  saveRDS(fit,file.path(out,dir,'coloc_result.rds'))
  pp<-as.numeric(fit$summary[paste0('PP.H',0:4,'.abf')])
  stopifnot(all(is.finite(pp)),all(pp>=0),abs(sum(pp)-1)<1e-10,fit$summary['nsnps']==nrow(h))
  sn<-as.data.table(fit$results);setnames(sn,'snp','SNP');sn[,position:=h$BP[match(SNP,h$SNP)]]
  wt(sn,file.path(dir,'coloc_SNP_results.tsv'));top<-which.max(sn$SNP.PP.H4)
  wr<-if(length(warnrows))rbindlist(warnrows) else data.table(Gene=character(),Probe=character(),stage=character(),warning=character())
  row<-data.table(Gene=gene,Probe=probe,Probe_BP=pb,BrainMeta_min_BP=sq$Min_BP,BrainMeta_max_BP=sq$Max_BP,N_BrainMeta_valid=nrow(eq),N_final_shared=nrow(h),PP_H0=pp[1],PP_H1=pp[2],PP_H2=pp[3],PP_H3=pp[4],PP_H4=pp[5],H4_over_H3=pp[5]/pp[4],H4_given_H3_or_H4=pp[5]/sum(pp[4:5]),Dominant_hypothesis=paste0('H',which.max(pp)-1),Study_defined_strong_coloc=pp[5]>=.8&&pp[5]>3*pp[4],Top_H4_SNP=sn$SNP[top],Top_SNP_PP_H4=sn$SNP.PP.H4[top],Check_dataset_warning=paste(unique(wr[grepl('^check_dataset',stage),warning]),collapse=' | '),Coloc_warning=paste(unique(wr[stage=='coloc.abf',warning]),collapse=' | '),Status='COMPLETED')
  results[[probe]]<-row;wt(row,file.path(dir,'coloc_summary.tsv'));lg('PROBE COMPLETED',gene,'n=',nrow(h),'H4=',pp[5])
 },error=function(e){
  msg<-conditionMessage(e);lg('PROBE FAILED',probe,stage,msg)
  fails[[probe]]<<-data.table(Gene=gene,Probe=probe,stage=stage,error=msg)
  results[[probe]]<<-data.table(Gene=gene,Probe=probe,Probe_BP=pb,Status='FAILED',Error=msg)
 })
 if(length(warnrows))warnings[[probe]]<-rbindlist(warnrows)
 lg('PROBE END',gene,'elapsed_seconds=',as.numeric(difftime(Sys.time(),t0,units='secs')))
}
res<-rbindlist(results,fill=TRUE);qc<-rbindlist(qclist,fill=TRUE)
wt(res,'04_summary/01_final_23probe_coloc_results.tsv');wt(qc,'04_summary/02_final_23probe_harmonization_QC.tsv')
wt(if(length(fails))rbindlist(fails)else data.table(Gene=character(),Probe=character(),stage=character(),error=character()),'04_summary/03_failed_probes.tsv')
wt(res[order(-PP_H4)],'04_summary/04_results_sorted_by_PP_H4.tsv')
wt(rbindlist(sourceqc),'01_inventory/02_full_cis_source_QC.tsv')
wt(if(length(warnings))rbindlist(warnings)else data.table(Gene=character(),Probe=character(),stage=character(),warning=character()),'05_logs/warnings.tsv')
pilotfile<-'outputs/MAPT_AS1_ABF_validation_20260924/01_ABF_results.tsv'
pilot<-fread(pilotfile)[analysis=='full_cis_same_scale_as_previous'];batch<-res[Probe=='ENSG00000264589.4']
stopifnot(nrow(pilot)==1,nrow(batch)==1)
oldval<-c(pilot$n,as.numeric(pilot[,paste0('PP.H',0:4),with=FALSE]));newval<-c(batch$N_final_shared,as.numeric(batch[,paste0('PP_H',0:4),with=FALSE]))
cp<-data.table(metric=c('nsnps',paste0('H',0:4)),pilot_value=oldval,batch_value=newval,absolute_difference=abs(oldval-newval),tolerance=c(0,rep(1e-12,5)),pilot_file=pilotfile)
cp[,PASS:=is.finite(absolute_difference)&absolute_difference<=tolerance]
wt(cp,'04_summary/06_MAPT_AS1_checkpoint_comparison.tsv')
stopifnot(all(cp$PASS),nrow(res)==23,uniqueN(res$Probe)==23)
done<-res[Status=='COMPLETED'];hi<-done[which.max(PP_H4)]
wt(data.table(N_probes_requested=23,N_probes_completed=nrow(done),N_probes_failed=length(fails),N_coloc_abf_calls=calls,N_strong_coloc_by_study_rule=sum(done$Study_defined_strong_coloc),Highest_PP_H4=hi$PP_H4,Gene_of_highest_PP_H4=hi$Gene,Probe_of_highest_PP_H4=hi$Probe,MAPT_AS1_checkpoint_PASS=all(cp$PASS),start_time=as.character(start),end_time=as.character(Sys.time())),'04_summary/07_batch_metrics.tsv')
lg('BATCH END','completed=',nrow(done),'failed=',length(fails),'coloc_abf_calls=',calls)
