# Existing SMR outputs only. No SMR execution, no significance redefinition.
# Usage: Rscript enrichment.R [output_directory]
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
library(data.table)
options(digits=17)
args=commandArgs(TRUE)
out=if(length(args)) args[1] else 'outputs/2026-09-22_Hatton_EUR_mQTL_rebuild/Enrichment'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base='__USER_HOME__/Desktop/TWAS/2026-09-22_Hatton_EUR_mQTL_rebuild'
eurfile=file.path(base,'01_Hatton_EUR_all_chr_merged_withFDR.tsv')
bsfile='__USER_HOME__/Desktop/TWAS/2026-06-18 修稿/数据源附表/smr-new/有意义-eBrain-mqtl-fdr注释.csv'
bf='__USER_HOME__/Desktop/TWAS/结果存档/smr-new/eBrain-mqtl-fdr.csv'
gplfile='__USER_HOME__/Desktop/smr/甲基化注释/GPL13534_HumanMethylation450_15017482_v.1.1.csv'
brainfiles=file.path('__SMR_ROOT__/Brain-mMeta-asd',sprintf('Brain-mMeta-asd_chr%d.smr',1:22))
wr=function(x,n)fwrite(x,file.path(out,n),sep='\t',na='NA')
# Invert the conditional noncentral hypergeometric distribution on a log-OR
# scale. This avoids absolute-root-tolerance loss in fisher.test's very large
# OR confidence bounds (reciprocal-OR roots close to zero).
exact_or=function(m){
 N=sum(m);K=sum(m[,1]);n=sum(m[1,]);x=m[1,1]
 support=seq(max(0,n-(N-K)),min(n,K));base=dhyper(support,K,N-K,n,log=TRUE)
 prob=function(logor){z=base+support*logor;z=exp(z-max(z));z/sum(z)}
 root=function(f)exp(uniroot(f,c(-40,40),tol=1e-11)$root)
 c(OR=root(function(t)sum(support*prob(t))-x),lower=root(function(t)sum(prob(t)[support>=x])-.025),upper=root(function(t)sum(prob(t)[support<=x])-.025))
}
eur=fread(eurfile);eur=eur[is.finite(p_SMR)&p_SMR>=0&p_SMR<=1]
brainall=rbindlist(lapply(brainfiles,fread)); brain=brainall[is.finite(p_SMR)&p_SMR>=0&p_SMR<=1]
bs=fread(bsfile);archive=fread(bf)
stopifnot(nrow(eur)==137310,!anyDuplicated(eur$probeID),!anyDuplicated(brain$probeID),nrow(bs)==21)
stopifnot(setequal(brainall$probeID,archive$probeID),setequal(archive[fdr_dan<.05]$probeID,bs$probeID))
stopifnot(max(abs(brainall$p_SMR-archive$p_SMR[match(brainall$probeID,archive$probeID)]))<1e-14)
es=eur[FDR_EUR<.05]$probeID;stopifnot(length(es)==30)
# The manifest's [Controls] footer is not part of the assay table.
gpl=fread(gplfile,skip='IlmnID,Name,',select=c('IlmnID','Genome_Build','CHR','MAPINFO'),colClasses='character')
stopifnot(!anyDuplicated(gpl$IlmnID))
i=match(eur$probeID,gpl$IlmnID)
eur[,`:=`(hg19_CHR=gpl$CHR[i],hg19_BP=as.numeric(gpl$MAPINFO[i]),Genome_Build=gpl$Genome_Build[i])]
print(eur[,.(N=.N,Missing_chr=sum(is.na(hg19_CHR)),Missing_bp=sum(is.na(hg19_BP)),Non37=sum(is.na(Genome_Build)|Genome_Build!='37'))])
stopifnot(all(eur$Genome_Build=='37'),all(!is.na(eur$hg19_CHR)),all(is.finite(eur$hg19_BP)))
eur[,In_core:=hg19_CHR=='17'&hg19_BP>=43492082&hg19_BP<=44913631]
eur[,Significant:=probeID%in%es]
wr(eur[,.(probeID,p_SMR,FDR_EUR,Significant,ProbeChr,Probe_bp,hg19_CHR,hg19_BP,In_core)],'01_EUR_background_region_classification.tsv')
region=matrix(c(sum(eur$Significant&eur$In_core),sum(eur$Significant&!eur$In_core),sum(!eur$Significant&eur$In_core),sum(!eur$Significant&!eur$In_core)),2,byrow=TRUE,dimnames=list(c('EUR_significant','EUR_not_significant'),c('Inside_core','Outside_core')))
fr=fisher.test(region,alternative='two.sided',conf.level=.95)
er=exact_or(region);fr$estimate[]=er[1];fr$conf.int[]=er[2:3]
wr(data.table(Group=rownames(region),region),'02_region_2x2.tsv')
universe=sort(intersect(eur$probeID,brain$probeID))
u=data.table(probeID=universe,Brain_significant=universe%in%bs$probeID,EUR_significant=universe%in%es)
u[,`:=`(Brain_p_SMR=brain$p_SMR[match(probeID,brain$probeID)],EUR_p_SMR=eur$p_SMR[match(probeID,eur$probeID)])]
wr(u,'03_shared_testable_universe.tsv')
ov=matrix(c(sum(u$Brain_significant&u$EUR_significant),sum(u$Brain_significant&!u$EUR_significant),sum(!u$Brain_significant&u$EUR_significant),sum(!u$Brain_significant&!u$EUR_significant)),2,byrow=TRUE,dimnames=list(c('Brain_significant','Brain_not_significant'),c('EUR_significant','EUR_not_significant')))
fo=fisher.test(ov,alternative='greater');fo2=fisher.test(ov,alternative='two.sided')
eo=exact_or(ov);fo$estimate[]=eo[1];fo2$estimate[]=eo[1];fo2$conf.int[]=eo[2:3]
wr(data.table(Group=rownames(ov),ov),'04_overlap_2x2.tsv')
excluded=rbindlist(list(data.table(Dataset='Brain-mMeta',probeID=setdiff(bs$probeID,universe)),data.table(Dataset='EUR',probeID=setdiff(es,universe))))
wr(excluded,'05_significant_probes_outside_shared_universe.tsv')
expected=sum(u$Brain_significant)*sum(u$EUR_significant)/nrow(u)
results=data.table(Analysis=c('17q21.31 core enrichment','Brain-mMeta EUR overlap enrichment'),Alternative=c('two.sided','greater'),Universe_N=c(nrow(eur),nrow(u)),Observed=c(region[1,1],ov[1,1]),Expected=c(length(es)*sum(eur$In_core)/nrow(eur),expected),Conditional_odds_ratio=c(unname(fr$estimate),unname(fo$estimate)),OR_95CI_lower_two_sided=c(fr$conf.int[1],fo2$conf.int[1]),OR_95CI_upper_two_sided=c(fr$conf.int[2],fo2$conf.int[2]),P_value=c(fr$p.value,fo$p.value),Observed_expected_ratio=c(region[1,1]/(length(es)*sum(eur$In_core)/nrow(eur)),ov[1,1]/expected))
wr(results,'06_enrichment_statistics.tsv')
summary=c('EXISTING OUTPUTS ONLY: no SMR rerun; original significance sets unchanged.',sprintf('EUR valid tested probes: %d; Brain valid tested probes: %d; Brain raw rows: %d',nrow(eur),nrow(brain),nrow(brainall)),sprintf('All %d EUR probes mapped to hg19 GPL coordinates; missing coordinates: 0.',nrow(eur)),'','REGION ENRICHMENT: two-sided Fisher exact test',capture.output(print(region)),capture.output(print(fr)),'','OVERLAP ENRICHMENT: one-sided Fisher exact test (greater)',sprintf('Shared-testable universe: %d',nrow(u)),sprintf('Original significant sets: Brain=21, EUR=30; within shared universe: Brain=%d, EUR=%d',sum(u$Brain_significant),sum(u$EUR_significant)),sprintf('Observed overlap: %d; expected overlap: %.12g; observed/expected: %.12g',ov[1,1],expected,ov[1,1]/expected),capture.output(print(ov)),sprintf('Primary one-sided Fisher P: %.12g',fo$p.value),sprintf('Two-sided 95%% CI for overlap odds ratio: %.12g to %.12g',fo2$conf.int[1],fo2$conf.int[2]),sprintf('Supplementary two-sided overlap P: %.12g',fo2$p.value),'Significant probes excluded only from overlap universe (not reclassified):',capture.output(print(excluded)),'','METHODS','Testable means present in existing SMR output with finite p_SMR in [0,1]; exact CpG IDs are used, no HEIDI filter.','Original EUR FDR<0.05 set and finalized Brain-mMeta 21-probe set are retained. Overlap test margins use only significant probes present in the shared-testable universe.','Region classification uses GPL Genome_Build 37 CHR/MAPINFO, identical to prior rebuild; inclusive interval chr17:43492082-44913631.','OR estimates are conditional maximum-likelihood estimates; confidence intervals are exact conditional two-sided 95% intervals, solved on a log-OR scale. Fisher P values use R fisher.test.','The overlap primary test is one-sided (greater), matching the enrichment question. Region test is two-sided as requested. P values are nominal, not adjusted across these two analyses.','CpG correlations and regional clustering are not modeled by Fisher tests. These are probe-level enrichments; overlap does not establish independent biological replication.','Full archived Brain-mMeta probe IDs and p_SMR agree with chr1-22 raw outputs; archived significant IDs agree with the finalized 21-probe set.')
writeLines(summary,file.path(out,'07_enrichment_summary.txt'))
wr(data.table(Source=c(eurfile,bsfile,bf,gplfile,brainfiles),MD5=unname(tools::md5sum(c(eurfile,bsfile,bf,gplfile,brainfiles)))),'08_input_checksums.tsv')
capture.output(sessionInfo(),file=file.path(out,'09_sessionInfo.txt'))
print(results);cat(paste(summary[1:4],collapse='\n'),'\n');print(excluded);print(region);print(ov)
write('Numerical note: conditional OR and exact two-sided CI were solved on the log-OR scale with tolerance 1e-11 to avoid loss of precision in default fisher.test reciprocal-root CI bounds for very large odds ratios. Fisher P values are unchanged.',file=file.path(out,'07_enrichment_summary.txt'),append=TRUE)

