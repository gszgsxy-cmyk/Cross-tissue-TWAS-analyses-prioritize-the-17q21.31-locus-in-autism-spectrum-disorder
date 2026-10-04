# Hannon fetal-brain sensitivity analysis: postprocess existing SMR output only.
# Usage: Rscript 12_analysis_script.R [output_directory]
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
library(data.table)
options(digits=17)
arg=commandArgs(TRUE);out=if(length(arg))arg[1] else '__SMR_ROOT__/taierbrain_ASD_sensitivity'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src='__SMR_ROOT__/taierbrain_ASD_sensitivity/04_Hannon_FetalBrain_ASD.smr'
num=function(x)suppressWarnings(as.numeric(x))
miss=function(x)is.na(x)|trimws(x)%in%c('','NA','NaN')
wr=function(d,n)fwrite(d,file.path(out,n),sep='\t',na='NA',quote='auto')
a=fread(src,colClasses='character',na.strings=NULL)
stopifnot(nrow(a)>0,!anyDuplicated(a$probeID))
a[,p_SMR_numeric:=num(p_SMR)]
a[,Valid_p_SMR:=is.finite(p_SMR_numeric)&p_SMR_numeric>=0&p_SMR_numeric<=1]
a[,FDR_FetalBrain:=NA_real_];a[Valid_p_SMR==TRUE,FDR_FetalBrain:=p.adjust(p_SMR_numeric,'BH')]
annroot='__USER_HOME__/Desktop/smr/甲基化注释'
short=fread(file.path(annroot,'注释.csv'),colClasses='character',na.strings=NULL)
# Manifest [Controls] footer can produce an expected field-count warning.
gpl=fread(file.path(annroot,'GPL13534_HumanMethylation450_15017482_v.1.1.csv'),skip='IlmnID,Name,',colClasses='character',na.strings=NULL)
stopifnot(!anyDuplicated(short$probeID),!anyDuplicated(gpl$IlmnID))
standardize=function(ann,gene){
 x=trimws(unlist(strsplit(c(ann,gene)[!miss(c(ann,gene))],';',fixed=TRUE)));x=x[!miss(x)];x[x=='KIAA1267']='KANSL1'
 if(length(x))paste(unique(x),collapse=';') else NA_character_
}
annotate=function(x){
 z=copy(x);if(!'Gene'%in%names(z))z[,Gene:=NA_character_]
 i=match(z$probeID,short$probeID);j=match(z$probeID,gpl$IlmnID)
 for(n in setdiff(names(short),'probeID'))z[,(paste0('Annotation_',n)):=short[[n]][i]]
 for(n in names(gpl))z[,(paste0('GPL_',n)):=gpl[[n]][j]]
 z[,Original_annotation_gene:=Annotation_UCSC_RefGene_Name]
 gene=z$Original_annotation_gene;fallback=miss(gene);gene[fallback]=z$GPL_UCSC_RefGene_Name[fallback]
 z[,Standardized_gene:=mapply(standardize,gene,Gene,USE.NAMES=FALSE)]
 z[,`:=`(hg19_CHR=GPL_CHR,hg19_BP=num(GPL_MAPINFO))]
 z[is.na(GPL_Genome_Build)|GPL_Genome_Build!='37',`:=`(hg19_CHR=NA_character_,hg19_BP=NA_real_)]
 z[,In_17q21.31_core:=ifelse(is.na(hg19_CHR)|is.na(hg19_BP),NA,hg19_CHR=='17'&hg19_BP>=43492082&hg19_BP<=44913631)]
 z[,Region_assignment_source:='GPL Genome_Build 37 CHR/MAPINFO']
 if('ProbeChr'%in%names(z)){
  outside=is.na(z$In_17q21.31_core)&!miss(z$ProbeChr)&z$ProbeChr!='17'&z$ProbeChr==z$GPL_CHR
  outside[is.na(outside)]=FALSE
  z[outside,`:=`(In_17q21.31_core=FALSE,Region_assignment_source='Concordant non-chr17 chromosome in original resource and GPL; hg19 position not asserted')]
 }
 z
}
a=annotate(a)
wr(a[is.na(hg19_BP),.(probeID,ProbeChr,Probe_bp,GPL_Genome_Build,GPL_CHR,GPL_MAPINFO,In_17q21.31_core,Region_assignment_source)],'coordinate_annotation_exceptions.tsv')
a[,HEIDI_assessment:=ifelse(miss(p_HEIDI),'Unavailable',ifelse(num(p_HEIDI)>.05,'No heterogeneity detected','Heterogeneity detected'))]
stopifnot(all(miss(a$p_HEIDI)|(is.finite(num(a$p_HEIDI))&num(a$p_HEIDI)>=0&num(a$p_HEIDI)<=1)))
wr(a,'05_Hannon_FetalBrain_all_withFDR.tsv')
s=a[!is.na(FDR_FetalBrain)&FDR_FetalBrain<.05];setorder(s,FDR_FetalBrain,p_SMR_numeric,probeID)
front=c('probeID','Standardized_gene','In_17q21.31_core','FDR_FetalBrain','p_SMR','b_SMR','se_SMR','p_HEIDI','nsnp_HEIDI','HEIDI_assessment','Original_annotation_gene','ProbeChr','Probe_bp')
setcolorder(s,c(front,setdiff(names(s),front)));wr(s,'06_Hannon_FetalBrain_FDRsig_annotated.tsv')
brain=fread('__USER_HOME__/Desktop/TWAS/2026-06-18 修稿/数据源附表/smr-new/有意义-eBrain-mqtl-fdr注释.csv')
eur=fread('__USER_HOME__/Desktop/TWAS/2026-09-22_Hatton_EUR_mQTL_rebuild/02_Hatton_EUR_FDRsig_annotated.tsv')
stopifnot(nrow(brain)==21,nrow(eur)==30)
related=function(z,g)sum(vapply(strsplit(as.character(ifelse(is.na(z$Standardized_gene),'',z$Standardized_gene)),';',fixed=TRUE),function(v)g%in%v,logical(1)))
compare=function(other,label,file){
 ids=sort(union(s$probeID,other$probeID));d=data.table(probeID=ids,Fetal_significant=ids%in%s$probeID,Comparison_significant=ids%in%other$probeID)
 d[,Category:=ifelse(Fetal_significant&Comparison_significant,'Overlap',ifelse(Fetal_significant,'Fetal-only',paste0(label,'-only')))]
 d[,`:=`(Fetal_tested=probeID%in%a[Valid_p_SMR==TRUE]$probeID,Fetal_p_SMR=num(a$p_SMR[match(probeID,a$probeID)]),Fetal_FDR=a$FDR_FetalBrain[match(probeID,a$probeID)],Comparison_p_SMR=other$p_SMR[match(probeID,other$probeID)],Comparison_FDR=other[[if(label=='BrainMeta')'fdr_dan' else 'FDR_EUR']][match(probeID,other$probeID)])]
 d=annotate(d);wr(d,file)
 ov=d[Category=='Overlap']
 data.table(Comparison=label,Overlap=nrow(ov),Fetal_only=sum(d$Category=='Fetal-only'),Comparison_only=sum(!d$Fetal_significant),Overlap_in_core=sum(ov$In_17q21.31_core,na.rm=TRUE),Overlap_MAPT=related(ov,'MAPT'),Overlap_KANSL1=related(ov,'KANSL1'),Overlap_CRHR1=related(ov,'CRHR1'))
}
cmp=rbind(compare(brain,'BrainMeta','08_FetalBrain_vs_BrainMeta_overlap.tsv'),compare(eur,'HattonEUR','09_FetalBrain_vs_HattonEUR_overlap.tsv'))
wr(cmp,'cross_dataset_counts.tsv')
v=a[Valid_p_SMR==TRUE];stopifnot(!anyNA(v$In_17q21.31_core))
v[,Significant:=!is.na(FDR_FetalBrain)&FDR_FetalBrain<.05]
tab=matrix(c(sum(v$In_17q21.31_core&v$Significant),sum(v$In_17q21.31_core&!v$Significant),sum(!v$In_17q21.31_core&v$Significant),sum(!v$In_17q21.31_core&!v$Significant)),2,byrow=TRUE,dimnames=list(c('Inside_core','Outside_core'),c('Significant','Not_significant')))
wr(data.table(Region=rownames(tab),tab),'regional_2x2.tsv')
if(all(rowSums(tab)>0)&&all(colSums(tab)>0)){
 ft=fisher.test(tab,alternative='two.sided');N=sum(tab);K=sum(tab[,1]);n=sum(tab[1,]);x=tab[1,1];support=seq(max(0,n-(N-K)),min(n,K));base=dhyper(support,K,N-K,n,log=TRUE)
 prob=function(t){z=base+support*t;z=exp(z-max(z));z/sum(z)}
 root=function(f)exp(uniroot(f,c(-50,50),tol=1e-11)$root)
 est=if(x==min(support))0 else if(x==max(support))Inf else root(function(t)sum(support*prob(t))-x)
 lo=if(x==min(support))0 else root(function(t)sum(prob(t)[support>=x])-.025)
 hi=if(x==max(support))Inf else root(function(t)sum(prob(t)[support<=x])-.025)
 enrich=data.table(Test='two-sided Fisher exact',OR=est,CI_lower=lo,CI_upper=hi,P=ft$p.value)
}else enrich=data.table(Test='Not estimable: zero marginal count',OR=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,P=NA_real_)
wr(enrich,'10_FetalBrain_17q21.31_enrichment.tsv')
mapping=fread('__SMR_ROOT__/taierbrain_ASD_sensitivity/02_fetalbrain_SNP_mapping.tsv')
nprobe=nrow(fread('__SMR_ROOT__/taierbrain/FB_Brain_2.epi',header=FALSE))
summary=c(sprintf('Total fetal-brain probes in resource: %d',nprobe),sprintf('Total mQTL SNPs: %d',nrow(mapping)),sprintf('SNPs mapped to 1000G rsIDs: %d',sum(mapping$mapping_status=='mapped')),sprintf('Mapping percentage: %.6f%%',100*mean(mapping$mapping_status=='mapped')),sprintf('Total valid SMR tests: %d',nrow(v)),sprintf('FDR-significant probes: %d',nrow(s)),sprintf('17q21.31 tested probes: %d',sum(v$In_17q21.31_core)),sprintf('17q21.31 significant probes: %d',sum(s$In_17q21.31_core)),sprintf('Proportion of significant probes in core: %s',if(nrow(s))sprintf('%.6f%%',100*mean(s$In_17q21.31_core)) else 'NA (no significant probes)'),sprintf('MAPT-related: %d',related(s,'MAPT')),sprintf('KANSL1/KIAA1267-related: %d',related(s,'KANSL1')),sprintf('CRHR1-related: %d',related(s,'CRHR1')),sprintf('HEIDI no heterogeneity: %d',sum(s$HEIDI_assessment=='No heterogeneity detected')),sprintf('HEIDI heterogeneity: %d',sum(s$HEIDI_assessment=='Heterogeneity detected')),sprintf('HEIDI unavailable: %d',sum(s$HEIDI_assessment=='Unavailable')),sprintf('Overlap with Brain-mMeta: %d',cmp[Comparison=='BrainMeta']$Overlap),sprintf('Overlap with Hatton EUR: %d',cmp[Comparison=='HattonEUR']$Overlap),'',capture.output(print(cmp)),capture.output(print(tab)),capture.output(print(enrich)),'','METHODS / INTERPRETATION','One genome-wide SMR/HEIDI run, SMR 1.3.1, existing ASD GWAS, 1000G EUR LD, HEIDI min 3 / max 20, 10 threads; other defaults as in archived Brain-mMeta command.','Fetal ESI has no allele-frequency column. SMR uses reference/GWAS frequency QC, but fetal-versus-GWAS frequency concordance cannot be assessed.','One BH correction across all finite p_SMR in [0,1]. No regional, gene, chromosome or HEIDI prefilter. No nsnp_HEIDI >=10 criterion.','Mapping preserves ESI row order, all allele coding, probe order and original BESD bytes. Only unique position and identical unordered allele pairs are mapped; no complement inference. Failed mappings excluded explicitly.','Original annotation columns and SMR fields preserved. Historical KIAA1267 standardized only in Standardized_gene. Region assignment uses GPL hg19 coordinates.','Hannon fetal-brain data: n=166; hg19; resource documentation states mQTL P<=1e-10. However, SMR query reports P=9.9382e-10 for rs28444699/cg08128007. This discrepancy is preserved in QC; no additional P-based filtering was applied. The source .summary also lists generic extraction thresholds.','Developmental-brain sensitivity analysis. Hannon fetal brain is a component of Brain-mMeta; overlaps are descriptive, not independent replication. No overlap enrichment test performed.','Hatton EUR comparison is developmental brain versus whole blood, using fixed finalized 30-probe set.','Regional enrichment is probe-level and does not model correlation among CpGs; P is nominal. Exact conditional OR/CI solved on log-OR scale for numerical precision.','Resource documentation: https://yanglab.westlake.edu.cn/software/smr/','Primary article: https://www.nature.com/articles/nn.4182')
summary=c(summary,sprintf('Coordinate exceptions: %d tested probes lack GPL hg19 positions. Concordant original/GPL non-chr17 chromosomes classify these as outside the core without asserting hg19 positions. See coordinate_annotation_exceptions.tsv.',sum(is.na(v$hg19_BP))))
writeLines(summary,file.path(out,'11_Hannon_FetalBrain_summary.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat(paste(summary[1:17],collapse='\n'),'\n')



