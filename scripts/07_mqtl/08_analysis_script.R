# Hatton EUR mQTL-ASD result rebuild. Reads existing SMR outputs only.
# No SMR binary, association test, or HEIDI test is invoked.
# Usage: Rscript 08_analysis_script.R [output_directory]
# Requires R >= 4.1 and data.table. TSV/QC outputs are produced by this script.
# The two XLSX files are formatted from these TSVs by 09_export_workbooks.mjs.
if (.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
library(data.table)
options(digits=17)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else '__USER_HOME__/Desktop/TWAS/2026-09-22_Hatton_EUR_mQTL_rebuild'
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,'QC'),showWarnings=FALSE)
smr_root <- '__SMR_ROOT__/mqtl1'
ann_root <- '__USER_HOME__/Desktop/smr/甲基化注释'
short_file <- file.path(ann_root,'注释.csv')
gpl_file <- file.path(ann_root,'GPL13534_HumanMethylation450_15017482_v.1.1.csv')
brain_file <- '__USER_HOME__/Desktop/TWAS/2026-06-18 修稿/数据源附表/smr-new/有意义-eBrain-mqtl-fdr注释.csv'
brain_backup <- '__USER_HOME__/Desktop/TWAS/结果存档/smr-new/有意义-eBrain-mqtl-fdr注释.csv'
old_file <- file.path(ann_root,'注释-eur-fdr_mqtl.csv')
num <- function(x) suppressWarnings(as.numeric(x))
missing <- function(x) is.na(x)|trimws(x)%in%c('','NA','NaN','N/A')
write_tsv <- function(d,n) fwrite(d,file.path(out,n),sep='\t',na='NA',quote='auto')
inputs <- file.path(smr_root,sprintf('msqtl_asd_chr%d.smr',1:22))
stopifnot(all(file.exists(inputs)),all(file.info(inputs)$size>0))
raw <- lapply(seq_along(inputs),function(i){
 z=fread(inputs[i],colClasses='character',na.strings=NULL)
 stopifnot(all(c('probeID','ProbeChr','p_SMR','p_HEIDI','nsnp_HEIDI')%in%names(z)),all(num(z$ProbeChr)==i))
 z[,`:=`(Source_chr=i,Source_row=seq_len(.N),Source_file=inputs[i])];z
})
original_cols <- setdiff(names(raw[[1]]),c('Source_chr','Source_row','Source_file'))
stopifnot(all(vapply(raw,function(z)identical(setdiff(names(z),c('Source_chr','Source_row','Source_file')),original_cols),logical(1))))
a <- rbindlist(raw,use.names=TRUE)
a[,p_SMR_numeric:=num(p_SMR)]
a[,Valid_p_SMR:=is.finite(p_SMR_numeric)&p_SMR_numeric>=0&p_SMR_numeric<=1]
a[,FDR_EUR:=NA_real_]
a[Valid_p_SMR==TRUE,FDR_EUR:=p.adjust(p_SMR_numeric,method='BH')]
# Chromosome-specific FDR is calculated only to audit the OLD result, never for formal selection.
a[,FDR_chr_audit:=NA_real_]
a[Valid_p_SMR==TRUE,FDR_chr_audit:=p.adjust(p_SMR_numeric,'BH'),by=Source_chr]
qc <- a[,.(Rows=.N,Valid_p_SMR=sum(Valid_p_SMR),Invalid_p_SMR=sum(!Valid_p_SMR),Zero_p_SMR=sum(p_SMR_numeric==0,na.rm=TRUE),Formal_significant_rows=sum(FDR_EUR<.05,na.rm=TRUE)),by=Source_chr]
dups <- a[,.(Rows=.N,N_chromosomes=uniqueN(Source_chr),Chromosomes=paste(sort(unique(Source_chr)),collapse=';')),by=probeID][Rows>1]
write_tsv(qc,'QC/chromosome_QC.tsv');write_tsv(dups,'QC/duplicate_probeID.tsv')
write_tsv(a[Valid_p_SMR==FALSE],'QC/invalid_p_SMR.tsv')
write_tsv(a,'01_Hatton_EUR_all_chr_merged_withFDR.tsv')
# Exact probeID joins only; no row expansion or gene-based overlap.
short <- fread(short_file,colClasses='character',na.strings=NULL)
stopifnot(!anyDuplicated(short$probeID))
# This Illumina manifest ends its assay table with a separate [Controls] section.
# fread may warn that this footer has 5 rather than 33 fields; it retains the assay table.
# Annotation matching/coordinates are independently checked in QC outputs.
gpl <- fread(gpl_file,skip='IlmnID,Name,',colClasses='character',na.strings=NULL)
stopifnot(!anyDuplicated(gpl$IlmnID))
brain <- fread(brain_file,colClasses='character',na.strings=NULL)
stopifnot(nrow(brain)==21,uniqueN(brain$probeID)==21,all(num(brain$fdr_dan)<.05))
if(file.exists(brain_backup)){
 b2=fread(brain_backup,colClasses='character',na.strings=NULL)
 stopifnot(setequal(brain$probeID,b2$probeID),identical(brain[order(probeID)],b2[order(probeID)]))
}
old <- fread(old_file,colClasses='character',na.strings=NULL)
stopifnot(!anyDuplicated(old$probeID),'fdr'%in%names(old))
old_sig <- old[is.finite(num(fdr))&num(fdr)<.05]
standardize <- function(annot,smrgene){
 x=c(if(!missing(annot))strsplit(annot,';',fixed=TRUE)[[1]] else character(),if(!missing(smrgene))strsplit(smrgene,';',fixed=TRUE)[[1]] else character())
 x=trimws(x);x=x[!missing(x)];x[x=='KIAA1267']='KANSL1';x=unique(x)
 if(length(x))paste(x,collapse=';') else NA_character_
}
annotate <- function(x){
 z=copy(x);is=match(z$probeID,short$probeID);ig=match(z$probeID,gpl$IlmnID)
 for(n in setdiff(names(short),'probeID'))z[,(paste0('Annotation_',n)):=short[[n]][is]]
 for(n in names(gpl))z[,(paste0('GPL_',n)):=gpl[[n]][ig]]
 z[,Original_annotation_gene:=Annotation_UCSC_RefGene_Name]
 use_gpl=missing(z$Original_annotation_gene)&!missing(z$GPL_UCSC_RefGene_Name)
 z[,Annotation_gene_for_standardization:=Original_annotation_gene]
 z[use_gpl,Annotation_gene_for_standardization:=GPL_UCSC_RefGene_Name]
 z[,Standardized_gene:=mapply(standardize,Annotation_gene_for_standardization,Gene,USE.NAMES=FALSE)]
 z[,`:=`(hg19_CHR=as.character(GPL_CHR),hg19_BP=num(GPL_MAPINFO))]
 good=num(z$GPL_Genome_Build)==37&!missing(z$GPL_CHR)&is.finite(z$hg19_BP)
 good[is.na(good)]=FALSE
 z[!good,`:=`(hg19_CHR=NA_character_,hg19_BP=NA_real_)]
 z[,In_17q21.31_core:=ifelse(is.na(hg19_CHR)|is.na(hg19_BP),NA,hg19_CHR=='17'&hg19_BP>=43492082&hg19_BP<=44913631)]
 z[,SMR_minus_GPL_bp:=num(Probe_bp)-hg19_BP]
 z[,Coordinate_source:='GPL13534 Genome_Build=37 CHR/MAPINFO; original SMR coordinates preserved']
 z[,HEIDI_assessment:=ifelse(missing(p_HEIDI),'Unavailable',ifelse(is.finite(num(p_HEIDI))&num(p_HEIDI)>=0&num(p_HEIDI)<=1,ifelse(num(p_HEIDI)>.05,'No heterogeneity detected','Heterogeneity detected'),'Invalid'))]
 z
}
sig <- annotate(a[!is.na(FDR_EUR)&FDR_EUR<.05]);setorder(sig,Source_chr,p_SMR_numeric,probeID)
stopifnot(!anyDuplicated(sig$probeID),!any(sig$HEIDI_assessment=='Invalid'))
front=c('probeID','Standardized_gene','In_17q21.31_core','FDR_EUR','p_SMR','b_SMR','se_SMR','p_HEIDI','nsnp_HEIDI','HEIDI_assessment','Original_annotation_gene','ProbeChr','Probe_bp','hg19_CHR','hg19_BP')
setcolorder(sig,c(front,setdiff(names(sig),front)))
write_tsv(sig,'02_Hatton_EUR_FDRsig_annotated.tsv')
brain_ann=annotate(brain)
# Include EUR non-significant values for Brain-only probes wherever the EUR test exists.
ids=sort(union(sig$probeID,brain$probeID));ei=match(ids,a$probeID);bi=match(ids,brain$probeID)
overlap=data.table(probeID=ids,EUR_significant=ids%in%sig$probeID,Brain_mMeta_significant=ids%in%brain$probeID)
overlap[,Category:=ifelse(EUR_significant&Brain_mMeta_significant,'Overlap',ifelse(EUR_significant,'EUR-only','Brain-mMeta-only'))]
overlap[,`:=`(EUR_p_SMR=a$p_SMR_numeric[ei],EUR_FDR=a$FDR_EUR[ei],EUR_p_HEIDI=num(a$p_HEIDI[ei]),EUR_nsnp_HEIDI=num(a$nsnp_HEIDI[ei]),Brain_mMeta_p_SMR=num(brain$p_SMR[bi]),Brain_mMeta_FDR=num(brain$fdr_dan[bi]),Brain_mMeta_p_HEIDI=num(brain$p_HEIDI[bi]),Brain_mMeta_nsnp_HEIDI=num(brain$nsnp_HEIDI[bi]))]
for(n in c('Standardized_gene','Original_annotation_gene','In_17q21.31_core','HEIDI_assessment')){
 overlap[,(paste0('EUR_',n)):=sig[[n]][match(ids,sig$probeID)]]
 overlap[,(paste0('Brain_mMeta_',n)):=brain_ann[[n]][bi]]
}
setorder(overlap,Category,probeID);write_tsv(overlap,'04_EUR_vs_BrainMeta_probe_overlap.tsv')
# Old/new comparison includes the union of both significant sets and the new whole-genome FDR for lost probes.
ci=sort(union(old_sig$probeID,sig$probeID));oi=match(ci,old$probeID);ni=match(ci,a$probeID)
cmp=data.table(probeID=ci,Old_significant=ci%in%old_sig$probeID,New_significant=ci%in%sig$probeID)
cmp[,Status:=ifelse(Old_significant&New_significant,'Retained',ifelse(Old_significant,'Lost after genome-wide FDR','Newly retained'))]
cmp[,`:=`(Old_FDR=num(old$fdr[oi]),New_FDR_EUR=a$FDR_EUR[ni],Reconstructed_chromosome_FDR=a$FDR_chr_audit[ni],Old_p_SMR=num(old$p_SMR[oi]),New_p_SMR=a$p_SMR_numeric[ni],ProbeChr=a$ProbeChr[ni],Probe_bp=a$Probe_bp[ni])]
annoall=annotate(a[match(ci,probeID)])
cmp[,`:=`(Original_annotation_gene=annoall$Original_annotation_gene,Standardized_gene=annoall$Standardized_gene,In_17q21.31_core=annoall$In_17q21.31_core)]
setorder(cmp,Status,New_FDR_EUR,probeID);write_tsv(cmp,'06_old_vs_new_EUR_FDR_comparison.tsv')
# Annotation and exact-source checks are emitted separately without deleting or replacing annotation fields.
anno_qc=rbindlist(list(sig[,.(probeID,Dataset='EUR significant',hg19_CHR,hg19_BP,ProbeChr,Probe_bp,SMR_minus_GPL_bp)],brain_ann[,.(probeID,Dataset='Brain-mMeta significant',hg19_CHR,hg19_BP,ProbeChr,Probe_bp,SMR_minus_GPL_bp)]))
write_tsv(anno_qc,'QC/annotation_coordinate_QC.tsv')
maxdiff=function(x,y){d=abs(x-y);if(all(is.na(d)))NA_real_ else max(d,na.rm=TRUE)}
old_fdr_delta=maxdiff(cmp$Old_FDR,cmp$Reconstructed_chromosome_FDR)
old_p_delta=maxdiff(cmp$Old_p_SMR,cmp$New_p_SMR)
related=function(g)sum(vapply(strsplit(ifelse(is.na(sig$Standardized_gene),'',sig$Standardized_gene),';',fixed=TRUE),function(v)g%in%v,logical(1)))
n=nrow(sig);k=sum(sig$In_17q21.31_core,na.rm=TRUE);inter=intersect(sig$probeID,brain$probeID)
summary=c(sprintf('Total SMR tests: %d',nrow(a)),sprintf('Valid p_SMR: %d',sum(a$Valid_p_SMR)),sprintf('FDR-significant probes: %d',n),sprintf('17q21.31 probes: %d (%.6f%%)',k,100*k/n),sprintf('MAPT-related: %d',related('MAPT')),sprintf('KANSL1/KIAA1267-related: %d',related('KANSL1')),sprintf('CRHR1-related: %d',related('CRHR1')),sprintf('HEIDI no heterogeneity: %d',sum(sig$HEIDI_assessment=='No heterogeneity detected')),sprintf('HEIDI heterogeneity: %d',sum(sig$HEIDI_assessment=='Heterogeneity detected')),sprintf('HEIDI unavailable: %d',sum(sig$HEIDI_assessment=='Unavailable')),sprintf('Overlap with Brain-mMeta: %d',length(inter)),sprintf('Brain-mMeta-only: %d',length(setdiff(brain$probeID,sig$probeID))),sprintf('EUR-only: %d',length(setdiff(sig$probeID,brain$probeID))),sprintf('Old EUR significant count: %d',nrow(old_sig)),sprintf('New EUR significant count: %d',n),sprintf('Brain-mMeta significant probes reproduced in EUR: %d/21 (%.6f%%)',length(inter),100*length(inter)/21),sprintf('EUR significant probes also present in Brain-mMeta: %d/%d (%.6f%%)',length(inter),n,100*length(inter)/n),'',sprintf('Duplicate probeID groups: %d; across chromosomes: %d',nrow(dups),sum(dups$N_chromosomes>1)),sprintf('Old/new retained: %d; lost: %d; newly retained: %d',sum(cmp$Status=='Retained'),sum(cmp$Status=='Lost after genome-wide FDR'),sum(cmp$Status=='Newly retained')),sprintf('Old significant probes in core: %d/%d',sum(cmp[Old_significant==TRUE]$In_17q21.31_core,na.rm=TRUE),nrow(old_sig)),sprintf('Old p_SMR versus original EUR maximum absolute difference: %.17g',old_p_delta),sprintf('Old FDR versus reconstructed chromosome-specific BH maximum difference: %.17g',old_fdr_delta),sprintf('New FDR range among significant probes: %.17g to %.17g',min(sig$FDR_EUR),max(sig$FDR_EUR)),'','Exact overlap probes:',paste(sort(inter),collapse='; '),'Brain-mMeta-only probes:',paste(sort(setdiff(brain$probeID,sig$probeID)),collapse='; '),'EUR-only probes:',paste(sort(setdiff(sig$probeID,brain$probeID)),collapse='; '),'','Rows per chromosome:',capture.output(print(qc)),'','METHODS AND QUALIFICATIONS','No SMR or HEIDI analysis was rerun. Existing chr1-22 files were merged with all original columns preserved.','Only finite numeric p_SMR values in [0,1] enter genome-wide BH; all tests enter before any locus/gene/HEIDI filtering.','Formal selection is FDR_EUR < 0.05 with unrounded values. The chromosome-specific audit FDR is never used for formal selection.','No nsnp_HEIDI threshold is used. Missing HEIDI is Unavailable, not no heterogeneity.','All original reduced and GPL annotations are retained in Annotation_* and GPL_* columns. Only Standardized_gene replaces exact KIAA1267 tokens by KANSL1 and deduplicates gene tokens.','Standardized_gene combines annotation gene tokens with original SMR Gene tokens; when reduced annotation is missing, GPL RefGene annotation is used. No new positional gene assignment is inferred.','Core assignment uses GPL Genome_Build=37 CHR/MAPINFO and the inclusive interval chr17:43492082-44913631; original SMR ProbeChr/Probe_bp and coordinate offsets are retained.','Probe overlap means exact CpG identifiers significant in both datasets, not a test of independent replication or effect-direction agreement.','17q21.31 concentration persists descriptively if the reported core count remains high; no enrichment test was requested or performed.',paste('Brain-mMeta reference:',brain_file),paste('Old EUR reference:',old_file))
writeLines(summary,file.path(out,'07_Hatton_EUR_summary.txt'),useBytes=TRUE)
write_tsv(data.table(Source=c(inputs,short_file,gpl_file,brain_file,brain_backup,old_file),MD5=unname(tools::md5sum(c(inputs,short_file,gpl_file,brain_file,brain_backup,old_file)))),'QC/input_checksums.tsv')
capture.output(sessionInfo(),file=file.path(out,'QC/R_sessionInfo.txt'))
cat(paste(summary[1:24],collapse='\n'),'\n')
