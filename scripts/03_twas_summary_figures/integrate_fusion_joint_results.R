# ============================================================
# 整合17q21.31区域FUSION条件/联合分析结果与基因注释。
# 原始记录日期：2026-09-11。
# ============================================================
# 输入为既有扩展区间post-process、included/dropped及边际关联结果。
# 生成会话中的注释表和基因汇总，并输出描述性统计。

root <- "__USER_HOME__/Desktop/TWAS/FUSION/fusion_twas-master"

# 1. 读取29种组织的GTEx v8 Whole Blood及GENCODE v26更新结果。
fusion_final_file <- file.path(
  root,
  "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
)

fusion_all <- read.delim(
  fusion_final_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. 读取既有扩展区间post-process输入文件。
extended_post_file <- file.path(
  root,
  "ASD.FUSION.17q21.31.GRCh38.extended.alltissueFDR005.postprocess.dat"
)

extended_post <- read.table(
  extended_post_file,
  header = TRUE,
  as.is = TRUE
)

# 3. 读取FUSION条件/联合分析的included与dropped结果。
inc_file <- file.path(
  root,
  "POST_GRCh38",
  "ASD.FUSION.17q21.31.extended.joint_included.dat"
)

drop_file <- file.path(
  root,
  "POST_GRCh38",
  "ASD.FUSION.17q21.31.extended.joint_dropped.dat"
)

joint_inc <- read.delim(
  inc_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

joint_drop <- read.delim(
  drop_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 4. 统一保留状态及条件/联合统计量字段，保留原有字段赋值规则。
joint_inc$STATUS <- "INCLUDED"
joint_drop$STATUS <- "DROPPED"

joint_inc$POST.Z <- joint_inc$JOINT.Z
joint_inc$POST.P <- joint_inc$JOINT.P

joint_drop$POST.Z <- joint_drop$COND.Z
joint_drop$POST.P <- joint_drop$COND.P

joint_all <- rbind(
  joint_inc[, c(
    "FILE", "ID",
    "TWAS.Z", "TWAS.P",
    "POST.Z", "POST.P",
    "STATUS"
  )],
  joint_drop[, c(
    "FILE", "ID",
    "TWAS.Z", "TWAS.P",
    "POST.Z", "POST.P",
    "STATUS"
  )]
)

cat("Joint result rows:", nrow(joint_all), "\n")
print(table(joint_all$STATUS))

# ============================================================
# 5. 按FILE和ID连接输入模型的PANEL及坐标信息。
# ============================================================

post_key <- extended_post[, c(
  "FILE",
  "ID",
  "PANEL",
  "CHR",
  "P0",
  "P1"
)]

joint_all2 <- merge(
  joint_all,
  post_key,
  by = c("FILE", "ID"),
  all.x = TRUE
)

cat(
  "Missing PANEL after merge:",
  sum(is.na(joint_all2$PANEL)),
  "\n"
)
# ============================================================
# 6. 按PANEL和ID连接边际关联统计量及基因注释。
# ============================================================

annotation <- fusion_all[, c(
  "PANEL",
  "ID",
  "ENSEMBL_ID",
  "GENE_SYMBOL",
  "GENE_SYMBOL_GENCODE",
  "GENE_NAME",
  "TWAS.Z",
  "TWAS.P",
  "FDR_ALL",
  "FDR_BRAIN",
  "FDR_NONBRAIN"
)]

# 检查注释表中PANEL与ID组合的唯一性。
cat(
  "Duplicated PANEL+ID in annotation:",
  sum(duplicated(paste(annotation$PANEL, annotation$ID))),
  "\n"
)

joint_annot <- merge(
  joint_all2,
  annotation,
  by = c("PANEL", "ID"),
  all.x = TRUE,
  suffixes = c(".POSTFILE", ".MARGINAL")
)

cat("Final merged rows:", nrow(joint_annot), "\n")

cat(
  "Missing gene annotation:",
  sum(
    is.na(joint_annot$GENE_SYMBOL_GENCODE)
  ),
  "\n"
)

print(table(joint_annot$STATUS))


target_genes <- c(
  "KANSL1-AS1",
  "KANSL1",
  "MAPT",
  "MAPT-AS1",
  "CRHR1",
  "LRRC37A",
  "LRRC37A2"
)

target_joint <- joint_annot[
  joint_annot$GENE_SYMBOL_GENCODE %in% target_genes,
  ,
  drop = FALSE
]

target_joint <- target_joint[
  order(
    target_joint$GENE_SYMBOL_GENCODE,
    target_joint$POST.P
  ),
]

target_show <- target_joint[, c(
  "GENE_SYMBOL_GENCODE",
  "PANEL",
  "ID",
  "TWAS.Z.MARGINAL",
  "TWAS.P.MARGINAL",
  "FDR_ALL",
  "FDR_BRAIN",
  "FDR_NONBRAIN",
  "POST.Z",
  "POST.P",
  "STATUS"
)]

print(
  target_show,
  row.names = FALSE
)

# 按关注基因统计进入条件/联合分析的显著基因–组织模型数。
table(
  target_joint$GENE_SYMBOL_GENCODE
)

# 汇总各关注基因的最小条件/联合分析P值。
aggregate(
  POST.P ~ GENE_SYMBOL_GENCODE,
  data = target_joint,
  FUN = min,
  na.rm = TRUE
)
setdiff(
  target_genes,
  unique(target_joint$GENE_SYMBOL_GENCODE)
)

# ============================================================
# 提取与17q21.31核心区间重叠的全部FUSION模型及基因。
# 此处保留所有可检验模型，不按显著性筛选。
# ============================================================

core_start <- 45414716
core_end   <- 46836265

core_allmodels <- fusion_all[
  fusion_all$CHR == 17 &
    !is.na(fusion_all$P0) &
    !is.na(fusion_all$P1) &
    fusion_all$P0 <= core_end &
    fusion_all$P1 >= core_start,
  ,
  drop = FALSE
]

core_genes_all <- sort(unique(
  core_allmodels$GENE_SYMBOL_GENCODE[
    !is.na(core_allmodels$GENE_SYMBOL_GENCODE) &
      core_allmodels$GENE_SYMBOL_GENCODE != ""
  ]
))

cat(
  "Core region FUSION-testable genes:",
  length(core_genes_all),
  "\n"
)

core_genes_all

core_genes_sig <- sort(unique(
  core_allmodels$GENE_SYMBOL_GENCODE[
    !is.na(core_allmodels$FDR_ALL) &
      core_allmodels$FDR_ALL < 0.05 &
      !is.na(core_allmodels$GENE_SYMBOL_GENCODE)
  ]
))

cat(
  "Core FDR-significant genes:",
  length(core_genes_sig),
  "\n"
)

core_genes_sig

setdiff(
  core_genes_all,
  c(
    "KANSL1-AS1",
    "KANSL1",
    "MAPT",
    "MAPT-AS1",
    "CRHR1",
    "LRRC37A",
    "LRRC37A2"
  )
)

# ============================================================
# 按基因汇总核心区间可检验模型、全组织FDR显著模型及最小边际P值。
# ============================================================

core_gene_summary <- do.call(
  rbind,
  lapply(core_genes_all, function(g) {
    
    x <- core_allmodels[
      core_allmodels$GENE_SYMBOL_GENCODE == g,
      ,
      drop = FALSE
    ]
    
    valid <- !is.na(x$TWAS.P)
    
    sig <- valid &
      !is.na(x$FDR_ALL) &
      x$FDR_ALL < 0.05
    
    if (any(valid)) {
      
      best_i <- which.min(x$TWAS.P)
      
      data.frame(
        GENE = g,
        N_MODELS = sum(valid),
        N_FDR_SIG = sum(sig),
        BEST_TISSUE = x$PANEL[best_i],
        BEST_TWAS_Z = x$TWAS.Z[best_i],
        BEST_TWAS_P = x$TWAS.P[best_i],
        BEST_FDR_ALL = x$FDR_ALL[best_i],
        ENTERED_JOINT = any(sig),
        stringsAsFactors = FALSE
      )
      
    } else {
      
      data.frame(
        GENE = g,
        N_MODELS = 0,
        N_FDR_SIG = 0,
        BEST_TISSUE = NA,
        BEST_TWAS_Z = NA,
        BEST_TWAS_P = NA,
        BEST_FDR_ALL = NA,
        ENTERED_JOINT = FALSE,
        stringsAsFactors = FALSE
      )
    }
  })
)

core_gene_summary <- core_gene_summary[
  order(
    !core_gene_summary$ENTERED_JOINT,
    core_gene_summary$BEST_TWAS_P
  ),
]

print(
  core_gene_summary,
  row.names = FALSE
)
table(core_gene_summary$ENTERED_JOINT)

sum(core_gene_summary$N_FDR_SIG)
core_gene_summary[
  !core_gene_summary$ENTERED_JOINT,
]
