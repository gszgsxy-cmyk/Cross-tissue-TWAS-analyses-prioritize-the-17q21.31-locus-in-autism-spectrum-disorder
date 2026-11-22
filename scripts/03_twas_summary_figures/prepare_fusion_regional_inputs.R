# 准备17q21.31区域FUSION条件/联合分析输入。
# 输入包括历史合并表及GTEx v8 Whole Blood、GENCODE v26更新结果。
# 保留区域筛选、权重路径核查及原有输入文件导出规则。
# 本文件为历史交互式处理记录，重复读取及会话清理步骤均保留。

root <- "__USER_HOME__/Desktop/TWAS/FUSION/fusion_twas-master"

# 定位按29种组织和22条常染色体汇总的历史FUSION结果文件。
fusion_file <- list.files(
  root,
  pattern = "^ASD\\.FUSION\\.all29\\.allchr\\.withFDR\\.dat$",
  full.names = TRUE,
  recursive = TRUE
)

length(fusion_file)
fusion_file

stopifnot(length(fusion_file) == 1)

fusion_all <- read.table(
  fusion_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(fusion_all)
colnames(fusion_all)
required_cols <- c(
  "PANEL",
  "FILE",
  "ID",
  "CHR",
  "P0",
  "P1",
  "TWAS.Z",
  "TWAS.P",
  "FDR_ALL"
)

setdiff(required_cols, colnames(fusion_all))


cat(
  "Number of tissues:",
  length(unique(fusion_all$PANEL)),
  "\n"
)

sort(unique(fusion_all$PANEL))
sig_all <- fusion_all[
  !is.na(fusion_all$FDR_ALL) &
    fusion_all$FDR_ALL < 0.05 &
    !is.na(fusion_all$TWAS.P),
]

cat(
  "All-tissue FDR<0.05 associations:",
  nrow(sig_all),
  "\n"
)

core_start <- 45414716
core_end   <- 46836265

core_sig <- sig_all[
  sig_all$CHR == 17 &
    sig_all$P0 <= core_end &
    sig_all$P1 >= core_start,
]


cat("Core significant associations:", nrow(core_sig), "\n")
cat("Core tissues represented:", length(unique(core_sig$PANEL)), "\n")
cat("Core unique Ensembl IDs:", length(unique(core_sig$ID)), "\n")

cat(
  "Duplicated PANEL+FILE:",
  sum(duplicated(
    paste(core_sig$PANEL, core_sig$FILE, sep = "|")
  )),
  "\n"
)
core_sig2 <- core_sig[
  order(core_sig$FDR_ALL, core_sig$TWAS.P),
]

core_sig2[, c(
  "PANEL",
  "ID",
  "CHR",
  "P0",
  "P1",
  "TWAS.Z",
  "TWAS.P",
  "FDR_ALL"
)]


ext_start <- 44414716
ext_end   <- 47836265

extended_sig <- sig_all[
  sig_all$CHR == 17 &
    sig_all$P0 <= ext_end &
    sig_all$P1 >= ext_start,
]

cat("Extended significant associations:", nrow(extended_sig), "\n")
cat("Extended tissues represented:", length(unique(extended_sig$PANEL)), "\n")
cat("Extended unique Ensembl IDs:", length(unique(extended_sig$ID)), "\n")

cat(
  "Duplicated PANEL+FILE:",
  sum(duplicated(
    paste(extended_sig$PANEL, extended_sig$FILE, sep = "|")
  )),
  "\n"
)
extended_sig2 <- extended_sig[
  order(extended_sig$FDR_ALL, extended_sig$TWAS.P),
]

extended_sig2[, c(
  "PANEL",
  "ID",
  "CHR",
  "P0",
  "P1",
  "TWAS.Z",
  "TWAS.P",
  "FDR_ALL"
)]
core_sig2[, c("PANEL","ID","P0","TWAS.Z","TWAS.P","FDR_ALL")]


dim(fusion_all)
setdiff(required_cols, colnames(fusion_all))


# -------------------------------------------------------------------------

# 切换至FUSION项目根目录，按原有FILE字段检查相对路径。
setwd(root)

# 汇总核心区间的显著关联、基因、组织及重复模型。
cat("Core associations:", nrow(core_sig), "\n")
cat("Core tissues:", length(unique(core_sig$PANEL)), "\n")
cat("Core genes:", length(unique(core_sig$ID)), "\n")
cat(
  "Core duplicated PANEL+FILE:",
  sum(duplicated(paste(core_sig$PANEL, core_sig$FILE, sep="|"))),
  "\n"
)

# 汇总扩展区间的显著关联、基因和组织。
cat("Extended associations:", nrow(extended_sig), "\n")
cat("Extended tissues:", length(unique(extended_sig$PANEL)), "\n")
cat("Extended genes:", length(unique(extended_sig$ID)), "\n")

# 核查FUSION post-process输入所需字段。
needed_post <- c(
  "PANEL","FILE","ID","CHR","P0","P1",
  "MODEL","TWAS.Z","TWAS.P"
)

cat(
  "Missing post-process columns:",
  paste(setdiff(needed_post, colnames(fusion_all)), collapse=", "),
  "\n"
)

# 检查核心区间与扩展区间的权重文件是否存在。
cat(
  "Core weight files found:",
  sum(file.exists(core_sig$FILE)),
  "/",
  nrow(core_sig),
  "\n"
)

cat(
  "Extended weight files found:",
  sum(file.exists(extended_sig$FILE)),
  "/",
  nrow(extended_sig),
  "\n"
)

# 查看核心区间模型的权重文件路径示例。
head(core_sig$FILE, 5)

missing_core <- core_sig[
  !file.exists(core_sig$FILE),
  c("PANEL","FILE","ID","CHR","P0","P1","TWAS.Z","TWAS.P","FDR_ALL")
]

missing_extended <- extended_sig[
  !file.exists(extended_sig$FILE),
  c("PANEL","FILE","ID","CHR","P0","P1","TWAS.Z","TWAS.P","FDR_ALL")
]

missing_core
missing_extended

cat("Missing core:", nrow(missing_core), "\n")
cat("Missing extended:", nrow(missing_extended), "\n")

weights_root <- file.path(root, "WEIGHTS")

for (f in unique(missing_extended$FILE)) {
  
  target <- basename(f)
  
  hit <- list.files(
    weights_root,
    pattern = paste0("^", gsub("\\.", "\\\\.", target), "$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  cat("\nMissing path:\n", f, "\n")
  cat("Found elsewhere:\n")
  print(hit)
}


# 导出权重路径修复前的扩展区间中间表。
write.table(
  extended_sig,
  file.path(root, "TEMP_17q21.31_extended_before_pathfix.dat"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------------------------------


# 检查历史合并结果是否仍含GTEx.Whole_Blood旧版本记录。
old_wb <- fusion_all[
  fusion_all$PANEL == "GTEx.Whole_Blood" |
    grepl("GTEx\\.Whole_Blood", fusion_all$FILE),
]

cat("Old Whole Blood rows in fusion_all:", nrow(old_wb), "\n")

table(old_wb$CHR)

old_wb_17 <- old_wb[old_wb$CHR == 17, ]

cat("Old Whole Blood rows on chr17:", nrow(old_wb_17), "\n")

old_wb_region <- old_wb_17[
  old_wb_17$P0 <= 47836265 &
    old_wb_17$P1 >= 44414716,
]

cat(
  "Old Whole Blood rows in extended 17q21.31:",
  nrow(old_wb_region),
  "\n"
)

old_wb_region[, c(
  "PANEL","FILE","ID","P0","P1",
  "TWAS.Z","TWAS.P","FDR_ALL"
)]

new_wb_region <- fusion_all[
  fusion_all$PANEL == "GTExv8.ALL.Whole_Blood" &
    fusion_all$CHR == 17 &
    fusion_all$P0 <= 47836265 &
    fusion_all$P1 >= 44414716,
]

cat(
  "New GTExv8 Whole Blood rows in extended region:",
  nrow(new_wb_region),
  "\n"
)

new_wb_region[, c(
  "PANEL","FILE","ID","P0","P1",
  "TWAS.Z","TWAS.P","FDR_ALL"
)]



# 核查Whole Blood版本标识并保留原有旧版本排除步骤。

# 列出PANEL字段中包含Whole_Blood的组织面板标识。
sort(unique(
  fusion_all$PANEL[
    grepl("Whole_Blood", fusion_all$PANEL)
  ]
))

# 统计GTEx.Whole_Blood旧版本记录数。
sum(fusion_all$PANEL == "GTEx.Whole_Blood")

# 统计GTExv8.ALL.Whole_Blood面板记录数。
sum(fusion_all$PANEL == "GTExv8.ALL.Whole_Blood")

fusion_no_oldWB <- fusion_all[
  fusion_all$PANEL != "GTEx.Whole_Blood" &
    !grepl("GTEx\\.Whole_Blood", fusion_all$FILE),
]

cat("Rows before:", nrow(fusion_all), "\n")
cat("Rows after removing old WB:", nrow(fusion_no_oldWB), "\n")
cat(
  "Old WB remaining:",
  sum(
    fusion_no_oldWB$PANEL == "GTEx.Whole_Blood" |
      grepl("GTEx\\.Whole_Blood", fusion_no_oldWB$FILE)
  ),
  "\n"
)

cat(
  "Tissues after removing old WB:",
  length(unique(fusion_no_oldWB$PANEL)),
  "\n"
)
wb_new_files <- list.files(
  root,
  pattern = "GTExv8\\.ALL\\.Whole_Blood.*\\.dat$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

length(wb_new_files)
basename(wb_new_files)



# 读取GTEx v8 Whole Blood及GENCODE v26更新后的合并结果。


fusion_final_file <- file.path(
  root,
  "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
)

file.exists(fusion_final_file)
fusion_all <- read.table(
  fusion_final_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(fusion_all)

sort(unique(
  fusion_all$PANEL[
    grepl("Whole_Blood", fusion_all$PANEL)
  ]
))

cat(
  "Old GTEx Whole Blood:",
  sum(fusion_all$PANEL == "GTEx.Whole_Blood"),
  "\n"
)

cat(
  "New GTExv8 Whole Blood:",
  sum(fusion_all$PANEL == "GTExv8.ALL.Whole_Blood"),
  "\n"
)

cat(
  "Number of tissues:",
  length(unique(fusion_all$PANEL)),
  "\n"
)


# 按历史处理顺序清理会话对象并重新读取结果。

rm(
  fusion_all,
  sig_all,
  core_sig,
  extended_sig,
  core_sig2,
  extended_sig2,
  old_wb,
  old_wb_17,
  old_wb_region,
  new_wb_region
)

gc()
rm(list = intersect(
  c(
    "fusion_all",
    "sig_all",
    "core_sig",
    "extended_sig",
    "core_sig2",
    "extended_sig2",
    "old_wb",
    "old_wb_17",
    "old_wb_region",
    "new_wb_region"
  ),
  ls()
))

gc()
root <- "__USER_HOME__/Desktop/TWAS/FUSION/fusion_twas-master"

fusion_final_file <- file.path(
  root,
  "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
)


file.exists(fusion_final_file)
fusion_all <- read.delim(
  fusion_final_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dim(fusion_all)

colnames(fusion_all)

cat("Total rows:", nrow(fusion_all), "\n")
cat("Total columns:", ncol(fusion_all), "\n")

cat(
  "Number of panels:",
  length(unique(fusion_all$PANEL)),
  "\n"
)

sort(unique(
  fusion_all$PANEL[
    grepl("Whole_Blood", fusion_all$PANEL)
  ]
))
grep(
  "FDR",
  colnames(fusion_all),
  value = TRUE
)
colnames(fusion_all)

# 清理区域中间对象，准备最终输入核查。
rm(list = intersect(
  c(
    "sig_all",
    "core_sig",
    "extended_sig",
    "core_sig2",
    "extended_sig2",
    "missing_core",
    "missing_extended"
  ),
  ls()
))

setwd(root)

# ============================================================
# QC-1：核查最终合并结果的行列数、29种组织及Whole Blood版本。
# ============================================================

cat("Rows:", nrow(fusion_all), "\n")
cat("Columns:", ncol(fusion_all), "\n")
cat("Panels:", length(unique(fusion_all$PANEL)), "\n")

cat(
  "Old GTEx.Whole_Blood:",
  sum(fusion_all$PANEL == "GTEx.Whole_Blood"),
  "\n"
)

cat(
  "New GTExv8.ALL.Whole_Blood:",
  sum(fusion_all$PANEL == "GTExv8.ALL.Whole_Blood"),
  "\n"
)

cat(
  "FDR columns:",
  paste(grep("FDR", colnames(fusion_all), value = TRUE), collapse = ", "),
  "\n"
)

stopifnot(
  nrow(fusion_all) > 0,
  length(unique(fusion_all$PANEL)) == 29,
  sum(fusion_all$PANEL == "GTEx.Whole_Blood") == 0,
  sum(fusion_all$PANEL == "GTExv8.ALL.Whole_Blood") > 0,
  "FDR_ALL" %in% colnames(fusion_all)
)

valid_twas <- (
  is.finite(fusion_all$TWAS.P) &
    fusion_all$TWAS.P >= 0 &
    fusion_all$TWAS.P <= 1
)

cat(
  "Valid TWAS tests:",
  sum(valid_twas),
  "\n"
)

cat(
  "FDR_ALL < 0.05:",
  sum(
    valid_twas &
      is.finite(fusion_all$FDR_ALL) &
      fusion_all$FDR_ALL < 0.05
  ),
  "\n"
)
# ============================================================
# QC-2：按既有全组织BH-FDR阈值提取显著FUSION关联。
# ============================================================

sig_all <- fusion_all[
  valid_twas &
    is.finite(fusion_all$FDR_ALL) &
    fusion_all$FDR_ALL < 0.05,
  ,
  drop = FALSE
]

cat("All-tissue FDR-significant associations:", nrow(sig_all), "\n")

# ============================================================
# QC-3：按预定义GRCh38核心区间及两侧各扩展1 Mb的区间筛选模型。
# ============================================================

core_start <- 45414716
core_end   <- 46836265

ext_start  <- 44414716
ext_end    <- 47836265

core_sig <- sig_all[
  sig_all$CHR == 17 &
    is.finite(sig_all$P0) &
    is.finite(sig_all$P1) &
    sig_all$P0 <= core_end &
    sig_all$P1 >= core_start,
  ,
  drop = FALSE
]

extended_sig <- sig_all[
  sig_all$CHR == 17 &
    is.finite(sig_all$P0) &
    is.finite(sig_all$P1) &
    sig_all$P0 <= ext_end &
    sig_all$P1 >= ext_start,
  ,
  drop = FALSE
]

cat(
  "\nCORE:",
  nrow(core_sig), "associations /",
  length(unique(core_sig$ID)), "genes /",
  length(unique(core_sig$PANEL)), "tissues\n"
)

cat(
  "EXTENDED:",
  nrow(extended_sig), "associations /",
  length(unique(extended_sig$ID)), "genes /",
  length(unique(extended_sig$PANEL)), "tissues\n"
)

cat(
  "Old WB in core:",
  sum(core_sig$PANEL == "GTEx.Whole_Blood"),
  "\n"
)

cat(
  "Old WB in extended:",
  sum(extended_sig$PANEL == "GTEx.Whole_Blood"),
  "\n"
)

cat(
  "New GTExv8 WB in core:",
  sum(core_sig$PANEL == "GTExv8.ALL.Whole_Blood"),
  "\n"
)

cat(
  "New GTExv8 WB in extended:",
  sum(extended_sig$PANEL == "GTExv8.ALL.Whole_Blood"),
  "\n"
)

cat(
  "Core duplicated PANEL+FILE:",
  sum(duplicated(paste(core_sig$PANEL, core_sig$FILE))),
  "\n"
)

cat(
  "Extended duplicated PANEL+FILE:",
  sum(duplicated(paste(extended_sig$PANEL, extended_sig$FILE))),
  "\n"
)
# ============================================================
# QC-4：解析模型权重路径并保留对应组织目录内的检索规则。
# ============================================================

make_absolute_file <- function(f) {
  
  f2 <- gsub("\\\\", "/", f)
  
  # 保留已包含Windows盘符的绝对路径。
  if (grepl("^[A-Za-z]:/", f2)) {
    return(f2)
  }
  
  # 去除相对路径开头的 ./。
  f2 <- sub("^\\./", "", f2)
  
  file.path(root, f2)
}


resolve_weight_paths <- function(dat) {
  
  out <- dat
  
  # 将其余相对路径连接至FUSION项目根目录。
  candidate <- vapply(
    out$FILE,
    make_absolute_file,
    FUN.VALUE = character(1)
  )
  
  for (i in seq_len(nrow(out))) {
    
    # 优先使用已存在的候选权重路径。
    if (file.exists(candidate[i])) {
      out$FILE[i] <- normalizePath(
        candidate[i],
        winslash = "/",
        mustWork = TRUE
      )
      next
    }
    
    # 候选路径不存在时，仅在对应组织的权重目录内检索同名文件。
    panel_dir <- file.path(
      root,
      "WEIGHTS",
      out$PANEL[i]
    )
    
    target <- basename(
      gsub("\\\\", "/", out$FILE[i])
    )
    
    hits <- character(0)
    
    if (dir.exists(panel_dir)) {
      
      all_files <- list.files(
        panel_dir,
        recursive = TRUE,
        full.names = TRUE
      )
      
      hits <- all_files[
        basename(all_files) == target
      ]
    }
    
    if (length(hits) == 1) {
      
      out$FILE[i] <- normalizePath(
        hits,
        winslash = "/",
        mustWork = TRUE
      )
      
    } else {
      
      out$FILE[i] <- candidate[i]
    }
  }
  
  out
}


core_input <- resolve_weight_paths(core_sig)

extended_input <- resolve_weight_paths(extended_sig)


cat(
  "Core weight files found:",
  sum(file.exists(core_input$FILE)),
  "/",
  nrow(core_input),
  "\n"
)

cat(
  "Extended weight files found:",
  sum(file.exists(extended_input$FILE)),
  "/",
  nrow(extended_input),
  "\n"
)

missing_core <- core_input[
  !file.exists(core_input$FILE),
  c("PANEL","ID","FILE","TWAS.Z","TWAS.P","FDR_ALL"),
  drop = FALSE
]

missing_extended <- extended_input[
  !file.exists(extended_input$FILE),
  c("PANEL","ID","FILE","TWAS.Z","TWAS.P","FDR_ALL"),
  drop = FALSE
]

missing_core
missing_extended

# 按原有文件名导出核心与扩展区间输入，并重新读取核查。
# -------------------------------------------------------------------------

core_input_file <- file.path(
  root,
  "ASD.FUSION.17q21.31.GRCh38.core.alltissueFDR005.input.dat"
)

extended_input_file <- file.path(
  root,
  "ASD.FUSION.17q21.31.GRCh38.extended.alltissueFDR005.input.dat"
)

write.table(
  core_input,
  core_input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  extended_input,
  extended_input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Core input:", core_input_file, "\n")
cat("Extended input:", extended_input_file, "\n")

core_check <- read.delim(
  core_input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

extended_check <- read.delim(
  extended_input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Core saved rows:", nrow(core_check), "\n")
cat("Extended saved rows:", nrow(extended_check), "\n")

cat(
  "Core files existing:",
  sum(file.exists(core_check$FILE)),
  "/",
  nrow(core_check),
  "\n"
)

cat(
  "Extended files existing:",
  sum(file.exists(extended_check$FILE)),
  "/",
  nrow(extended_check),
  "\n"
)

cat(
  "Old WB in saved core:",
  sum(core_check$PANEL == "GTEx.Whole_Blood"),
  "\n"
)

cat(
  "Old WB in saved extended:",
  sum(extended_check$PANEL == "GTEx.Whole_Blood"),
  "\n"
)
write.table(
  core_input,
  file.path(
    root,
    "ASD.FUSION.17q21.31.GRCh38.core.alltissueFDR005.input.dat"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
# 按FUSION post-process字段顺序构建核心区间输入表。
fusion_cols <- c(
  "PANEL",
  "FILE",
  "ID",
  "CHR",
  "P0",
  "P1",
  "HSQ",
  "BEST.GWAS.ID",
  "BEST.GWAS.Z",
  "EQTL.ID",
  "EQTL.R2",
  "EQTL.Z",
  "EQTL.GWAS.Z",
  "NSNP",
  "NWGT",
  "MODEL",
  "MODELCV.R2",
  "MODELCV.PV",
  "TWAS.Z",
  "TWAS.P"
)

setdiff(fusion_cols, colnames(core_input))
core_post <- core_input[, fusion_cols]

core_post_file <- file.path(
  root,
  "ASD.FUSION.17q21.31.GRCh38.core.alltissueFDR005.postprocess.dat"
)

write.table(
  core_post,
  core_post_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

dim(core_post)

sum(file.exists(core_post$FILE))

head(core_post[, c("PANEL","ID","CHR","P0","P1","TWAS.Z","TWAS.P")])

test_post <- read.table(
  core_post_file,
  as.is = TRUE,
  header = TRUE
)

dim(test_post)
colnames(test_post)

# 按带版本及不带版本的Ensembl ID查看MAPT-AS1关联记录。


fusion_all[
  fusion_all$ID == "ENSG00000264589.2",
  c(
    "PANEL",
    "ID",
    "GENE_SYMBOL",
    "GENE_SYMBOL_GENCODE",
    "GENE_NAME",
    "TWAS.Z",
    "TWAS.P",
    "FDR_ALL"
  )
]
fusion_all[
  fusion_all$ENSEMBL_ID == "ENSG00000264589",
  c(
    "PANEL",
    "ID",
    "GENE_SYMBOL",
    "GENE_SYMBOL_GENCODE",
    "GENE_NAME",
    "TWAS.Z",
    "TWAS.P",
    "FDR_ALL"
  )
]


# 使用相同字段顺序构建扩展区间post-process输入表。
extended_post <- extended_input[, fusion_cols]

extended_post_file <- file.path(
  root,
  "ASD.FUSION.17q21.31.GRCh38.extended.alltissueFDR005.postprocess.dat"
)

write.table(
  extended_post,
  extended_post_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# 使用FUSION post-process对应的read.table参数复核扩展区间输入。
extended_test <- read.table(
  extended_post_file,
  as.is = TRUE,
  header = TRUE
)

dim(extended_test)

sum(file.exists(extended_test$FILE))
