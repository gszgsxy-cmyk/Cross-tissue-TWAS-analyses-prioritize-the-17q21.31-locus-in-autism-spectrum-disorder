jti_dir <- "__USER_HOME__/Desktop/TWAS/MR-JTI-master/model_training/predixcan/src"

jti_files <- list.files(
  jti_dir,
  pattern = "^result_ASD_JTI_.*\\.txt$",
  full.names = TRUE
)

length(jti_files)
basename(jti_files)
header_list <- lapply(
  jti_files,
  function(f) {
    names(
      read.delim(
        f,
        header = TRUE,
        nrows = 2,
        check.names = FALSE
      )
    )
  }
)

table(
  sapply(
    header_list,
    paste,
    collapse = "|"
  )
)
x_test <- read.delim(
  jti_files[1],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(x_test)
head(x_test)


jti_list <- lapply(
  jti_files,
  function(f) {
    
    x <- read.delim(
      f,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    tissue <- basename(f)
    
    tissue <- sub(
      "^result_ASD_JTI_",
      "",
      tissue
    )
    
    tissue <- sub(
      "\\.txt$",
      "",
      tissue
    )
    
    x$TISSUE <- tissue
    
    if (grepl("^Brain_", tissue)) {
      x$TISSUE_GROUP <- "Brain"
    } else {
      x$TISSUE_GROUP <- "Non-brain"
    }
    
    x$SOURCE_FILE <- basename(f)
    
    x
  }
)

jti_all <- do.call(
  rbind,
  jti_list
)

rownames(jti_all) <- NULL
cat(
  "Total rows:",
  nrow(jti_all),
  "\n"
)

cat(
  "Number of tissues:",
  length(unique(jti_all$TISSUE)),
  "\n"
)

table(jti_all$TISSUE_GROUP)


jti_tissue_qc <- as.data.frame(
  table(
    jti_all$TISSUE,
    jti_all$TISSUE_GROUP
  )
)

jti_tissue_qc <- jti_tissue_qc[
  jti_tissue_qc$Freq > 0,
]

jti_tissue_qc


jti_all$GENE_ID_ORIGINAL <- as.character(
  jti_all$gene
)


jti_all$ENSEMBL_ID <- as.character(
  jti_all$gene
)

jti_all$GENE_SYMBOL <- as.character(
  jti_all$genename
)
head(
  jti_all[
    ,
    c(
      "gene",
      "genename",
      "ENSEMBL_ID",
      "GENE_SYMBOL"
    )
  ]
)
jti_all$zscore <- suppressWarnings(
  as.numeric(jti_all$zscore)
)

jti_all$effect_size <- suppressWarnings(
  as.numeric(jti_all$effect_size)
)

jti_all$pvalue <- suppressWarnings(
  as.numeric(jti_all$pvalue)
)
summary(jti_all$zscore)
summary(jti_all$effect_size)
summary(jti_all$pvalue)


valid_jti <- (
  is.finite(jti_all$pvalue) &
    jti_all$pvalue >= 0 &
    jti_all$pvalue <= 1
)
cat(
  "Total JTI tests:",
  nrow(jti_all),
  "\n"
)

cat(
  "Valid JTI tests:",
  sum(valid_jti),
  "\n"
)

cat(
  "Invalid/missing JTI tests:",
  sum(!valid_jti),
  "\n"
)

# 4计算FDR ------------------------------------------------------------------


jti_all$FDR_ALL <- NA_real_

jti_all$FDR_ALL[
  valid_jti
] <- p.adjust(
  jti_all$pvalue[
    valid_jti
  ],
  method = "BH"
)
valid_brain_jti <- (
  valid_jti &
    jti_all$TISSUE_GROUP == "Brain"
)
jti_all$FDR_BRAIN <- NA_real_

jti_all$FDR_BRAIN[
  valid_brain_jti
] <- p.adjust(
  jti_all$pvalue[
    valid_brain_jti
  ],
  method = "BH"
)

valid_nonbrain_jti <- (
  valid_jti &
    jti_all$TISSUE_GROUP == "Non-brain"
)
jti_all$FDR_NONBRAIN <- NA_real_

jti_all$FDR_NONBRAIN[
  valid_nonbrain_jti
] <- p.adjust(
  jti_all$pvalue[
    valid_nonbrain_jti
  ],
  method = "BH"
)
jti_sig_all <- jti_all[
  !is.na(jti_all$FDR_ALL) &
    jti_all$FDR_ALL < 0.05,
  ,
  drop = FALSE
]
jti_sig_brain <- jti_all[
  !is.na(jti_all$FDR_BRAIN) &
    jti_all$FDR_BRAIN < 0.05,
  ,
  drop = FALSE
]
jti_sig_nonbrain <- jti_all[
  !is.na(jti_all$FDR_NONBRAIN) &
    jti_all$FDR_NONBRAIN < 0.05,
  ,
  drop = FALSE
]
cat(
  "\n========== JTI UPDATED SUMMARY ==========\n"
)

cat(
  "Total tests:",
  nrow(jti_all),
  "\n"
)

cat(
  "Valid tests:",
  sum(valid_jti),
  "\n"
)

cat(
  "Invalid/missing tests:",
  sum(!valid_jti),
  "\n"
)

cat(
  "\nAll-tissue valid tests:",
  sum(valid_jti),
  "\n"
)

cat(
  "All-tissue FDR<0.05 features:",
  nrow(jti_sig_all),
  "\n"
)

cat(
  "All-tissue unique genes:",
  length(unique(jti_sig_all$ENSEMBL_ID)),
  "\n"
)

cat(
  "\nBrain valid tests:",
  sum(valid_brain_jti),
  "\n"
)

cat(
  "Brain FDR<0.05 features:",
  nrow(jti_sig_brain),
  "\n"
)

cat(
  "Brain unique genes:",
  length(unique(jti_sig_brain$ENSEMBL_ID)),
  "\n"
)

cat(
  "\nNon-brain valid tests:",
  sum(valid_nonbrain_jti),
  "\n"
)

cat(
  "Non-brain FDR<0.05 features:",
  nrow(jti_sig_nonbrain),
  "\n"
)

cat(
  "Non-brain unique genes:",
  length(unique(jti_sig_nonbrain$ENSEMBL_ID)),
  "\n"
)

cat(
  "=========================================\n"
)

# 9 -----------------------------------------------------------------------
jti_brain_genes <- unique(
  jti_sig_brain$ENSEMBL_ID
)

jti_nonbrain_genes <- unique(
  jti_sig_nonbrain$ENSEMBL_ID
)

jti_shared_genes <- intersect(
  jti_brain_genes,
  jti_nonbrain_genes
)
cat(
  "JTI brain/non-brain shared genes:",
  length(jti_shared_genes),
  "\n"
)
jti_shared_table <- unique(
  jti_all[
    jti_all$ENSEMBL_ID %in% jti_shared_genes,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL"
    )
  ]
)
jti_shared_table <- jti_shared_table[
  !duplicated(
    jti_shared_table$ENSEMBL_ID
  ),
  ,
  drop = FALSE
]
jti_shared_table <- jti_shared_table[
  order(jti_shared_table$GENE_SYMBOL),
  ,
  drop = FALSE
]
print(
  jti_shared_table,
  row.names = FALSE
)


jti_symbol_qc <- aggregate(
  GENE_SYMBOL ~ ENSEMBL_ID,
  data = unique(
    jti_all[
      ,
      c(
        "ENSEMBL_ID",
        "GENE_SYMBOL"
      )
    ]
  ),
  FUN = function(x) {
    length(unique(x))
  }
)

table(jti_symbol_qc$GENE_SYMBOL)
jti_symbol_qc[
  jti_symbol_qc$GENE_SYMBOL > 1,
]



# 12 和fusion联合交集 ----------------------------------------------------------
fusion_shared_file <- paste0(
  "__USER_HOME__/Desktop/TWAS/FUSION/",
  "fusion_twas-master/",
  "ASD.FUSION.GENCODEv26.BRAIN_NONBRAIN_shared32_genes.tsv"
)

fusion_shared <- read.delim(
  fusion_shared_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

nrow(fusion_shared)

fusion_shared_genes <- unique(
  fusion_shared$ENSEMBL_ID
)

cross_shared_ids <- intersect(
  jti_shared_genes,
  fusion_shared_genes
)


names(fusion_shared)

head(fusion_shared)

cat(
  "JTI shared genes:",
  length(jti_shared_genes),
  "\n"
)

cat(
  "FUSION shared genes:",
  length(fusion_shared_genes),
  "\n"
)

cat(
  "JTI ∩ FUSION shared genes:",
  length(cross_shared_ids),
  "\n"
)


cross_table <- merge(
  jti_shared_table,
  fusion_shared[
    ,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL_GENCODE",
      "GENE_TYPE_GENCODE"
    )
  ],
  by = "ENSEMBL_ID",
  all = FALSE
)


names(cross_table)[
  names(cross_table) == "GENE_SYMBOL"
] <- "JTI_GENE_SYMBOL"

names(cross_table)[
  names(cross_table) == "GENE_SYMBOL_GENCODE"
] <- "FUSION_GENE_SYMBOL"


cross_table <- cross_table[
  order(cross_table$FUSION_GENE_SYMBOL),
  ,
  drop = FALSE
]
print(
  cross_table,
  row.names = FALSE
)


cross_table$SYMBOL_MATCH <- (
  cross_table$JTI_GENE_SYMBOL ==
    cross_table$FUSION_GENE_SYMBOL
)
table(
  cross_table$SYMBOL_MATCH,
  useNA = "ifany"
)
cross_table[
  !cross_table$SYMBOL_MATCH,
  ,
  drop = FALSE
]


write.table(
  jti_all,
  file.path(
    jti_dir,
    "ASD.JTI.29tissues.updated.all_withFDR.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  jti_sig_all,
  file.path(
    jti_dir,
    "ASD.JTI.29tissues.updated.ALL_FDR005.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  jti_sig_brain,
  file.path(
    jti_dir,
    "ASD.JTI.29tissues.updated.BRAIN_FDR005_exploratory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  jti_sig_nonbrain,
  file.path(
    jti_dir,
    "ASD.JTI.29tissues.updated.NONBRAIN_FDR005_exploratory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  jti_shared_table,
  file.path(
    jti_dir,
    "ASD.JTI.BRAIN_NONBRAIN.shared_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  cross_table,
  file.path(
    jti_dir,
    "ASD.JTI_FUSION.BRAIN_NONBRAIN.shared_crossmethod_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 矫正方向-------------------------------------------------------------------------
display_map <- setNames(
  cross_table$JTI_GENE_SYMBOL,
  cross_table$ENSEMBL_ID
)

jti_cross_long <- jti_all[
  jti_all$ENSEMBL_ID %in% cross_shared_ids,
  c(
    "ENSEMBL_ID",
    "GENE_SYMBOL",
    "TISSUE",
    "TISSUE_GROUP",
    "zscore",
    "effect_size",
    "pvalue",
    "FDR_ALL",
    "FDR_BRAIN",
    "FDR_NONBRAIN"
  ),
  drop = FALSE
]

jti_cross_long$DISPLAY_SYMBOL <-
  display_map[jti_cross_long$ENSEMBL_ID]

jti_cross_long$DIRECTION <- ifelse(
  jti_cross_long$zscore > 0,
  "Positive",
  ifelse(
    jti_cross_long$zscore < 0,
    "Negative",
    "Zero"
  )
)
jti_cross_long$FDR_GROUP <- ifelse(
  jti_cross_long$TISSUE_GROUP == "Brain",
  jti_cross_long$FDR_BRAIN,
  jti_cross_long$FDR_NONBRAIN
)

jti_cross_long$SIG_GROUP <- (
  !is.na(jti_cross_long$FDR_GROUP) &
    jti_cross_long$FDR_GROUP < 0.05
)
table(
  jti_cross_long$DISPLAY_SYMBOL,
  jti_cross_long$TISSUE_GROUP
)


jti_direction_summary <- do.call(
  rbind,
  lapply(
    split(
      jti_cross_long,
      interaction(
        jti_cross_long$ENSEMBL_ID,
        jti_cross_long$TISSUE_GROUP,
        drop = TRUE
      )
    ),
    function(x) {
      
      data.frame(
        ENSEMBL_ID = x$ENSEMBL_ID[1],
        GENE_SYMBOL = x$DISPLAY_SYMBOL[1],
        TISSUE_GROUP = x$TISSUE_GROUP[1],
        
        N_MODELS = nrow(x),
        
        N_POSITIVE = sum(
          x$zscore > 0,
          na.rm = TRUE
        ),
        
        N_NEGATIVE = sum(
          x$zscore < 0,
          na.rm = TRUE
        ),
        
        N_GROUP_FDR005 = sum(
          x$SIG_GROUP,
          na.rm = TRUE
        ),
        
        N_SIG_POSITIVE = sum(
          x$SIG_GROUP &
            x$zscore > 0,
          na.rm = TRUE
        ),
        
        N_SIG_NEGATIVE = sum(
          x$SIG_GROUP &
            x$zscore < 0,
          na.rm = TRUE
        ),
        
        MIN_Z = min(
          x$zscore,
          na.rm = TRUE
        ),
        
        MAX_Z = max(
          x$zscore,
          na.rm = TRUE
        ),
        
        ALL_POSITIVE = all(
          x$zscore > 0,
          na.rm = TRUE
        ),
        
        ALL_NEGATIVE = all(
          x$zscore < 0,
          na.rm = TRUE
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(jti_direction_summary) <- NULL
jti_direction_summary <- jti_direction_summary[
  order(
    jti_direction_summary$GENE_SYMBOL,
    jti_direction_summary$TISSUE_GROUP
  ),
  ,
  drop = FALSE
]
print(
  jti_direction_summary,
  row.names = FALSE
)

jti_cross_long[
  jti_cross_long$ENSEMBL_ID ==
    "ENSG00000214401",
  c(
    "TISSUE",
    "TISSUE_GROUP",
    "zscore",
    "effect_size",
    "pvalue",
    "FDR_ALL",
    "FDR_GROUP",
    "DIRECTION"
  )
]


write.table(
  jti_cross_long,
  file.path(
    jti_dir,
    "ASD.JTI.crossmethod11_genes.29tissues.full_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
print(jti_direction_summary, row.names = FALSE)


# 合并 ----------------------------------------------------------------------

fusion_file <- paste0(
  "__USER_HOME__/Desktop/TWAS/FUSION/",
  "fusion_twas-master/",
  "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
)

fusion_all <- read.delim(
  fusion_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cat("FUSION rows:", nrow(fusion_all), "\n")
cat("Panels:", length(unique(fusion_all$PANEL)), "\n")
names(fusion_all)


fusion_cross_long <- fusion_all[
  fusion_all$ENSEMBL_ID %in% cross_shared_ids,
  ,
  drop = FALSE
]

cat(
  "Cross-method genes:",
  length(unique(fusion_cross_long$ENSEMBL_ID)),
  "\n"
)

table(
  fusion_cross_long$GENE_SYMBOL_GENCODE,
  fusion_cross_long$TISSUE_GROUP
)


display_map <- setNames(
  cross_table$JTI_GENE_SYMBOL,
  cross_table$ENSEMBL_ID
)

fusion_cross_long$DISPLAY_SYMBOL <-
  display_map[fusion_cross_long$ENSEMBL_ID]


fusion_cross_long$TWAS.Z <- suppressWarnings(
  as.numeric(fusion_cross_long$TWAS.Z)
)

fusion_cross_long$TWAS.P <- suppressWarnings(
  as.numeric(fusion_cross_long$TWAS.P)
)
fusion_cross_long$DIRECTION <- ifelse(
  fusion_cross_long$TWAS.Z > 0,
  "Positive",
  ifelse(
    fusion_cross_long$TWAS.Z < 0,
    "Negative",
    "Zero"
  )
)
fusion_cross_long$FDR_GROUP <- ifelse(
  fusion_cross_long$TISSUE_GROUP == "Brain",
  fusion_cross_long$FDR_BRAIN,
  fusion_cross_long$FDR_NONBRAIN
)

fusion_cross_long$SIG_GROUP <- (
  !is.na(fusion_cross_long$FDR_GROUP) &
    fusion_cross_long$FDR_GROUP < 0.05
)
fusion_direction_summary <- do.call(
  rbind,
  lapply(
    split(
      fusion_cross_long,
      interaction(
        fusion_cross_long$ENSEMBL_ID,
        fusion_cross_long$TISSUE_GROUP,
        drop = TRUE
      )
    ),
    function(x) {
      
      z <- x$TWAS.Z
      
      data.frame(
        ENSEMBL_ID = x$ENSEMBL_ID[1],
        GENE_SYMBOL = x$DISPLAY_SYMBOL[1],
        TISSUE_GROUP = x$TISSUE_GROUP[1],
        
        N_MODELS = sum(is.finite(z)),
        
        N_POSITIVE = sum(
          z > 0,
          na.rm = TRUE
        ),
        
        N_NEGATIVE = sum(
          z < 0,
          na.rm = TRUE
        ),
        
        N_GROUP_FDR005 = sum(
          x$SIG_GROUP,
          na.rm = TRUE
        ),
        
        N_SIG_POSITIVE = sum(
          x$SIG_GROUP &
            z > 0,
          na.rm = TRUE
        ),
        
        N_SIG_NEGATIVE = sum(
          x$SIG_GROUP &
            z < 0,
          na.rm = TRUE
        ),
        
        MIN_Z = if (
          any(is.finite(z))
        ) min(z, na.rm = TRUE) else NA,
        
        MAX_Z = if (
          any(is.finite(z))
        ) max(z, na.rm = TRUE) else NA,
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(fusion_direction_summary) <- NULL

fusion_direction_summary <-
  fusion_direction_summary[
    order(
      fusion_direction_summary$GENE_SYMBOL,
      fusion_direction_summary$TISSUE_GROUP
    ),
    ,
    drop = FALSE
  ]
print(
  fusion_direction_summary,
  row.names = FALSE
)
classify_direction <- function(pos, neg) {
  
  if (pos > 0 & neg == 0) {
    return("Positive")
  }
  
  if (pos == 0 & neg > 0) {
    return("Negative")
  }
  
  if (pos > 0 & neg > 0) {
    return("Mixed")
  }
  
  return("No significant feature")
}
jti_direction_summary$SIG_DIRECTION <- mapply(
  classify_direction,
  jti_direction_summary$N_SIG_POSITIVE,
  jti_direction_summary$N_SIG_NEGATIVE
)
fusion_direction_summary$SIG_DIRECTION <- mapply(
  classify_direction,
  fusion_direction_summary$N_SIG_POSITIVE,
  fusion_direction_summary$N_SIG_NEGATIVE
)
jti_dir_wide <- reshape(
  jti_direction_summary[
    ,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL",
      "TISSUE_GROUP",
      "N_GROUP_FDR005",
      "N_SIG_POSITIVE",
      "N_SIG_NEGATIVE",
      "SIG_DIRECTION"
    )
  ],
  idvar = c(
    "ENSEMBL_ID",
    "GENE_SYMBOL"
  ),
  timevar = "TISSUE_GROUP",
  direction = "wide"
)
fusion_dir_wide <- reshape(
  fusion_direction_summary[
    ,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL",
      "TISSUE_GROUP",
      "N_GROUP_FDR005",
      "N_SIG_POSITIVE",
      "N_SIG_NEGATIVE",
      "SIG_DIRECTION"
    )
  ],
  idvar = c(
    "ENSEMBL_ID",
    "GENE_SYMBOL"
  ),
  timevar = "TISSUE_GROUP",
  direction = "wide"
)
names(jti_dir_wide)[
  names(jti_dir_wide) == "GENE_SYMBOL"
] <- "JTI_GENE_SYMBOL"

names(fusion_dir_wide)[
  names(fusion_dir_wide) == "GENE_SYMBOL"
] <- "FUSION_DISPLAY_SYMBOL"


cross_direction_final <- merge(
  jti_dir_wide,
  fusion_dir_wide,
  by = "ENSEMBL_ID",
  all = FALSE
)
print(
  cross_direction_final,
  row.names = FALSE
)
cross_direction_final$CONSISTENT_POSITIVE <- (
  cross_direction_final$SIG_DIRECTION.Brain.x == "Positive" &
    cross_direction_final$`SIG_DIRECTION.Non-brain.x` == "Positive" &
    cross_direction_final$SIG_DIRECTION.Brain.y == "Positive" &
    cross_direction_final$`SIG_DIRECTION.Non-brain.y` == "Positive"
)
names(cross_direction_final)
print(cross_direction_final, row.names = FALSE)


write.table(
  cross_direction_final,
  file.path(
    jti_dir,
    "ASD.JTI_FUSION.crossmethod11.direction_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
consistent_positive7 <- cross_direction_final[
  cross_direction_final$CONSISTENT_POSITIVE,
  ,
  drop = FALSE
]

write.table(
  consistent_positive7,
  file.path(
    jti_dir,
    "ASD.JTI_FUSION.crossmethod.consistent_positive7.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

consistent_positive7[
  ,
  c(
    "ENSEMBL_ID",
    "JTI_GENE_SYMBOL"
  )
]
