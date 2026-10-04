jti_heatmap21 <- jti_all[
  jti_all$ENSEMBL_ID %in% jti_shared_genes,
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
jti_heatmap21$FDR_GROUP <- ifelse(
  jti_heatmap21$TISSUE_GROUP == "Brain",
  jti_heatmap21$FDR_BRAIN,
  jti_heatmap21$FDR_NONBRAIN
)

jti_heatmap21$SIG_GROUP <- (
  !is.na(jti_heatmap21$FDR_GROUP) &
    jti_heatmap21$FDR_GROUP < 0.05
)
cat(
  "JTI heatmap genes:",
  length(unique(jti_heatmap21$ENSEMBL_ID)),
  "\n"
)

cat(
  "JTI tissues:",
  length(unique(jti_heatmap21$TISSUE)),
  "\n"
)

nrow(jti_heatmap21)


brain_tissues <- c(
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra"
)

nonbrain_tissues <- c(
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum",
  "Colon_Sigmoid",
  "Colon_Transverse",
  "Esophagus_Gastroesophageal_Junction",
  "Esophagus_Mucosa",
  "Esophagus_Muscularis",
  "Minor_Salivary_Gland",
  "Pancreas",
  "Prostate",
  "Small_Intestine_Terminal_Ileum",
  "Stomach",
  "Testis",
  "Uterus",
  "Vagina",
  "Whole_Blood"
)

tissue_order <- c(
  brain_tissues,
  nonbrain_tissues
)
jti_grid21 <- expand.grid(
  ENSEMBL_ID = jti_shared_table$ENSEMBL_ID,
  TISSUE = tissue_order,
  stringsAsFactors = FALSE
)

jti_heatmap21_full <- merge(
  jti_grid21,
  jti_heatmap21,
  by = c(
    "ENSEMBL_ID",
    "TISSUE"
  ),
  all.x = TRUE,
  sort = FALSE
)

symbol_map_jti <- setNames(
  jti_shared_table$GENE_SYMBOL,
  jti_shared_table$ENSEMBL_ID
)

jti_heatmap21_full$DISPLAY_SYMBOL <-
  symbol_map_jti[
    jti_heatmap21_full$ENSEMBL_ID
  ]
jti_heatmap21_full$TISSUE_GROUP <- ifelse(
  jti_heatmap21_full$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)
dim(jti_heatmap21_full)

fusion_heatmap32 <- fusion_all[
  fusion_all$ENSEMBL_ID %in% fusion_shared_genes,
  ,
  drop = FALSE
]

sort(
  unique(
    fusion_heatmap32$PANEL
  )
)

# 开始画图 --------------------------------------------------------------------

tissue_label_map <- c(
  
  "Brain_Amygdala" =
    "Amygdala",
  
  "Brain_Anterior_cingulate_cortex_BA24" =
    "ACC (BA24)",
  
  "Brain_Caudate_basal_ganglia" =
    "Caudate",
  
  "Brain_Cerebellar_Hemisphere" =
    "Cerebellar hemisphere",
  
  "Brain_Cerebellum" =
    "Cerebellum",
  
  "Brain_Cortex" =
    "Cortex",
  
  "Brain_Frontal_Cortex_BA9" =
    "Frontal cortex (BA9)",
  
  "Brain_Hippocampus" =
    "Hippocampus",
  
  "Brain_Hypothalamus" =
    "Hypothalamus",
  
  "Brain_Nucleus_accumbens_basal_ganglia" =
    "Nucleus accumbens",
  
  "Brain_Putamen_basal_ganglia" =
    "Putamen",
  
  "Brain_Spinal_cord_cervical_c-1" =
    "Spinal cord",
  
  "Brain_Substantia_nigra" =
    "Substantia nigra",
  
  "Adipose_Subcutaneous" =
    "Adipose subcutaneous",
  
  "Adipose_Visceral_Omentum" =
    "Adipose visceral",
  
  "Colon_Sigmoid" =
    "Colon sigmoid",
  
  "Colon_Transverse" =
    "Colon transverse",
  
  "Esophagus_Gastroesophageal_Junction" =
    "Esophagus GE junction",
  
  "Esophagus_Mucosa" =
    "Esophagus mucosa",
  
  "Esophagus_Muscularis" =
    "Esophagus muscularis",
  
  "Minor_Salivary_Gland" =
    "Minor salivary gland",
  
  "Pancreas" =
    "Pancreas",
  
  "Prostate" =
    "Prostate",
  
  "Small_Intestine_Terminal_Ileum" =
    "Small intestine",
  
  "Stomach" =
    "Stomach",
  
  "Testis" =
    "Testis",
  
  "Uterus" =
    "Uterus",
  
  "Vagina" =
    "Vagina",
  
  "Whole_Blood" =
    "Whole blood"
)
jti_heatmap21_full$TISSUE_LABEL <-
  tissue_label_map[
    jti_heatmap21_full$TISSUE
  ]
sum(is.na(jti_heatmap21_full$TISSUE_LABEL))

tissue_label_order <-
  unname(
    tissue_label_map[
      tissue_order
    ]
  )

jti_heatmap21_full$TISSUE_LABEL <-
  factor(
    jti_heatmap21_full$TISSUE_LABEL,
    levels = tissue_label_order
  )

gene_order <- sort(
  unique(
    jti_heatmap21_full$DISPLAY_SYMBOL
  )
)
jti_heatmap21_full$DISPLAY_SYMBOL <-
  factor(
    jti_heatmap21_full$DISPLAY_SYMBOL,
    levels = rev(gene_order)
  )
SIG_GROUP


jti_heatmap21_full$SIG_LABEL <- ifelse(
  !is.na(jti_heatmap21_full$SIG_GROUP) &
    jti_heatmap21_full$SIG_GROUP,
  "*",
  ""
)
table(jti_heatmap21_full$SIG_LABEL)


jti_heatmap21_full$STAR_COLOR <- ifelse(
  abs(jti_heatmap21_full$zscore) >= 3,
  "white",
  "black"
)


max_abs_z <- max(
  abs(
    c(
      jti_heatmap21_full$zscore,
      suppressWarnings(
        as.numeric(
          fusion_heatmap32$TWAS.Z
        )
      )
    )
  ),
  na.rm = TRUE
)

max_abs_z

z_limit <- ceiling(
  max_abs_z * 2
) / 2

z_limit

# 画图 ----------------------------------------------------------------------


library(ggplot2)
library(scales)

p_jti21 <- ggplot(
  jti_heatmap21_full,
  aes(
    x = TISSUE_LABEL,
    y = DISPLAY_SYMBOL
  )
) +
  
  geom_tile(
    aes(
      fill = zscore
    ),
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = jti_heatmap21_full[
      jti_heatmap21_full$SIG_LABEL == "*",
      ,
      drop = FALSE
    ],
    aes(
      label = SIG_LABEL,
      color = STAR_COLOR
    ),
    size = 3.4,
    fontface = "bold",
    vjust = 0.55
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradient2(
    name = "TWAS Z-score",
    low = "#3980b8",
    mid = "white",
    high = "#ef3b3c",
    midpoint = 0,
    limits = c(
      -z_limit,
      z_limit
    ),
    oob = scales::squish,
    na.value = "grey85"
  ) +
  
  facet_grid(
    cols = vars(TISSUE_GROUP),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "JTI cross-tissue ASD-associated signals",
    subtitle =
      "* Group-specific BH-FDR < 0.05"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "black"
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5
    ),
    
    axis.text.y = element_text(
      size = 10.5
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right"
  )
p_jti21

# 保存 ----------------------------------------------------------------------

ggsave(
  filename = file.path(
    jti_dir,
    "JTI_shared21_heatmap_reviewer_revision.pdf"
  ),
  plot = p_jti21,
  width = 18,
  height = 9,
  units = "in"
)
ggsave(
  filename = file.path(
    jti_dir,
    "JTI_shared21_heatmap_reviewer_revision.tiff"
  ),
  plot = p_jti21,
  width = 18,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

write.table(
  jti_heatmap21_full,
  file.path(
    jti_dir,
    "Supplementary_JTI_shared21_heatmap_full_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 画fusion的 ----------------------------------------------------------------

fusion_heatmap32$TISSUE <- sub(
  "^GTExv8\\.ALL\\.",
  "",
  fusion_heatmap32$PANEL
)

fusion_heatmap32$TISSUE_GROUP <- ifelse(
  fusion_heatmap32$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)
cat(
  "FUSION heatmap genes:",
  length(unique(fusion_heatmap32$ENSEMBL_ID)),
  "\n"
)

cat(
  "FUSION tissues:",
  length(unique(fusion_heatmap32$TISSUE)),
  "\n"
)

setdiff(
  tissue_order,
  unique(fusion_heatmap32$TISSUE)
)

setdiff(
  unique(fusion_heatmap32$TISSUE),
  tissue_order
)

fusion_dup_qc <- table(
  paste(
    fusion_heatmap32$ENSEMBL_ID,
    fusion_heatmap32$TISSUE,
    sep = "|||"
  )
)

cat(
  "Duplicated gene-tissue pairs:",
  sum(fusion_dup_qc > 1),
  "\n"
)

fusion_heatmap32$TWAS.Z <- suppressWarnings(
  as.numeric(fusion_heatmap32$TWAS.Z)
)

fusion_heatmap32$TWAS.P <- suppressWarnings(
  as.numeric(fusion_heatmap32$TWAS.P)
)

fusion_heatmap32$FDR_GROUP <- ifelse(
  fusion_heatmap32$TISSUE_GROUP == "Brain",
  fusion_heatmap32$FDR_BRAIN,
  fusion_heatmap32$FDR_NONBRAIN
)

fusion_heatmap32$SIG_GROUP <- (
  !is.na(fusion_heatmap32$FDR_GROUP) &
    fusion_heatmap32$FDR_GROUP < 0.05
)

fusion_gene_table <- unique(
  fusion_shared[
    ,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL_GENCODE"
    )
  ]
)

fusion_grid32 <- expand.grid(
  ENSEMBL_ID = fusion_gene_table$ENSEMBL_ID,
  TISSUE = tissue_order,
  stringsAsFactors = FALSE
)
fusion_heatmap32_full <- merge(
  fusion_grid32,
  fusion_heatmap32,
  by = c(
    "ENSEMBL_ID",
    "TISSUE"
  ),
  all.x = TRUE,
  sort = FALSE
)
dim(fusion_heatmap32_full)


fusion_symbol_map <- setNames(
  fusion_gene_table$GENE_SYMBOL_GENCODE,
  fusion_gene_table$ENSEMBL_ID
)

fusion_heatmap32_full$DISPLAY_SYMBOL <-
  fusion_symbol_map[
    fusion_heatmap32_full$ENSEMBL_ID
  ]
fusion_heatmap32_full$TISSUE_GROUP <- ifelse(
  fusion_heatmap32_full$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)

fusion_heatmap32_full$TISSUE_GROUP <- factor(
  fusion_heatmap32_full$TISSUE_GROUP,
  levels = c(
    "Brain",
    "Non-brain"
  )
)


fusion_heatmap32_full$TISSUE_LABEL <-
  tissue_label_map[
    fusion_heatmap32_full$TISSUE
  ]
sum(
  is.na(
    fusion_heatmap32_full$TISSUE_LABEL
  )
)
fusion_heatmap32_full$TISSUE_LABEL <- factor(
  fusion_heatmap32_full$TISSUE_LABEL,
  levels = tissue_label_order
)


fusion_gene_order <- sort(
  unique(
    fusion_heatmap32_full$DISPLAY_SYMBOL
  )
)

fusion_heatmap32_full$DISPLAY_SYMBOL <- factor(
  fusion_heatmap32_full$DISPLAY_SYMBOL,
  levels = rev(
    fusion_gene_order
  )
)

fusion_heatmap32_full$SIG_LABEL <- ifelse(
  !is.na(fusion_heatmap32_full$SIG_GROUP) &
    fusion_heatmap32_full$SIG_GROUP,
  "*",
  ""
)

fusion_heatmap32_full$STAR_COLOR <- ifelse(
  abs(fusion_heatmap32_full$TWAS.Z) >= 3,
  "white",
  "black"
)

table(
  fusion_heatmap32_full$SIG_LABEL
)
z_limit

p_fusion32 <- ggplot(
  fusion_heatmap32_full,
  aes(
    x = TISSUE_LABEL,
    y = DISPLAY_SYMBOL
  )
) +
  
  geom_tile(
    aes(
      fill = TWAS.Z
    ),
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = fusion_heatmap32_full[
      fusion_heatmap32_full$SIG_LABEL == "*",
      ,
      drop = FALSE
    ],
    aes(
      label = SIG_LABEL,
      color = STAR_COLOR
    ),
    size = 3.4,
    fontface = "bold",
    vjust = 0.55
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradient2(
    name = "TWAS Z-score",
    low = "#3980b8",
    mid = "white",
    high = "#ef3b3c",
    midpoint = 0,
    limits = c(
      -z_limit,
      z_limit
    ),
    oob = scales::squish,
    na.value = "grey85"
  ) +
  
  facet_grid(
    cols = vars(TISSUE_GROUP),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "FUSION cross-tissue ASD-associated signals",
    subtitle =
      "* Group-specific BH-FDR < 0.05"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      color = "black"
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5
    ),
    
    axis.text.y = element_text(
      size = 10.5
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right"
  )

p_fusion32

fusion_dir <- "__USER_HOME__/Desktop/TWAS/FUSION/fusion_twas-master"
ggsave(
  filename = file.path(
    fusion_dir,
    "FUSION_shared32_heatmap_reviewer_revision.pdf"
  ),
  plot = p_fusion32,
  width = 18,
  height = 11,
  units = "in"
)
ggsave(
  filename = file.path(
    fusion_dir,
    "FUSION_shared32_heatmap_reviewer_revision.tiff"
  ),
  plot = p_fusion32,
  width = 18,
  height = 11,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


write.table(
  fusion_heatmap32_full,
  file.path(
    fusion_dir,
    "Supplementary_FUSION_shared32_heatmap_full_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 11合并 --------------------------------------------------------------------

cross_summary_long <- rbind(
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "JTI",
    TISSUE_GROUP = "Brain",
    N_POSITIVE = cross_direction_final$N_SIG_POSITIVE.Brain.x,
    N_NEGATIVE = cross_direction_final$N_SIG_NEGATIVE.Brain.x,
    DIRECTION = cross_direction_final$SIG_DIRECTION.Brain.x,
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "JTI",
    TISSUE_GROUP = "Non-brain",
    N_POSITIVE = cross_direction_final$`N_SIG_POSITIVE.Non-brain.x`,
    N_NEGATIVE = cross_direction_final$`N_SIG_NEGATIVE.Non-brain.x`,
    DIRECTION = cross_direction_final$`SIG_DIRECTION.Non-brain.x`,
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "FUSION",
    TISSUE_GROUP = "Brain",
    N_POSITIVE = cross_direction_final$N_SIG_POSITIVE.Brain.y,
    N_NEGATIVE = cross_direction_final$N_SIG_NEGATIVE.Brain.y,
    DIRECTION = cross_direction_final$SIG_DIRECTION.Brain.y,
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "FUSION",
    TISSUE_GROUP = "Non-brain",
    N_POSITIVE = cross_direction_final$`N_SIG_POSITIVE.Non-brain.y`,
    N_NEGATIVE = cross_direction_final$`N_SIG_NEGATIVE.Non-brain.y`,
    DIRECTION = cross_direction_final$`SIG_DIRECTION.Non-brain.y`,
    stringsAsFactors = FALSE
  )
)
dim(cross_summary_long)

cross_summary_long$COLUMN <- paste(
  cross_summary_long$METHOD,
  cross_summary_long$TISSUE_GROUP,
  sep = "\n"
)

column_order <- c(
  "JTI\nBrain",
  "JTI\nNon-brain",
  "FUSION\nBrain",
  "FUSION\nNon-brain"
)

cross_summary_long$COLUMN <- factor(
  cross_summary_long$COLUMN,
  levels = column_order
)

cross_summary_long$COUNT_LABEL <- paste0(
  cross_summary_long$N_POSITIVE,
  "+/",
  cross_summary_long$N_NEGATIVE,
  "-"
)

head(cross_summary_long)

cross_gene_order <- c(
  "AC022784.1",
  "FAM215B",
  "KANSL1-AS1",
  "LRRC37A2",
  "MAPT-IT1",
  "SLC35G5",
  "SPPL2C",
  "AC091132.2",
  "AL133325.3",
  "ARHGAP27",
  "CRHR1"
)
cross_summary_long$GENE_SYMBOL <- factor(
  cross_summary_long$GENE_SYMBOL,
  levels = rev(cross_gene_order)
)
table(cross_summary_long$DIRECTION)
library(ggplot2)
p_cross11 <- ggplot(
  cross_summary_long,
  aes(
    x = COLUMN,
    y = GENE_SYMBOL
  )
) +
  
  geom_tile(
    aes(
      fill = DIRECTION
    ),
    color = "white",
    linewidth = 0.8,
    width = 0.95,
    height = 0.95
  ) +
  
  geom_text(
    aes(
      label = COUNT_LABEL
    ),
    size = 4,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    name = "Direction among\nsignificant features",
    values = c(
      "Positive" = "#ef3b3c",
      "Negative" = "#3980b8",
      "Mixed" = "grey70",
      "No significant feature" = "white"
    ),
    breaks = c(
      "Positive",
      "Negative",
      "Mixed"
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Cross-method direction patterns of 11 shared genes",
    subtitle =
      "Numbers indicate significant positive/negative gene–tissue features (group-specific BH-FDR < 0.05)"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 9.5
    ),
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right"
  )

p_cross11
p_cross11 <- p_cross11 +
  
  geom_hline(
    yintercept = c(
      2.5,
      4.5
    ),
    linewidth = 0.7
  )

p_cross11

# 保存 ----------------------------------------------------------------------


ggsave(
  filename = file.path(
    jti_dir,
    "JTI_FUSION_crossmethod11_direction_summary.pdf"
  ),
  plot = p_cross11,
  width = 8,
  height = 7,
  units = "in"
)

ggsave(
  filename = file.path(
    jti_dir,
    "JTI_FUSION_crossmethod11_direction_summary.tiff"
  ),
  plot = p_cross11,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
write.table(
  cross_summary_long,
  file.path(
    jti_dir,
    "ASD.JTI_FUSION.crossmethod11.direction_plot_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 改好看点 --------------------------------------------------------------------
cross_summary_long$COUNT_LABEL <- paste0(
  cross_summary_long$N_POSITIVE,
  "/",
  cross_summary_long$N_NEGATIVE
)

cross_summary_long$METHOD <- factor(
  cross_summary_long$METHOD,
  levels = c("JTI", "FUSION")
)

cross_summary_long$TISSUE_GROUP <- factor(
  cross_summary_long$TISSUE_GROUP,
  levels = c("Brain", "Non-brain")
)

cross_gene_order <- c(
  "AC022784.1",
  "FAM215B",
  "KANSL1-AS1",
  "LRRC37A2",
  "MAPT-IT1",
  "SLC35G5",
  "SPPL2C",
  "AC091132.2",
  "AL133325.3",
  "ARHGAP27",
  "CRHR1"
)

cross_summary_long$GENE_SYMBOL <- factor(
  cross_summary_long$GENE_SYMBOL,
  levels = rev(cross_gene_order)
)
cross_summary_long$TEXT_COLOR <- ifelse(
  cross_summary_long$DIRECTION %in% c(
    "Positive",
    "Negative"
  ),
  "white",
  "black"
)
library(ggplot2)
library(grid)

p_cross11_v2 <- ggplot(
  cross_summary_long,
  aes(
    x = TISSUE_GROUP,
    y = GENE_SYMBOL
  )
) +
  
  geom_tile(
    aes(
      fill = DIRECTION
    ),
    width = 0.94,
    height = 0.88,
    color = "white",
    linewidth = 1
  ) +
  
  geom_text(
    aes(
      label = COUNT_LABEL,
      color = TEXT_COLOR
    ),
    size = 3.8,
    fontface = "bold"
  ) +
  
  scale_color_identity() +
  
  scale_fill_manual(
    name = "Direction",
    values = c(
      "Positive" = "#ef3b3c",
      "Negative" = "#3980b8",
      "Mixed" = "grey72",
      "No significant feature" = "grey95"
    ),
    breaks = c(
      "Positive",
      "Negative",
      "Mixed"
    )
  ) +
  
  facet_grid(
    cols = vars(METHOD),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    
    title =
      "Cross-method directional patterns of 11 shared genes",
    
    caption =
      "Cell labels indicate positive/negative significant gene–tissue features; significance was defined as group-specific BH-FDR < 0.05."
  ) +
  
  guides(
    fill = guide_legend(
      title.position = "left",
      nrow = 1
    )
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      5,
      "mm"
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      color = NA
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      margin = margin(
        6,
        6,
        6,
        6
      )
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      face = "bold",
      color = "black",
      margin = margin(
        t = 5
      )
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      color = "black",
      margin = margin(
        r = 5
      )
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    plot.caption = element_text(
      size = 9,
      color = "grey30",
      hjust = 0,
      margin = margin(
        t = 8
      )
    ),
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

p_cross11_v2

ggsave(
  file.path(
    jti_dir,
    "JTI_FUSION_crossmethod11_direction_summary_v2.pdf"
  ),
  p_cross11_v2,
  width = 7.5,
  height = 6.7,
  units = "in"
)
ggsave(
  file.path(
    jti_dir,
    "JTI_FUSION_crossmethod11_direction_summary_v2.tiff"
  ),
  p_cross11_v2,
  width = 7.5,
  height = 6.7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# 重新读一下数据 -----------------------------------------------------------------

# ============================================================
# 0. Paths
# ============================================================

jti_dir <- paste0(
  "__USER_HOME__/Desktop/TWAS/",
  "MR-JTI-master/model_training/predixcan/src"
)

fusion_dir <- paste0(
  "__USER_HOME__/Desktop/TWAS/FUSION/",
  "fusion_twas-master"
)


# ============================================================
# 1. Check formal files
# ============================================================

# ============================================================
# 0. Paths
# ============================================================

jti_dir <- paste0(
  "__USER_HOME__/Desktop/TWAS/",
  "MR-JTI-master/model_training/predixcan/src"
)

fusion_dir <- paste0(
  "__USER_HOME__/Desktop/TWAS/FUSION/",
  "fusion_twas-master"
)


# ============================================================
# 1. Check formal files
# ============================================================

core_files <- c(
  
  JTI_ALL =
    file.path(
      jti_dir,
      "ASD.JTI.29tissues.updated.all_withFDR.tsv"
    ),
  
  JTI_SHARED21 =
    file.path(
      jti_dir,
      "ASD.JTI.BRAIN_NONBRAIN.shared_genes.tsv"
    ),
  
  CROSS11 =
    file.path(
      jti_dir,
      "ASD.JTI_FUSION.crossmethod11.direction_summary.tsv"
    ),
  
  FUSION_ALL =
    file.path(
      fusion_dir,
      "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
    ),
  
  FUSION_SHARED32 =
    file.path(
      fusion_dir,
      "ASD.FUSION.GENCODEv26.BRAIN_NONBRAIN_shared32_genes.tsv"
    )
)

print(
  data.frame(
    FILE = names(core_files),
    EXISTS = file.exists(core_files),
    PATH = unname(core_files),
    row.names = NULL
  )
)
stopifnot(
  all(
    file.exists(core_files)
  )
)

# ============================================================
# 2. Reload formal datasets
# ============================================================

jti_all <- read.delim(
  core_files["JTI_ALL"],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

jti_shared_table <- read.delim(
  core_files["JTI_SHARED21"],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cross_direction_final <- read.delim(
  core_files["CROSS11"],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fusion_all <- read.delim(
  core_files["FUSION_ALL"],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fusion_shared <- read.delim(
  core_files["FUSION_SHARED32"],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)



jti_shared_genes <- unique(
  jti_shared_table$ENSEMBL_ID
)

fusion_shared_genes <- unique(
  fusion_shared$ENSEMBL_ID
)

cross_shared_ids <- intersect(
  jti_shared_genes,
  fusion_shared_genes
)

cat(
  "JTI shared:",
  length(jti_shared_genes),
  "\n"
)

cat(
  "FUSION shared:",
  length(fusion_shared_genes),
  "\n"
)

cat(
  "Cross-method shared:",
  length(cross_shared_ids),
  "\n"
)
# ============================================================
# 3. Tissue order
# ============================================================

brain_tissues <- c(
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra"
)

nonbrain_tissues <- c(
  "Adipose_Subcutaneous",
  "Adipose_Visceral_Omentum",
  "Colon_Sigmoid",
  "Colon_Transverse",
  "Esophagus_Gastroesophageal_Junction",
  "Esophagus_Mucosa",
  "Esophagus_Muscularis",
  "Minor_Salivary_Gland",
  "Pancreas",
  "Prostate",
  "Small_Intestine_Terminal_Ileum",
  "Stomach",
  "Testis",
  "Uterus",
  "Vagina",
  "Whole_Blood"
)

tissue_order <- c(
  brain_tissues,
  nonbrain_tissues
)

cat(
  "Brain tissues:",
  length(brain_tissues),
  "\n"
)

cat(
  "Non-brain tissues:",
  length(nonbrain_tissues),
  "\n"
)

##########3

tissue_label_map <- c(
  
  "Brain_Amygdala" =
    "Amygdala",
  
  "Brain_Anterior_cingulate_cortex_BA24" =
    "ACC (BA24)",
  
  "Brain_Caudate_basal_ganglia" =
    "Caudate",
  
  "Brain_Cerebellar_Hemisphere" =
    "Cerebellar hemisphere",
  
  "Brain_Cerebellum" =
    "Cerebellum",
  
  "Brain_Cortex" =
    "Cortex",
  
  "Brain_Frontal_Cortex_BA9" =
    "Frontal cortex (BA9)",
  
  "Brain_Hippocampus" =
    "Hippocampus",
  
  "Brain_Hypothalamus" =
    "Hypothalamus",
  
  "Brain_Nucleus_accumbens_basal_ganglia" =
    "Nucleus accumbens",
  
  "Brain_Putamen_basal_ganglia" =
    "Putamen",
  
  "Brain_Spinal_cord_cervical_c-1" =
    "Spinal cord",
  
  "Brain_Substantia_nigra" =
    "Substantia nigra",
  
  "Adipose_Subcutaneous" =
    "Adipose subcutaneous",
  
  "Adipose_Visceral_Omentum" =
    "Adipose visceral",
  
  "Colon_Sigmoid" =
    "Colon sigmoid",
  
  "Colon_Transverse" =
    "Colon transverse",
  
  "Esophagus_Gastroesophageal_Junction" =
    "Esophagus GE junction",
  
  "Esophagus_Mucosa" =
    "Esophagus mucosa",
  
  "Esophagus_Muscularis" =
    "Esophagus muscularis",
  
  "Minor_Salivary_Gland" =
    "Minor salivary gland",
  
  "Pancreas" =
    "Pancreas",
  
  "Prostate" =
    "Prostate",
  
  "Small_Intestine_Terminal_Ileum" =
    "Small intestine",
  
  "Stomach" =
    "Stomach",
  
  "Testis" =
    "Testis",
  
  "Uterus" =
    "Uterus",
  
  "Vagina" =
    "Vagina",
  
  "Whole_Blood" =
    "Whole blood"
)

tissue_label_order <- unname(
  tissue_label_map[
    tissue_order
  ]
)
# ============================================================
# 4. Rebuild JTI 21-gene heatmap dataset
# ============================================================

jti_heatmap21 <- jti_all[
  jti_all$ENSEMBL_ID %in% jti_shared_genes,
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

jti_heatmap21$FDR_GROUP <- ifelse(
  jti_heatmap21$TISSUE_GROUP == "Brain",
  jti_heatmap21$FDR_BRAIN,
  jti_heatmap21$FDR_NONBRAIN
)

jti_heatmap21$SIG_GROUP <- (
  !is.na(jti_heatmap21$FDR_GROUP) &
    jti_heatmap21$FDR_GROUP < 0.05
)

jti_grid21 <- expand.grid(
  ENSEMBL_ID = jti_shared_table$ENSEMBL_ID,
  TISSUE = tissue_order,
  stringsAsFactors = FALSE
)

jti_heatmap21_full <- merge(
  jti_grid21,
  jti_heatmap21,
  by = c(
    "ENSEMBL_ID",
    "TISSUE"
  ),
  all.x = TRUE,
  sort = FALSE
)



symbol_map_jti <- setNames(
  jti_shared_table$GENE_SYMBOL,
  jti_shared_table$ENSEMBL_ID
)

jti_heatmap21_full$DISPLAY_SYMBOL <-
  symbol_map_jti[
    jti_heatmap21_full$ENSEMBL_ID
  ]

jti_heatmap21_full$TISSUE_GROUP <- ifelse(
  jti_heatmap21_full$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)

jti_heatmap21_full$TISSUE_GROUP <- factor(
  jti_heatmap21_full$TISSUE_GROUP,
  levels = c(
    "Brain",
    "Non-brain"
  )
)

jti_heatmap21_full$TISSUE_LABEL <-
  tissue_label_map[
    jti_heatmap21_full$TISSUE
  ]

jti_heatmap21_full$TISSUE_LABEL <- factor(
  jti_heatmap21_full$TISSUE_LABEL,
  levels = tissue_label_order
)
jti_heatmap21_full$SIG_LABEL <- ifelse(
  !is.na(jti_heatmap21_full$SIG_GROUP) &
    jti_heatmap21_full$SIG_GROUP,
  "*",
  ""
)

jti_heatmap21_full$STAR_COLOR <- ifelse(
  abs(jti_heatmap21_full$zscore) >= 3,
  "white",
  "black"
)
gene_order <- sort(
  unique(
    jti_heatmap21_full$DISPLAY_SYMBOL
  )
)

jti_heatmap21_full$DISPLAY_SYMBOL <- factor(
  jti_heatmap21_full$DISPLAY_SYMBOL,
  levels = rev(
    gene_order
  )
)

# ============================================================
# 5. Rebuild FUSION 32-gene heatmap dataset
# ============================================================

fusion_heatmap32 <- fusion_all[
  fusion_all$ENSEMBL_ID %in% fusion_shared_genes,
  ,
  drop = FALSE
]

fusion_heatmap32$TISSUE <- sub(
  "^GTExv8\\.ALL\\.",
  "",
  fusion_heatmap32$PANEL
)

fusion_heatmap32$TISSUE_GROUP <- ifelse(
  fusion_heatmap32$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)

fusion_heatmap32$TWAS.Z <- suppressWarnings(
  as.numeric(
    fusion_heatmap32$TWAS.Z
  )
)

fusion_heatmap32$TWAS.P <- suppressWarnings(
  as.numeric(
    fusion_heatmap32$TWAS.P
  )
)

fusion_heatmap32$FDR_GROUP <- ifelse(
  fusion_heatmap32$TISSUE_GROUP == "Brain",
  fusion_heatmap32$FDR_BRAIN,
  fusion_heatmap32$FDR_NONBRAIN
)

fusion_heatmap32$SIG_GROUP <- (
  !is.na(fusion_heatmap32$FDR_GROUP) &
    fusion_heatmap32$FDR_GROUP < 0.05
)
###
fusion_dup_qc <- table(
  paste(
    fusion_heatmap32$ENSEMBL_ID,
    fusion_heatmap32$TISSUE,
    sep = "|||"
  )
)

cat(
  "Duplicated FUSION gene-tissue pairs:",
  sum(fusion_dup_qc > 1),
  "\n"
)
###
fusion_gene_table <- unique(
  fusion_shared[
    ,
    c(
      "ENSEMBL_ID",
      "GENE_SYMBOL_GENCODE"
    )
  ]
)

fusion_grid32 <- expand.grid(
  ENSEMBL_ID =
    fusion_gene_table$ENSEMBL_ID,
  TISSUE =
    tissue_order,
  stringsAsFactors = FALSE
)

fusion_heatmap32_full <- merge(
  fusion_grid32,
  fusion_heatmap32,
  by = c(
    "ENSEMBL_ID",
    "TISSUE"
  ),
  all.x = TRUE,
  sort = FALSE
)

########
fusion_symbol_map <- setNames(
  fusion_gene_table$GENE_SYMBOL_GENCODE,
  fusion_gene_table$ENSEMBL_ID
)

fusion_heatmap32_full$DISPLAY_SYMBOL <-
  fusion_symbol_map[
    fusion_heatmap32_full$ENSEMBL_ID
  ]

fusion_heatmap32_full$TISSUE_GROUP <- ifelse(
  fusion_heatmap32_full$TISSUE %in% brain_tissues,
  "Brain",
  "Non-brain"
)

fusion_heatmap32_full$TISSUE_GROUP <- factor(
  fusion_heatmap32_full$TISSUE_GROUP,
  levels = c(
    "Brain",
    "Non-brain"
  )
)

fusion_heatmap32_full$TISSUE_LABEL <-
  tissue_label_map[
    fusion_heatmap32_full$TISSUE
  ]

fusion_heatmap32_full$TISSUE_LABEL <- factor(
  fusion_heatmap32_full$TISSUE_LABEL,
  levels = tissue_label_order
)





fusion_heatmap32_full$SIG_LABEL <- ifelse(
  !is.na(fusion_heatmap32_full$SIG_GROUP) &
    fusion_heatmap32_full$SIG_GROUP,
  "*",
  ""
)

fusion_heatmap32_full$STAR_COLOR <- ifelse(
  abs(
    fusion_heatmap32_full$TWAS.Z
  ) >= 3,
  "white",
  "black"
)

fusion_gene_order <- sort(
  unique(
    fusion_heatmap32_full$DISPLAY_SYMBOL
  )
)

fusion_heatmap32_full$DISPLAY_SYMBOL <- factor(
  fusion_heatmap32_full$DISPLAY_SYMBOL,
  levels = rev(
    fusion_gene_order
  )
)

# ============================================================
# 6. Shared Z-score range
# ============================================================

max_abs_z <- max(
  abs(
    c(
      jti_heatmap21_full$zscore,
      fusion_heatmap32_full$TWAS.Z
    )
  ),
  na.rm = TRUE
)

z_limit <- ceiling(
  max_abs_z * 2
) / 2

cat(
  "Shared Z-score limit:",
  -z_limit,
  "to",
  z_limit,
  "\n"
)

# ============================================================
# 7. Figure data QC
# ============================================================

cat(
  "\n===== JTI HEATMAP =====\n"
)

cat(
  "Genes:",
  length(
    unique(
      jti_heatmap21_full$ENSEMBL_ID
    )
  ),
  "\n"
)

cat(
  "Tissues:",
  length(
    unique(
      jti_heatmap21_full$TISSUE
    )
  ),
  "\n"
)

cat(
  "Grid:",
  nrow(jti_heatmap21_full),
  "\n"
)

cat(
  "Available models:",
  sum(
    !is.na(
      jti_heatmap21_full$zscore
    )
  ),
  "\n"
)

cat(
  "Missing models:",
  sum(
    is.na(
      jti_heatmap21_full$zscore
    )
  ),
  "\n"
)

cat(
  "Significant cells:",
  sum(
    jti_heatmap21_full$SIG_GROUP,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\n===== FUSION HEATMAP =====\n"
)

cat(
  "Genes:",
  length(
    unique(
      fusion_heatmap32_full$ENSEMBL_ID
    )
  ),
  "\n"
)

cat(
  "Tissues:",
  length(
    unique(
      fusion_heatmap32_full$TISSUE
    )
  ),
  "\n"
)

cat(
  "Grid:",
  nrow(fusion_heatmap32_full),
  "\n"
)

cat(
  "Available Z:",
  sum(
    !is.na(
      fusion_heatmap32_full$TWAS.Z
    )
  ),
  "\n"
)

cat(
  "Missing Z:",
  sum(
    is.na(
      fusion_heatmap32_full$TWAS.Z
    )
  ),
  "\n"
)

cat(
  "Significant cells:",
  sum(
    fusion_heatmap32_full$SIG_GROUP,
    na.rm = TRUE
  ),
  "\n"
)
save(
  jti_heatmap21_full,
  fusion_heatmap32_full,
  cross_direction_final,
  brain_tissues,
  nonbrain_tissues,
  tissue_order,
  tissue_label_map,
  tissue_label_order,
  z_limit,
  
  file = file.path(
    jti_dir,
    "ASD_TWAS_Figure_BCD_ready.RData"
  )
)

# 继续跑 ---------------------------------------------------------------------


cat("Genes:", length(unique(jti_heatmap21_full$ENSEMBL_ID)), "\n")
cat("Tissues:", length(unique(jti_heatmap21_full$TISSUE)), "\n")
cat("Grid:", nrow(jti_heatmap21_full), "\n")
cat("Available models:", sum(!is.na(jti_heatmap21_full$zscore)), "\n")
cat("Missing models:", sum(is.na(jti_heatmap21_full$zscore)), "\n")
cat("Significant cells:", sum(jti_heatmap21_full$SIG_GROUP, na.rm = TRUE), "\n")

# ============================================================
# Final strict QC before plotting
# ============================================================

expected_cross11 <- c(
  "ENSG00000248538", # AC022784.1
  "ENSG00000236234", # AC091132.2
  "ENSG00000278041", # AL133325.3
  "ENSG00000159314", # ARHGAP27
  "ENSG00000120088", # CRHR1
  "ENSG00000232300", # FAM215B
  "ENSG00000214401", # KANSL1-AS1
  "ENSG00000238083", # LRRC37A2
  "ENSG00000279685", # MAPT-IT1
  "ENSG00000177710", # SLC35G5
  "ENSG00000185294"  # SPPL2C
)

stopifnot(
  
  # JTI
  length(unique(jti_heatmap21_full$ENSEMBL_ID)) == 21,
  length(unique(jti_heatmap21_full$TISSUE)) == 29,
  nrow(jti_heatmap21_full) == 609,
  sum(!is.na(jti_heatmap21_full$zscore)) == 335,
  sum(is.na(jti_heatmap21_full$zscore)) == 274,
  sum(jti_heatmap21_full$SIG_GROUP, na.rm = TRUE) == 188,
  
  # FUSION
  length(unique(fusion_heatmap32_full$ENSEMBL_ID)) == 32,
  length(unique(fusion_heatmap32_full$TISSUE)) == 29,
  nrow(fusion_heatmap32_full) == 928,
  sum(!is.na(fusion_heatmap32_full$TWAS.Z)) == 539,
  sum(is.na(fusion_heatmap32_full$TWAS.Z)) == 389,
  sum(fusion_heatmap32_full$SIG_GROUP, na.rm = TRUE) == 389,
  
  # Tissue structure
  length(brain_tissues) == 13,
  length(nonbrain_tissues) == 16,
  sum(is.na(jti_heatmap21_full$TISSUE_LABEL)) == 0,
  sum(is.na(fusion_heatmap32_full$TISSUE_LABEL)) == 0,
  
  # Cross-method
  length(cross_shared_ids) == 11,
  setequal(
    cross_shared_ids,
    expected_cross11
  )
)

cat(
  "\nALL FIGURE DATA QC CHECKS PASSED.\n"
)



# 核对数据范围 ------------------------------------------------------------------
# JTI实际进入热图的Z值范围
jti_z_range <- range(
  jti_heatmap21_full$zscore,
  na.rm = TRUE
)

# FUSION实际进入热图的Z值范围
fusion_z_range <- range(
  fusion_heatmap32_full$TWAS.Z,
  na.rm = TRUE
)

# 两张图合并后的真实总体范围
shared_z_range <- range(
  c(
    jti_heatmap21_full$zscore,
    fusion_heatmap32_full$TWAS.Z
  ),
  na.rm = TRUE
)

cat(
  "JTI range:",
  jti_z_range[1],
  "to",
  jti_z_range[2],
  "\n"
)

cat(
  "FUSION range:",
  fusion_z_range[1],
  "to",
  fusion_z_range[2],
  "\n"
)

cat(
  "Combined range:",
  shared_z_range[1],
  "to",
  shared_z_range[2],
  "\n"
)

# JTI minimum
jti_heatmap21_full[
  which.min(jti_heatmap21_full$zscore),
  c(
    "ENSEMBL_ID",
    "DISPLAY_SYMBOL",
    "TISSUE",
    "zscore"
  )
]

# JTI maximum
jti_heatmap21_full[
  which.max(jti_heatmap21_full$zscore),
  c(
    "ENSEMBL_ID",
    "DISPLAY_SYMBOL",
    "TISSUE",
    "zscore"
  )
]

# FUSION minimum
fusion_heatmap32_full[
  which.min(fusion_heatmap32_full$TWAS.Z),
  c(
    "ENSEMBL_ID",
    "DISPLAY_SYMBOL",
    "TISSUE",
    "TWAS.Z"
  )
]

# FUSION maximum
fusion_heatmap32_full[
  which.max(fusion_heatmap32_full$TWAS.Z),
  c(
    "ENSEMBL_ID",
    "DISPLAY_SYMBOL",
    "TISSUE",
    "TWAS.Z"
  )
]

# 重新绘图 --------------------------------------------------------------------

library(ggplot2)
library(grid)
library(scales)
z_limits <- c(-6.1, 6.1)
z_breaks <- seq(-6, 6, 2)
p_jti21 <- ggplot(
  jti_heatmap21_full,
  aes(
    x = TISSUE_LABEL,
    y = DISPLAY_SYMBOL
  )
) +
  
  geom_tile(
    aes(fill = zscore),
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = jti_heatmap21_full[
      jti_heatmap21_full$SIG_LABEL == "*",
      ,
      drop = FALSE
    ],
    aes(
      label = SIG_LABEL,
      color = STAR_COLOR
    ),
    size = 3.4,
    fontface = "bold",
    vjust = 0.55
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradient2(
    name = "TWAS Z-score",
    low = "#3980b8",
    mid = "white",
    high = "#ef3b3c",
    midpoint = 0,
    limits = z_limits,
    breaks = z_breaks,
    na.value = "grey85"
  ) +
  
  facet_grid(
    cols = vars(TISSUE_GROUP),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "B  JTI 21-gene cross-tissue heatmap",
    subtitle = "* Group-specific BH-FDR < 0.05"
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      color = "grey40",
      linewidth = 0.5
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      color = "black"
    ),
    
    plot.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 9.5,
      color = "grey25"
    ),
    
    legend.title = element_text(
      size = 10.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9.5
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

p_jti21


ggsave(
  file.path(
    jti_dir,
    "Figure_B_JTI_shared21_heatmap.pdf"
  ),
  p_jti21,
  width = 18,
  height = 9,
  units = "in"
)

ggsave(
  file.path(
    jti_dir,
    "Figure_B_JTI_shared21_heatmap.tiff"
  ),
  p_jti21,
  width = 18,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)



p_fusion32 <- ggplot(
  fusion_heatmap32_full,
  aes(
    x = TISSUE_LABEL,
    y = DISPLAY_SYMBOL
  )
) +
  
  geom_tile(
    aes(fill = TWAS.Z),
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = fusion_heatmap32_full[
      fusion_heatmap32_full$SIG_LABEL == "*",
      ,
      drop = FALSE
    ],
    aes(
      label = SIG_LABEL,
      color = STAR_COLOR
    ),
    size = 3.2,
    fontface = "bold",
    vjust = 0.55
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradient2(
    name = "TWAS Z-score",
    low = "#3980b8",
    mid = "white",
    high = "#ef3b3c",
    midpoint = 0,
    limits = z_limits,
    breaks = z_breaks,
    na.value = "grey85"
  ) +
  
  facet_grid(
    cols = vars(TISSUE_GROUP),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "C  FUSION 32-gene cross-tissue heatmap",
    subtitle = "* Group-specific BH-FDR < 0.05"
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      color = "grey40",
      linewidth = 0.5
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 9.5,
      color = "black"
    ),
    
    plot.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 9.5,
      color = "grey25"
    ),
    
    legend.title = element_text(
      size = 10.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9.5
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

p_fusion32


ggsave(
  file.path(
    fusion_dir,
    "Figure_C_FUSION_shared32_heatmap.pdf"
  ),
  p_fusion32,
  width = 18,
  height = 11,
  units = "in"
)

ggsave(
  file.path(
    fusion_dir,
    "Figure_C_FUSION_shared32_heatmap.tiff"
  ),
  p_fusion32,
  width = 18,
  height = 11,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# D -----------------------------------------------------------------------


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


cross_summary_long <- rbind(
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "JTI",
    TISSUE_GROUP = "Brain",
    N_POSITIVE = cross_direction_final$N_SIG_POSITIVE.Brain.x,
    N_NEGATIVE = cross_direction_final$N_SIG_NEGATIVE.Brain.x,
    DIRECTION = cross_direction_final$SIG_DIRECTION.Brain.x
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "JTI",
    TISSUE_GROUP = "Non-brain",
    N_POSITIVE =
      cross_direction_final$`N_SIG_POSITIVE.Non-brain.x`,
    N_NEGATIVE =
      cross_direction_final$`N_SIG_NEGATIVE.Non-brain.x`,
    DIRECTION =
      cross_direction_final$`SIG_DIRECTION.Non-brain.x`
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "FUSION",
    TISSUE_GROUP = "Brain",
    N_POSITIVE = cross_direction_final$N_SIG_POSITIVE.Brain.y,
    N_NEGATIVE = cross_direction_final$N_SIG_NEGATIVE.Brain.y,
    DIRECTION = cross_direction_final$SIG_DIRECTION.Brain.y
  ),
  
  data.frame(
    ENSEMBL_ID = cross_direction_final$ENSEMBL_ID,
    GENE_SYMBOL = cross_direction_final$JTI_GENE_SYMBOL,
    METHOD = "FUSION",
    TISSUE_GROUP = "Non-brain",
    N_POSITIVE =
      cross_direction_final$`N_SIG_POSITIVE.Non-brain.y`,
    N_NEGATIVE =
      cross_direction_final$`N_SIG_NEGATIVE.Non-brain.y`,
    DIRECTION =
      cross_direction_final$`SIG_DIRECTION.Non-brain.y`
  )
)


cross_summary_long$COUNT_LABEL <- paste0(
  cross_summary_long$N_POSITIVE,
  "/",
  cross_summary_long$N_NEGATIVE
)


cross_gene_order <- c(
  "AC022784.1",
  "FAM215B",
  "KANSL1-AS1",
  "LRRC37A2",
  "MAPT-IT1",
  "SLC35G5",
  "SPPL2C",
  "AC091132.2",
  "AL133325.3",
  "ARHGAP27",
  "CRHR1"
)

cross_summary_long$GENE_SYMBOL <- factor(
  cross_summary_long$GENE_SYMBOL,
  levels = rev(cross_gene_order)
)

cross_summary_long$METHOD <- factor(
  cross_summary_long$METHOD,
  levels = c("JTI", "FUSION")
)

cross_summary_long$TISSUE_GROUP <- factor(
  cross_summary_long$TISSUE_GROUP,
  levels = c(
    "Brain",
    "Non-brain"
  )
)
cross_summary_long$TEXT_COLOR <- ifelse(
  cross_summary_long$DIRECTION %in% c(
    "Positive",
    "Negative"
  ),
  "white",
  "black"
)


p_cross11 <- ggplot(
  cross_summary_long,
  aes(
    x = TISSUE_GROUP,
    y = GENE_SYMBOL
  )
) +
  
  geom_tile(
    aes(fill = DIRECTION),
    width = 0.94,
    height = 0.88,
    color = "white",
    linewidth = 1
  ) +
  
  geom_text(
    aes(
      label = COUNT_LABEL,
      color = TEXT_COLOR
    ),
    size = 3.8,
    fontface = "bold"
  ) +
  
  scale_color_identity() +
  
  scale_fill_manual(
    name = "Direction",
    values = c(
      "Positive" = "#ef3b3c",
      "Negative" = "#3980b8",
      "Mixed" = "grey72",
      "No significant feature" = "grey95"
    ),
    breaks = c(
      "Positive",
      "Negative",
      "Mixed"
    )
  ) +
  
  facet_grid(
    cols = vars(METHOD),
    scales = "free_x",
    space = "free_x"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title =
      "D  Direction summary for 11 cross-method genes",
    caption =
      "Cell labels indicate positive/negative significant gene–tissue features."
  ) +
  
  guides(
    fill = guide_legend(
      title.position = "left",
      nrow = 1
    )
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      5,
      "mm"
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      color = NA
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 10.5,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      color = "black"
    ),
    
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.caption = element_text(
      size = 9,
      color = "grey30",
      hjust = 0
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    )
  )

p_cross11


ggsave(
  file.path(
    jti_dir,
    "Figure_D_crossmethod11_direction_summary.pdf"
  ),
  p_cross11,
  width = 7.5,
  height = 6.7,
  units = "in"
)

ggsave(
  file.path(
    jti_dir,
    "Figure_D_crossmethod11_direction_summary.tiff"
  ),
  p_cross11,
  width = 7.5,
  height = 6.7,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
