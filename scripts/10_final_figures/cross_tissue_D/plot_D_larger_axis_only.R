options(stringsAsFactors=FALSE)
if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE)
script_file <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)),mustWork=TRUE)
module_dir <- dirname(script_file)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>1) stop("Usage: Rscript plot_D_larger_axis_only.R [output_directory]")
out <- if(length(args)) args[1] else file.path(module_dir,"output")
data_dir <- file.path(module_dir,"data")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
jti_heatmap21_full <- read.delim(file.path(data_dir,"Figure_2D_JTI_shared21_alltissueFDR_plotdata.tsv"),check.names=FALSE)
fusion_heatmap32_full <- read.delim(file.path(data_dir,"Figure_3D_FUSION_shared32_alltissueFDR_plotdata.tsv"),check.names=FALSE)
# Reuse archived plot blocks. No FDR recalculation.
library(ggplot2)
library(grid)
library(scales)
brain <- c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")
nonbrain <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Minor_Salivary_Gland","Pancreas","Prostate","Small_Intestine_Terminal_Ileum","Stomach","Testis","Uterus","Vagina","Whole_Blood")
tissue_order <- c(brain,nonbrain)
prepare <- function(d,zcol) {
  d$SIG_ALL <- !is.na(d$FDR_ALL) & d$FDR_ALL < 0.05
  d$SIG_LABEL <- ifelse(d$SIG_ALL,"*","")
  d$STAR_COLOR <- ifelse(abs(d[[zcol]]) >= 3,"white","black")
  d$TISSUE_GROUP <- factor(d$TISSUE_GROUP,levels=c("Brain","Non-brain"))
  labels <- unique(d[,c("TISSUE","TISSUE_LABEL")])
  labmap <- setNames(as.character(labels$TISSUE_LABEL),labels$TISSUE)
  d$TISSUE_LABEL <- factor(d$TISSUE_LABEL,levels=unname(labmap[tissue_order]))
  d$DISPLAY_SYMBOL <- factor(d$DISPLAY_SYMBOL,levels=rev(sort(unique(d$DISPLAY_SYMBOL))))
  d
}
jti_heatmap21_full <- prepare(jti_heatmap21_full,"zscore")
fusion_heatmap32_full <- prepare(fusion_heatmap32_full,"TWAS.Z")
stopifnot(nrow(jti_heatmap21_full)==609,nrow(fusion_heatmap32_full)==928)
z_limits <- c(-6.1,6.1)
z_breaks <- seq(-6,6,2)
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
    title = "D  JTI 21-gene cross-tissue heatmap",
    subtitle = "* All-tissue BH-FDR < 0.05"
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
    title = "D  FUSION 32-gene cross-tissue heatmap",
    subtitle = "* All-tissue BH-FDR < 0.05"
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

# Only the two axis label sizes change; all other ggplot parameters are retained.
qc <- list()
for (method in c("JTI","FUSION")) {
 original <- if(method=="JTI") p_jti21 else p_fusion32
 revised <- original + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=13))
 a<-ggplot_build(original);b<-ggplot_build(revised)
 stopifnot(identical(a$data,b$data),identical(original$data,revised$data),
           identical(original$labels,revised$labels),
           identical(a$plot$scales$get_scales("fill")$map(c(-6.1,-3,0,3,6.1,NA)),b$plot$scales$get_scales("fill")$map(c(-6.1,-3,0,3,6.1,NA))),
           identical(original$facet,revised$facet))
 oldtheme<-original$theme;newtheme<-revised$theme
 newtheme$axis.text.x$size<-oldtheme$axis.text.x$size
 newtheme$axis.text.y$size<-oldtheme$axis.text.y$size
 stopifnot(identical(oldtheme,newtheme))
 name<-if(method=="JTI") "Figure_2D_JTI_LargerAxis_v1" else "Figure_3D_FUSION_LargerAxis_v1"
 ht<-if(method=="JTI")9 else 11
 for(ext in c("pdf","tiff","png")) {
  target<-file.path(out,paste0(name,if(ext=="png")"_preview" else "",".",ext))
  stopifnot(!file.exists(target))
  if(ext=="pdf") ggsave(target,revised,width=18,height=ht,units="in",device=pdf,useDingbats=FALSE)
  if(ext=="tiff") ggsave(target,revised,width=18,height=ht,units="in",dpi=600,compression="lzw")
  if(ext=="png") ggsave(target,revised,width=18,height=ht,units="in",dpi=120)
 }
 qc[[method]]<-data.frame(method=method,cells=nrow(revised$data),stars=sum(revised$data$SIG_ALL),
  missing_Z=sum(is.na(revised$data[[if(method=="JTI")"zscore" else "TWAS.Z"]])),
  original_x_pt=original$theme$axis.text.x$size,new_x_pt=12,
  original_y_pt=original$theme$axis.text.y$size,new_y_pt=13,
  geom_data_identical=identical(a$data,b$data),other_theme_identical=identical(oldtheme,newtheme))
 cat(method,"completed; only axis font sizes changed.\n")
}
write.table(do.call(rbind,qc),file.path(out,"Axis_Font_Only_QC.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
