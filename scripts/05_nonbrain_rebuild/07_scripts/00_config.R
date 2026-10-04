if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
library(data.table)
# Run from 07_scripts; root is portable when the complete package is copied.
root <- normalizePath("..", winslash="/", mustWork=TRUE)
new_root <- "__SMR_ROOT__/GTEX_nonbrain_completion_2026-09-18"
old_map <- fread(file.path(root,"01_QC/old10_input_mapping.tsv"))
new_tissues <- c("Colon_Sigmoid","Colon_Transverse","Minor_Salivary_Gland","Pancreas","Small_Intestine_Terminal_Ileum","Whole_Blood")
tissues <- c(old_map$Tissue,new_tissues)
write_tsv <- function(x,path) fwrite(x,file.path(root,path),sep="\t",na="NA")
