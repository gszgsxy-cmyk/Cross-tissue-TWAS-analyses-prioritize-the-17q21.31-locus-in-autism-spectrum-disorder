options(stringsAsFactors = FALSE, warn = 1)

project_root <- "__USER_HOME__/Desktop/TWAS/FUSION/fusion_twas-master"
source_file <- file.path(
  project_root,
  "ASD.FUSION.all29.GTExv8WB.updated.GENCODEv26.allchr.withFDR.dat"
)
output_dir <- file.path(project_root, "POST_GRCh38_ALLMODEL_SENSITIVITY")

regions <- list(
  core = c(start = 45414716, end = 46836265),
  extended = c(start = 44414716, end = 47836265)
)

input_paths <- c(
  core = file.path(
    project_root,
    "ASD.FUSION.17q21.31.GRCh38.core.allmodels.postprocess.dat"
  ),
  extended = file.path(
    project_root,
    "ASD.FUSION.17q21.31.GRCh38.extended.allmodels.postprocess.dat"
  )
)

standard_columns <- c(
  "PANEL", "FILE", "ID", "CHR", "P0", "P1", "HSQ",
  "BEST.GWAS.ID", "BEST.GWAS.Z", "EQTL.ID", "EQTL.R2", "EQTL.Z",
  "EQTL.GWAS.Z", "NSNP", "NWGT", "MODEL", "MODELCV.R2",
  "MODELCV.PV", "TWAS.Z", "TWAS.P"
)

if (!file.exists(source_file)) {
  stop("Authoritative FUSION result file not found: ", source_file)
}

planned_files <- c(
  input_paths,
  file.path(output_dir, "ASD_FUSION_17q21.31_allmodel_input_qc.tsv"),
  file.path(output_dir, "ASD_FUSION_17q21.31_allmodel_weight_resolution.tsv"),
  file.path(output_dir, "ASD_FUSION_17q21.31_allmodel_excluded_models.tsv")
)
existing_planned <- planned_files[file.exists(planned_files)]
if (length(existing_planned) > 0) {
  stop(
    "Refusing to overwrite existing all-model sensitivity files: ",
    paste(existing_planned, collapse = " | ")
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}

cat("Reading authoritative FUSION results...\n")
dat <- read.delim(
  source_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "NaN", "")
)

missing_columns <- setdiff(standard_columns, names(dat))
if (length(missing_columns) > 0) {
  stop("Missing required FUSION columns: ", paste(missing_columns, collapse = ", "))
}

dat$CHR_NUM <- suppressWarnings(as.integer(dat$CHR))
dat$P0_NUM <- suppressWarnings(as.numeric(dat$P0))
dat$P1_NUM <- suppressWarnings(as.numeric(dat$P1))
dat$TWAS_Z_NUM <- suppressWarnings(as.numeric(dat$TWAS.Z))
dat$TWAS_P_NUM <- suppressWarnings(as.numeric(dat$TWAS.P))

gene_symbol <- function(x) {
  candidates <- c("GENE_SYMBOL_GENCODE", "GENE_SYMBOL", "ENSEMBL_ID", "ID")
  out <- rep(NA_character_, nrow(x))
  for (nm in candidates) {
    if (!(nm %in% names(x))) next
    value <- trimws(as.character(x[[nm]]))
    usable <- is.na(out) & !is.na(value) & nzchar(value) & value != "NA"
    out[usable] <- value[usable]
  }
  out
}

panel_cache <- new.env(parent = emptyenv())

panel_weight_index <- function(panel) {
  key <- make.names(panel)
  if (exists(key, envir = panel_cache, inherits = FALSE)) {
    return(get(key, envir = panel_cache, inherits = FALSE))
  }
  panel_dir <- file.path(project_root, "WEIGHTS", panel)
  paths <- if (dir.exists(panel_dir)) {
    list.files(
      panel_dir,
      pattern = "\\.wgt\\.RDat$",
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    character(0)
  }
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  index <- split(paths, basename(paths))
  assign(key, index, envir = panel_cache)
  index
}

resolve_weight <- function(file_entry, panel) {
  file_entry <- as.character(file_entry)
  panel <- as.character(panel)
  clean <- gsub("\\\\", "/", file_entry)
  relative <- sub("^\\./", "", clean)
  base <- basename(clean)

  direct_candidates <- unique(c(
    if (grepl("^[A-Za-z]:/", clean)) clean else character(0),
    file.path(project_root, relative),
    file.path(project_root, "WEIGHTS", panel, base),
    file.path(project_root, "WEIGHTS", panel, panel, base)
  ))
  direct_hits <- direct_candidates[file.exists(direct_candidates)]
  if (length(direct_hits) > 0) {
    return(list(
      status = "found",
      path = normalizePath(direct_hits[1], winslash = "/", mustWork = TRUE),
      match_count = length(unique(normalizePath(
        direct_hits,
        winslash = "/",
        mustWork = TRUE
      )))
    ))
  }

  index <- panel_weight_index(panel)
  recursive_hits <- index[[base]]
  if (is.null(recursive_hits)) recursive_hits <- character(0)
  recursive_hits <- unique(recursive_hits)
  if (length(recursive_hits) == 1) {
    return(list(status = "found", path = recursive_hits, match_count = 1L))
  }
  if (length(recursive_hits) > 1) {
    return(list(
      status = "ambiguous",
      path = paste(recursive_hits, collapse = " | "),
      match_count = length(recursive_hits)
    ))
  }
  list(status = "missing", path = NA_character_, match_count = 0L)
}

qc_rows <- list()
resolution_rows <- list()
excluded_rows <- list()

for (region_name in names(regions)) {
  bounds <- regions[[region_name]]
  in_region <- (
    dat$CHR_NUM == 17L &
    is.finite(dat$P0_NUM) &
    is.finite(dat$P1_NUM) &
    dat$P0_NUM <= bounds[["end"]] &
    dat$P1_NUM >= bounds[["start"]]
  )
  regional <- dat[in_region, , drop = FALSE]

  valid_stats <- (
    is.finite(regional$TWAS_Z_NUM) &
    is.finite(regional$TWAS_P_NUM) &
    regional$TWAS_P_NUM >= 0 &
    regional$TWAS_P_NUM <= 1
  )
  candidates <- regional[valid_stats, , drop = FALSE]

  resolution <- lapply(
    seq_len(nrow(candidates)),
    function(i) resolve_weight(candidates$FILE[i], candidates$PANEL[i])
  )
  resolution_status <- vapply(resolution, `[[`, character(1), "status")
  resolution_path <- vapply(
    resolution,
    function(x) if (is.na(x$path)) "" else x$path,
    character(1)
  )
  resolution_count <- vapply(resolution, `[[`, integer(1), "match_count")

  resolution_table <- data.frame(
    Analysis = region_name,
    PANEL = candidates$PANEL,
    ID = candidates$ID,
    Gene = gene_symbol(candidates),
    Original_FILE = candidates$FILE,
    Resolution_status = resolution_status,
    Resolution_match_count = resolution_count,
    Resolved_FILE = resolution_path,
    stringsAsFactors = FALSE
  )
  resolution_rows[[region_name]] <- resolution_table

  include <- resolution_status == "found"
  prepared <- candidates[include, standard_columns, drop = FALSE]
  prepared$FILE <- resolution_path[include]

  if (anyDuplicated(paste(prepared$PANEL, prepared$FILE, sep = "\t")) > 0) {
    stop("Duplicate PANEL+FILE keys detected in ", region_name, " input")
  }

  write.table(
    prepared,
    file = input_paths[[region_name]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  invalid_table <- regional[!valid_stats, , drop = FALSE]
  invalid_excluded <- data.frame(
    Analysis = rep(region_name, nrow(invalid_table)),
    PANEL = invalid_table$PANEL,
    ID = invalid_table$ID,
    Gene = gene_symbol(invalid_table),
    FILE = invalid_table$FILE,
    Exclusion_reason = "Invalid or unavailable TWAS.Z/TWAS.P",
    stringsAsFactors = FALSE
  )
  weight_excluded_base <- candidates[!include, , drop = FALSE]
  weight_excluded <- data.frame(
    Analysis = rep(region_name, nrow(weight_excluded_base)),
    PANEL = weight_excluded_base$PANEL,
    ID = weight_excluded_base$ID,
    Gene = gene_symbol(weight_excluded_base),
    FILE = weight_excluded_base$FILE,
    Exclusion_reason = paste0(
      "Prediction weight ", resolution_status[!include]
    ),
    stringsAsFactors = FALSE
  )
  excluded_rows[[region_name]] <- rbind(invalid_excluded, weight_excluded)

  prepared_gene <- gene_symbol(candidates[include, , drop = FALSE])
  duplicate_key <- paste(prepared$PANEL, prepared$FILE, sep = "\t")
  p_equal_one <- sum(candidates$TWAS_P_NUM[include] == 1)

  qc_rows[[region_name]] <- data.frame(
    Analysis = region_name,
    Region_start_GRCh38 = bounds[["start"]],
    Region_end_GRCh38 = bounds[["end"]],
    Regional_models = nrow(regional),
    Valid_TWAS_models = nrow(candidates),
    Input_models = nrow(prepared),
    Unique_genes = length(unique(prepared_gene)),
    Tissues = length(unique(prepared$PANEL)),
    KANSL1_models = sum(prepared_gene == "KANSL1", na.rm = TRUE),
    KANSL1_AS1_models = sum(prepared_gene == "KANSL1-AS1", na.rm = TRUE),
    MAPT_AS1_models = sum(prepared_gene == "MAPT-AS1", na.rm = TRUE),
    Duplicate_PANEL_FILE_rows = sum(duplicated(duplicate_key)),
    Missing_TWAS_Z_or_P = sum(!valid_stats),
    Weight_files_found = sum(include),
    Weight_files_missing_or_ambiguous = sum(!include),
    TWAS_P_equal_1 = p_equal_one,
    Input_file = input_paths[[region_name]],
    stringsAsFactors = FALSE
  )

  cat(
    region_name,
    ": regional=", nrow(regional),
    ", valid=", nrow(candidates),
    ", input=", nrow(prepared),
    ", missing_weights=", sum(!include),
    "\n",
    sep = ""
  )
}

qc <- do.call(rbind, qc_rows)
resolution_all <- do.call(rbind, resolution_rows)
excluded_all <- do.call(rbind, excluded_rows)

write.table(
  qc,
  file = file.path(output_dir, "ASD_FUSION_17q21.31_allmodel_input_qc.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  resolution_all,
  file = file.path(
    output_dir,
    "ASD_FUSION_17q21.31_allmodel_weight_resolution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  excluded_all,
  file = file.path(
    output_dir,
    "ASD_FUSION_17q21.31_allmodel_excluded_models.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\nInput QC:\n")
print(qc, row.names = FALSE)
cat("\nPreparation complete. FUSION has not been run by this script.\n")
