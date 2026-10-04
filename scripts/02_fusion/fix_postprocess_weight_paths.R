input_file  <- "ASD.17q21.31.test.input.dat"
output_file <- "ASD.17q21.31.test.input.fixed.dat"

d <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ???FUSION??????????????????
all_weights <- list.files(
  ".",
  pattern = "\\.wgt\\.RDat$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Total weight files found:", length(all_weights), "\n")

requested_names <- basename(d$FILE)
available_names <- basename(all_weights)

# ??????????
duplicate_names <- unique(
  available_names[duplicated(available_names)]
)

if (any(requested_names %in% duplicate_names)) {
  stop(
    "Duplicated weight filenames were found: ",
    paste(intersect(requested_names, duplicate_names), collapse = ", ")
  )
}

idx <- match(requested_names, available_names)

check_table <- data.frame(
  ID = d$ID,
  RequestedFile = requested_names,
  Matched = !is.na(idx),
  stringsAsFactors = FALSE
)

print(check_table)

if (any(is.na(idx))) {
  stop(
    "Some weight files were not found: ",
    paste(requested_names[is.na(idx)], collapse = ", ")
  )
}

# ????????
d$FILE <- normalizePath(
  all_weights[idx],
  winslash = "/",
  mustWork = TRUE
)

# ?????????????????????
d$PANEL[is.na(d$PANEL) | d$PANEL == "NA" | d$PANEL == ""] <-
  "GTEx.Whole_Blood"

write.table(
  d,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFixed input written to:", output_file, "\n")
cat("Existing weight files:", sum(file.exists(d$FILE)), "/", nrow(d), "\n")
