# Validation

Validated on 2026-11-11 using R 4.4.1, ggplot2 3.5.1, and circlize 0.4.16.

- The standalone script ran successfully from outside the package directory.
- All 12 bundled data/order files are byte-identical to the verified inputs;
  their SHA256 values are recorded in data_manifest.tsv.
- All four generated panels and the shared legend have identical decoded PDF
  drawing streams, page dimensions, and extracted text to the approved version.
  PDF metadata such as creation time is not part of this comparison.
- The script validates unique complete gene/tissue grids and existing stars
  against the stored group-specific FDR fields. It does not recalculate FDR.
- Counts match the table in README.md; missing Z values remain missing.
- The color range remains -6.1 to 6.1, with no changes to Z or star values.
- No analysis, association testing, model selection, or manuscript editing ran.
- Archive contents were scanned for personal usernames, absolute local paths,
  email addresses, and common credential patterns. No matches were found.
- Generated outputs, runtime logs, saved R workspaces, and build helpers are
  excluded from the upload archive.
