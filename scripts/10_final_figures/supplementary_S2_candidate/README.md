# Supplementary Figure S2 layout candidate

Author adoption of the revised layout in the final supplementary upload is unconfirmed.
This module is a CANDIDATE, not a claim about the final published supplementary figure.
The candidate TIFF is byte-identical to S2 embedded in the locally reviewed revised
supplementary-figure Word document. Older local supplementary files retain the previous
layout; journal-side adoption has not been independently verified. The candidate
directory and its plotting inputs are therefore retained unchanged.
`python assemble_candidate.py --check-only` reads the five-CpG source and archived panels.
To assemble the existing candidate: `python assemble_candidate.py --bold-font PATH_TO_ARIAL_BOLD_TTF`.
The script requires the approved heading font rather than silently substituting a font.
It reuses original vectors; forestploter statistics, points, CI, rows and colours are not
recalculated. Source R/Python records are in archived_sources/ with root placeholders;
they require historical inputs and are NOT the portable execution entry point.
Current entry point dependencies: pypdf and reportlab. No figure was regenerated here.
