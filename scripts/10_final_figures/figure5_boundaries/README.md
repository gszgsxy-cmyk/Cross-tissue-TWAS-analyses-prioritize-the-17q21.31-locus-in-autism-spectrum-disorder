# Figure 5 core-region boundaries and axis unit

`python annotate_boundaries.py --check-only` checks the archived vector base, six-tick
linear coordinate calibration, preserved primitives and in-memory MB-to-Mb text
correction without writing a figure.
`python annotate_boundaries.py --output-dir output` writes a vector PDF; optional
`--raster` produces a 4500-pixel-wide TIFF (600-DPI metadata) and PNG using PDFium.
Dependencies: numpy, pdfplumber, pypdf, reportlab, Pillow; optional pypdfium2.

The base is the author's existing glyph-compatible vector figure, not GWAS/LD inputs.
Its original MB label and SHA256 are retained unchanged. The output changes the
existing axis text literal to `chr 17 physical position (Mb)` using the same font;
it does not mask, rasterize or redraw the plot to correct the unit.
Boundaries are GRCh38 45.414716/46.836265 Mb, 0.6-pt grey dashed lines, under existing
labels/points. The v2 inter-track connection is retained. No FUSION statistics are run.
The R source in scripts/02_fusion mixes computation and drawing and is not the entry
point for this annotation. Plot text and all 1,056 circles remain vector.
The corrected output was rendered with Poppler and matched the approved TIFF
pixel for pixel. Before/after differences are confined to the unit-label region;
all boundary lines, scatter points and other labels are unchanged.
The optional PDFium raster can differ in antialiasing from historical Poppler TIFFs.
