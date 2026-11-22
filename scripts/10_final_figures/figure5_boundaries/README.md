# Figure 5 boundaries and axis unit

Annotate the existing vector Figure 5 with GRCh38 core boundaries at
45.414716 and 46.836265 Mb and correct the x-axis unit from MB to Mb.
The input is `data/Figure_5_Conditional_Joint_base.pdf`.

From this directory:

```sh
python annotate_boundaries.py --check-only
python annotate_boundaries.py --output-dir output
python annotate_boundaries.py --output-dir output/raster --raster
```

Dependencies are `numpy`, `pdfplumber`, `pypdf`, `reportlab` and `Pillow`;
`--raster` additionally requires `pypdfium2`. These are listed in the root
`requirements-plotting.txt`.

Outputs are `Figure_5_Core_Boundaries_Connected.pdf` and `annotation_QC.json`.
With `--raster`, the script also writes PNG and a 4500-pixel-wide TIFF with
600-DPI metadata. Choose a new output directory to retain earlier versions.

The script annotates the supplied vector result; it does not reconstruct the
regional analysis from GWAS or LD inputs. The original base is unchanged.
