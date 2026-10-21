# ASD TWAS Analysis Code

Study-specific code accompanying the ASD transcriptome-wide association study.
Version: 2 October 2026.

## 1. Scope

The collection documents TWAS result processing, regional conditional/joint analyses, SMR analyses, colocalization, external evidence integration, and selected visualizations.

## 2. Analysis modules

| Directory | Purpose |
| --- | --- |
| `scripts/02_fusion` | FUSION conditional/joint plotting and sensitivity-input preparation |
| `scripts/03_twas_summary_figures` | TWAS result processing and exploratory visualization |
| `scripts/04_smr_execution` | SMR execution scripts |
| `scripts/05_nonbrain_rebuild` | Non-brain eQTL-SMR aggregation and tissue-specific FDR processing |
| `scripts/06_coloc` | Final 23-probe BrainMeta full-cis colocalization |
| `scripts/07_mqtl` | Hatton EUR and fetal-brain mQTL-SMR processing and enrichment |
| `scripts/08_external_evidence` | External ASD evidence integration |
| `scripts/09_figures` | Figure 4 panel preparation and assembly |
| `third_party` | FUSION and MR-JTI source files and licenses |
| `config` and `tools` | Path configuration template and helper |

## 3. Analytical specifications

- FUSION conditional/joint analyses use GRCh38/hg38: core interval chr17:45,414,716-46,836,265; extended interval chr17:44,414,716-47,836,265.
- Brain mQTL regional summaries use GRCh37/hg19: chr17:43,492,082-44,913,631.
- The final colocalization analysis comprises 23 probes, with one full-cis analysis per probe, using R 4.4.1 and coloc 5.2.3.
- Hatton EUR significance is defined by genome-wide Benjamini-Hochberg correction across valid SMR tests. Cross-resource CpG comparisons retain the resource-specific significance definitions within the shared testable background.
- Walker developmental-brain evidence provides developmental context and is not counted as ASD replication.

## 4. Input data and software

Required external inputs include GWAS summary statistics, expression prediction weights, LD reference panels, QTL resources, genomic annotations, and published supplementary datasets. These datasets must be supplied separately. Script-specific imports identify the R and Python dependencies; shell scripts require the corresponding analysis executables.

Update path placeholders and working directories before use. The configuration helper can create a separate copy using the roots specified in `config/paths.example.json`. Some scripts additionally require intermediate files and the original analysis directory structure.

## 5. Archive coverage

The archive contains available analysis scripts and sequential processing records. The complete original JTI launch workflow and workbook-formatting steps are not included. Exploratory plotting scripts do not represent the latest author-prepared Figures 2 and 3. End-to-end execution has not been validated for this release.

## 6. Licensing

FUSION and derived FUSION scripts retain the upstream GPLv3 license. MR-JTI/PrediXcan-R retains the upstream MIT license. License texts are provided in `third_party`. A collection-wide license has not been assigned.


## 7. FUSION regional conditional/joint analysis reproducibility

The combined archive supports the FUSION checks described in Reviewer 3 Major Comments 1 and 2 during manuscript revision, covering the core primary (353 models), extended primary (358), core all-testable (559), and extended all-testable (1,011) analyses.

- [Combined reproducibility ZIP](reproducibility/FUSION_GRCh38_Reproducibility.zip) ([direct download](https://raw.githubusercontent.com/gszgsxy-cmyk/Cross-tissue-TWAS-analyses-prioritize-the-17q21.31-locus-in-autism-spectrum-disorder/main/reproducibility/FUSION_GRCh38_Reproducibility.zip)): 810,699 bytes; 58 files.
- Contents include recovered commands, diagnostic commands and parameters, software/script information, 1000 Genomes European-ancestry LD-reference provenance and GRCh38 coordinate checks, prediction-model P0/P1 checks, regional boundaries, selection diagnostics, numerical-instability records, and default versus common absolute-Z threshold (4.058183) comparisons.
- After extraction, start with `README_EN.md` and `Major2_numerical_stability_notes_EN.md`. `SHA256_manifest.tsv` verifies the 57 other files. Global and locus-specific selection are distinguished; the unstable 10/19 locus-specific selections do not establish reliable independent signals.

Recovered historical records are distinguished from diagnostic reruns generated on **16 October 2026**. FUSION diagnostic runs used **R 3.6.0**. The exact historical FUSION Git commit and complete independent historical stdout/stderr were not recovered. Bundled current local and instrumented scripts are identified by checksums; they are not claimed to be byte-identical historical source snapshots.

Complete external GWAS data, prediction weights, QTL resources and LD genotypes are not bundled and must be obtained separately under the providers' access and licensing conditions. The full 2,281-record model-membership TSV and Supplementary Tables workbook are separate article supplements. Companion-workbook identifiers inside the ZIP describe the archive-preparation snapshot; use the final article supplements for the updated workbook.

ZIP SHA256:

```text
f52792140f945eb3d124c367afabc166de5194f58fe3d864ad79de509106293b
```
