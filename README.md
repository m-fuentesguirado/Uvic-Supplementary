# UVic supplementary materials

This repository contains the supplementary documents, analysis code and one
standalone sequencing workflow for the study.

## Contents

| Folder | Contents |
| --- | --- |
| `analysis/` | R Markdown code for sequencing QC, downstream/clinical analysis and fragmentomics |
| `supplementary/` | Final supplementary figures and supplementary methods PDFs |
| `workflow/fastq_to_dx/` | Standalone Nextflow FASTQ-to-duplex-consensus workflow |

### Analysis files

- `qc_plots_R2.Rmd`: targeted-panel sequencing QC and Supplementary Figures S1-S6.
- `downstream_analysis_R1.Rmd`: variant filtering, ctDNA and clinical analyses, including Supplementary Figures S7-S11.
- `downstream_analysis_R2.Rmd`: fragmentomics analysis and figure generation.
- `unified_figures.R`: shared figure colours and plotting theme.

Patient-level input data, sequencing files and reference-genome files are not
included.
