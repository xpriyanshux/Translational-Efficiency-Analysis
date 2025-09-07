# Translational-Efficiency-Analysis

This workflow processes raw RNA-seq and Ribo-seq count files, merges them into a Ribolog-ready matrix, and runs logit-seq based TE differential analysis. The pipeline generates:

A combined RNA+RPF count matrix

A sample design table (metadata for Ribolog)

Differential TE results with adjusted p-values (BH correction)

Volcano plots for visualization

Lists of significantly upregulated and downregulated transcripts

⚙️ Requirements
Software

R (≥ 4.0)

Ribolog

tidyverse packages
 (dplyr, readr, stringr, ggplot2, purrr)

Input files

RNA-seq count files (*_RNA_RNA.out)

Tab-delimited, containing at least: GeneID and counts

Gene IDs may include version numbers (e.g., ENSG000001.1)

Ribo-seq count files (*_Ribo.out)

Same structure as RNA count files

⚠️ If you already have a combined RNA+RPF count table, you can skip the merging steps.
