# Translational Efficiency Analysis with Ribolog

**Author:** Priyanshu Panda  
**Purpose:** Identify transcripts with **differential translational efficiency (TE)** between treatment (**H**) and control (**N**) using RNA-seq and Ribo-seq count data.

---

## 📌 Overview
This repository contains an R-based workflow that processes RNA-seq and Ribo-seq counts, merges them into a Ribolog-ready matrix, and performs **logit-seq TE differential analysis**.  

The pipeline outputs:  
- Combined RNA+RPF matrix  
- Sample metadata table  
- Differential TE results (BH-adjusted p-values)  
- Volcano plots  
- Lists of significantly upregulated and downregulated transcripts  

---

## ⚙️ Requirements

- R (≥ 4.0)  
- [Ribolog](https://github.com/ohlerlab/Ribolog)  
- tidyverse packages:  
  - `dplyr`  
  - `readr`  
  - `stringr`  
  - `ggplot2`  
  - `purrr`  


## 🚀 Workflow

1. **Load libraries**  
   Import Ribolog + tidyverse dependencies.  

2. **Define file paths**  
   Set base directory and count file locations.  

3. **Merge RNA-seq counts**  
   - Reads all `*_RNA_RNA.out` files.  
   - Removes transcript version numbers.  

4. **Merge Ribo-seq counts**  
   - Reads all `*_Ribo.out` files.  
   - Same preprocessing as RNA.  

5. **Combine RNA + Ribo matrices**  
   - Creates Ribolog-ready matrix (`rr_matrix.tsv`).  

6. **Create sample design matrix**  
   - Defines **condition**: control (**N**) vs treatment (**H**).  
   - Defines **read type**: RNA or RPF.  

7. **Filter + run Ribolog**  
   - Applies minimum count filter.  
   - Runs `logit_seq()` to identify differential TE genes.  

8. **Visualize with volcano plot**  
   - X-axis: log2 TE fold change  
   - Y-axis: -log10 adjusted p-value  
   - Significant transcripts in **red**.  

9. **Extract significant transcripts**  
   - Thresholds: `padj < 0.05` and `|log2FC| > 1`  
   - Saves:  
     - `significant_TE_transcripts.csv`  
     - `TE_upregulated.csv`  
     - `TE_downregulated.csv`  

---

## 📊 Example Volcano Plot

- Gray points → non-significant genes  
- Red points → significant TE changes  
- Vertical dashed lines → log2FC ±1  
- Horizontal dashed line → adj. p-value = 0.05  

---

## 📝 Outputs

- `rr_matrix.tsv` → merged RNA+RPF matrix  
- `sample_info.tsv` → sample metadata  
- `significant_TE_transcripts.csv` → significant TE results  
- `TE_upregulated.csv` → TE up in treatment (H)  
- `TE_downregulated.csv` → TE down in treatment (H)  
- Volcano plot (in R session)  

---

## 🔍 Notes

- If you only analyze **RNA-seq**, skip the Ribo merge and Ribolog run.  
- If you already have a **combined RNA+RPF matrix**, begin from **Step 5**.  
- Interpretation:  
  - **Positive log2FC** → higher TE in treatment (**H**)  
  - **Negative log2FC** → lower TE in treatment (**H**)  



