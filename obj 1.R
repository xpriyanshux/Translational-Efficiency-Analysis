############################################################
# Translational Efficiency Analysis with Ribolog
# Author: Priyanshu Panda
# Purpose: Identify transcripts with differential translational efficiency
#          in treatment (H) vs control (N) using RNA-seq and Ribo-seq counts
############################################################

# -----------------------------------------------
# STEP 0: Load libraries
# -----------------------------------------------
library(Ribolog)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)   # for reduce()

# -----------------------------------------------
# STEP 1: Define file paths
# -----------------------------------------------
base_dir <- "D:/New folder/rome/rrrr/ROME ASSIGNMENT"
rna_dir  <- file.path(base_dir, "counts_rna")
ribo_dir <- file.path(base_dir, "counts_ribo")

# Output paths
rr_matrix_path   <- file.path(base_dir, "rr_matrix.tsv")
sample_info_path <- file.path(base_dir, "sample_info.tsv")

# -----------------------------------------------
# STEP 2: Read and merge RNA counts
# -----------------------------------------------
# Each RNA-seq count file is read, GeneID version is removed,
# and the count column is renamed with the sample ID.
# ⚠️ Note: If you already have a combined RNA matrix, 
#          you can skip this merging step.
rna_files <- list.files(rna_dir, pattern = "_RNA_RNA.out$", full.names = TRUE)

rna_list <- lapply(rna_files, function(file) {
  df <- read_tsv(file, col_names = TRUE)
  df$GeneID <- sub("\\..*", "", df$GeneID)  # remove transcript version
  sample_id <- str_extract(basename(file), "IB_\\d+")
  colnames(df)[2] <- paste0(sample_id, "_RNA")
  return(df)
})

rna_merged <- reduce(rna_list, full_join, by = "GeneID") %>%
  group_by(GeneID) %>%
  summarise(across(everything(), sum), .groups = "drop")

# -----------------------------------------------
# STEP 3: Read and merge Ribo counts
# -----------------------------------------------
# Same as RNA above, but for Ribo-seq counts.
# ⚠️ Not everyone will need this step:
#     - If you are analyzing *only* RNA data, skip it.
#     - If you already have RNA+RPF in one table, skip both merge steps.
ribo_files <- list.files(ribo_dir, pattern = "_Ribo.out$", full.names = TRUE)

ribo_list <- lapply(ribo_files, function(file) {
  df <- read_tsv(file, col_names = TRUE)
  df$GeneID <- sub("\\..*", "", df$GeneID)
  sample_id <- str_extract(basename(file), "IB_\\d+")
  colnames(df)[2] <- paste0(sample_id, "_RPF")
  return(df)
})

ribo_merged <- reduce(ribo_list, full_join, by = "GeneID") %>%
  group_by(GeneID) %>%
  summarise(across(everything(), sum), .groups = "drop")

# -----------------------------------------------
# STEP 4: Combine RNA + Ribo matrices
# -----------------------------------------------
# Ribolog expects one combined matrix: RNA and RPF counts
# for each sample in the same table.
rr_matrix <- full_join(rna_merged, ribo_merged, by = "GeneID") %>%
  arrange(GeneID) %>%
  filter(rowSums(select(., ends_with("_RNA"))) > 0)

# Save Ribolog-ready matrix
write_tsv(rr_matrix, rr_matrix_path)

# -----------------------------------------------
# STEP 5: Create sample design matrix
# -----------------------------------------------
# This tells Ribolog:
#   - which samples are control (N) vs treatment (H)
#   - whether the column is RNA or RPF
sample_ids    <- colnames(rr_matrix)[-1]  # exclude GeneID
sample_number <- as.numeric(str_extract(sample_ids, "(?<=IB_)(\\d{2})"))

sample_info <- data.frame(
  sample_id = sample_ids,
  condition = ifelse(sample_number <= 3, "N", "H"),
  read_type = ifelse(str_detect(sample_ids, "_RNA$"), "RNA", "RPF")
)

write_tsv(sample_info, sample_info_path)

# -----------------------------------------------
# STEP 6: Filter matrix and run Ribolog
# -----------------------------------------------
# Filter out rows with NA counts and apply Ribolog's min count filter.
rr_matrix <- rr_matrix %>%
  filter(!if_any(ends_with("_RNA") | ends_with("_RPF"), is.na))

cat("Remaining transcripts after NA removal:", nrow(rr_matrix), "\n")

rr_matrix <- Ribolog::min_count_filter(rr_matrix)

# Run TE differential analysis
fit <- Ribolog::logit_seq(
  x           = as.matrix(rr_matrix[, -1]),
  design      = sample_info,
  model       = read_type ~ condition,
  feature_list= rr_matrix$GeneID,
  adj_method  = "BH"
)

# -----------------------------------------------
# STEP 7: Volcano plot
# -----------------------------------------------
fit$padj <- fit$`BH_Pr(>|z|)_conditionH`

ggplot(fit, aes(x = Estimate_conditionH, y = -log10(padj))) +
  geom_point(alpha = 0.4, color = "gray") +
  geom_point(
    data = subset(fit, padj < 0.05 & abs(Estimate_conditionH) > 1),
    aes(x = Estimate_conditionH, y = -log10(padj)),
    color = "red"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Differential Translational Efficiency (H vs N)",
    x = "log2 TE Fold Change",
    y = "-log10 Adjusted p-value (BH)"
  ) +
  theme_minimal()

# -----------------------------------------------
# STEP 8: Extract significant transcripts
# -----------------------------------------------
padj_cutoff   <- 0.05
log2fc_cutoff <- 1

fit$GeneID <- rownames(fit)

te_up <- fit %>%
  filter(padj < padj_cutoff & Estimate_conditionH > log2fc_cutoff)

te_down <- fit %>%
  filter(padj < padj_cutoff & Estimate_conditionH < -log2fc_cutoff)

te_sig <- bind_rows(te_up, te_down)

# Save results
write.csv(te_sig,   "significant_TE_transcripts.csv", row.names = FALSE)
write.csv(te_up,    "TE_upregulated.csv", row.names = FALSE)
write.csv(te_down,  "TE_downregulated.csv", row.names = FALSE)
