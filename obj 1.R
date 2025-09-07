
library(Ribolog)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)  # <-- this provides reduce()


# -----------------------------------------------
# STEP 1: Define file paths
# -----------------------------------------------
base_dir <- "D:/New folder/rome/rrrr/ROME ASSIGNMENT"
rna_dir <- file.path(base_dir, "counts_rna")
ribo_dir <- file.path(base_dir, "counts_ribo")

# Output paths
rr_matrix_file <- file.path(base_dir, "rr_matrix.tsv")
sample_info_file <- file.path(base_dir, "sample_info.tsv")

# -----------------------------------------------
# STEP 2: Read and merge RNA counts
# Purpose: Create one matrix of RNA counts from all samples
# -----------------------------------------------
rna_files <- list.files(rna_dir, pattern = "_RNA_RNA.out$", full.names = TRUE)

rna_list <- lapply(rna_files, function(file) {
  df <- read_tsv(file, col_names = TRUE)
  df$GeneID <- sub("\\..*", "", df$GeneID)  # Strip version
  sample_name <- str_extract(basename(file), "IB_\\d+")
  colnames(df)[2] <- paste0(sample_name, "_RNA")
  return(df)
})

rna_merged <- reduce(rna_list, full_join, by = "GeneID")

# -----------------------------------------------
# STEP 3: Read and merge Ribo counts
# Purpose: Create one matrix of Ribo counts from all samples
# -----------------------------------------------
ribo_files <- list.files(ribo_dir, pattern = "_Ribo.out$", full.names = TRUE)

ribo_list <- lapply(ribo_files, function(file) {
  df <- read_tsv(file, col_names = TRUE)
  df$GeneID <- sub("\\..*", "", df$GeneID)
  sample_name <- str_extract(basename(file), "IB_\\d+")
  colnames(df)[2] <- paste0(sample_name, "_RPF")
  return(df)
})

ribo_merged <- reduce(ribo_list, full_join, by = "GeneID")

rna_merged <- rna_merged %>%
  group_by(GeneID) %>%
  summarise(across(everything(), sum), .groups = "drop")

ribo_merged <- ribo_merged %>%
  group_by(GeneID) %>%
  summarise(across(everything(), sum), .groups = "drop")


# -----------------------------------------------
# STEP 4: Combine RNA + Ribo matrices
# Purpose: Create Ribolog-ready matrix of RNA + RPF counts
# -----------------------------------------------
rr_matrix <- full_join(rna_merged, ribo_merged, by = "GeneID") %>%
  dplyr::arrange(GeneID) %>%
  dplyr::filter(rowSums(dplyr::select(., dplyr::ends_with("_RNA"))) > 0)

# Save matrix for use in logit_seq
write_tsv(rr_matrix, rr_matrix_file)

# -----------------------------------------------
# STEP 5: Create sample design matrix
# Purpose: Tell Ribolog which sample belongs to which group
# -----------------------------------------------
sample_ids <- colnames(rr_matrix)[-1]  # exclude GeneID
# Extract just the sample number (01 to 06) using lookbehind
sample_number <- as.numeric(str_extract(sample_ids, "(?<=IB_)(\\d{2})"))

sample_info <- data.frame(
  sample_id = sample_ids,
  condition = ifelse(sample_number <= 3, "N", "H"),
  read_type = ifelse(str_detect(sample_ids, "_RNA$"), "RNA", "RPF")
)


# Save sample info
write_tsv(sample_info, sample_info_file)

# -----------------------------------------------
# STEP 6: Run Ribolog analysis
# Purpose: Perform TE differential testing using logit_seq
# -----------------------------------------------

# Remove any transcript with NA in RNA or RPF counts
rr_matrix <- rr_matrix %>%
  dplyr::filter(!if_any(ends_with("_RNA") | ends_with("_RPF"), is.na))

# Optional: check how many were removed
nrow_before <- nrow(rr_matrix)
cat("Remaining transcripts after NA removal:", nrow_before, "\n")

rr_matrix <- Ribolog::min_count_filter(rr_matrix)

fit <- Ribolog::logit_seq(
  x = as.matrix(rr_matrix[, -1]),       
  design = sample_info,                 
  model = read_type ~ condition,        
  feature_list = rr_matrix$GeneID,      
  adj_method = "BH"                     # <- Required argument
)




library(ggplot2)

fit$padj <- fit$`BH_Pr(>|z|)_conditionH`  # optional, for clarity

ggplot(fit, aes(x = Estimate_conditionH, y = -log10(padj))) +
  geom_point(alpha = 0.4, color = "gray") +
  geom_point(data = subset(fit, padj < 0.05 & abs(Estimate_conditionH) > 1),
             aes(x = Estimate_conditionH, y = -log10(padj)), color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Differential Translational Efficiency (H vs N)",
    x = "log2 TE Fold Change",
    y = "-log10 Adjusted p-value (BH)"
  ) +
  theme_minimal()




# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Add a gene ID column for easy handling (row names aren't always preserved)
fit$GeneID <- rownames(fit)

# Upregulated TE (TE higher in treatment H vs control N)
te_up <- fit %>%
  filter(padj < padj_cutoff & Estimate_conditionH > log2fc_cutoff)

# Downregulated TE (TE lower in treatment H vs control N)
te_down <- fit %>%
  filter(padj < padj_cutoff & Estimate_conditionH < -log2fc_cutoff)

# Combined significant hits
te_sig <- bind_rows(te_up, te_down)

write.csv(te_sig, "significant_TE_transcripts.csv", row.names = FALSE)
write.csv(te_up, "TE_upregulated.csv", row.names = FALSE)
write.csv(te_down, "TE_downregulated.csv", row.names = FALSE)



