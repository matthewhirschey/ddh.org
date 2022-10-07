
library(tidyverse)
library(impute)
library(uwot)

# pca <- FALSE # Use UMAP for dimension reduction
pca <- TRUE # Use PCA for dimension reduction

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

# Load cell line metadata
cell_metadata <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_meta.Rds"))) %>% 
  dplyr::select(X1, cell_line, sex, primary_or_metastasis, age, lineage, lineage_subtype)

# DEPENDENCY --------------------------
# Load dependency data
achilles_long <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles_long.Rds")))

## Transpose matrix
achilles_wide <- achilles_long %>% 
  pivot_wider(names_from = gene, values_from = dep_score) %>% 
  column_to_rownames("X1")

## KNN imputation
# sum(apply(achilles_wide, 2, function(x){sum(is.na(x))}))

achilles_wide_trans <- t(achilles_wide)
achilles_wide_trans_knn <- impute::impute.knn(achilles_wide_trans)
achilles_wide_knn <- t(achilles_wide_trans_knn$data)

## Normalization - Autoscaling
achilles_wide_norm <- apply(achilles_wide_knn, 2, function(x) (x - mean(x, na.rm = FALSE))/sd(x, na.rm = FALSE))
  
## Dimension reduction (PCA/UMAP)
if(pca) {
  dep_reduced <- prcomp(achilles_wide_norm, scale. = FALSE, center = FALSE)
  dep_reduced <- dep_reduced$x
  
  } else {
  dep_reduced <- uwot::umap(achilles_wide_norm,
                            n_neighbors = 15, 
                            n_components = sqrt(nrow(achilles_wide_norm))
                            )
}

## Transpose PCA/UMAP components/embeddings
dep_reduced_t <- dep_reduced %>% 
  t() %>% 
  as.data.frame() %>% 
  janitor::clean_names()

# Define cell line names
cell_lines <- colnames(dep_reduced_t)

cell_line_combs <- data.frame(comb = apply(combn(cell_lines, 2), 2, paste, collapse = ";")) %>% 
  separate(comb, into = c("cell1", "cell2"), sep = ";")

# ONLY FOR TESTING
# cell_line_combs <- cell_line_combs %>%
#   dplyr::filter(cell1 == "ach_000739" | cell2 == "ach_000739" ) # HEPG2

# Run models
cell_model_p <- list()
for(i in 1:nrow(cell_line_combs)) {
  
  cell1 <- cell_line_combs %>%
    dplyr::slice(i) %>% 
    dplyr::pull(cell1)
  
  cell2 <- cell_line_combs %>% 
    dplyr::slice(i) %>% 
    dplyr::pull(cell2)

  # LM
  part1 <- noquote(cell1)
  part2 <- " ~ "
  frag1 <- paste(part1, part2)
  part3 <- noquote(cell2)
  frag2 <- paste0(frag1, part3)
  fmla <- as.formula(frag2)
  
  cell_lm <- lm(formula = fmla, data = dep_reduced_t)
  cell_lm_sum <- summary(cell_lm)
  cell_lm_p <- cell_lm_sum$coefficients[2, 4]
  cell_lm_coef <- cell_lm_sum$coefficients[2, 1]
  
  # FINAL RESULT
  cell_model_p[[i]] <- data.frame(cell1 = cell1,
                                  cell2 = cell2,
                                  coef = cell_lm_coef,
                                  pval = cell_lm_p)
  
  # TRACK LOOP
  print(paste0("Done ", i, " out of ", nrow(cell_line_combs)))
  
}

# Extract p-values
cell_line_dep_sim <- cell_model_p %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(fdr = p.adjust(pval, method = "fdr"),
                bonferroni = p.adjust(pval, method = "bonferroni"),
                cell1 = gsub("ach_", "ACH-", cell1),
                cell2 = gsub("ach_", "ACH-", cell2)
  ) %>% 
  dplyr::left_join(cell_metadata, by = c("cell1" = "X1")) %>% 
  dplyr::left_join(cell_metadata, by = c("cell2" = "X1")) %>% 
  dplyr::rename(cell1_name = cell_line.x,
                cell2_name = cell_line.y) %>% 
  dplyr::relocate(cell1_name, .after = cell1) %>% 
  dplyr::relocate(cell2_name, .after = cell2)

saveRDS(cell_line_dep_sim, file = here::here(app_data_dir, paste0(release, "_cell_line_dep_sim.Rds")))

# GENE EXPRESSION --------------------------
# Load expression data
expression_long <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_long.Rds")))

## Transpose matrix
gene_expression_wide <- expression_long %>%
  dplyr::select(-protein_expression) %>%
  pivot_wider(names_from = gene, values_from = gene_expression) %>%
  column_to_rownames("X1") %>%
  filter(rowSums(is.na(.)) != ncol(.))

## Data cleaning
### Remove columns that only have NAs
gene_expression_wide <- gene_expression_wide[, apply(gene_expression_wide, 2, function(x) !all(is.na(x)))]

### Remove columns that only have zeros
gene_expression_wide <- gene_expression_wide[, apply(gene_expression_wide, 2, function(x) !all(x == 0, na.rm = TRUE))]

### Remove columns with var = 0
gene_expression_wide <- gene_expression_wide[, !apply(gene_expression_wide, 2, function(x){var(x, na.rm = TRUE)}) == 0]

## KNN imputation
# sum(apply(gene_expression_wide, 2, function(x){sum(is.na(x))}))

gene_expression_wide_trans <- t(gene_expression_wide)
gene_expression_wide_trans_knn <- impute::impute.knn(gene_expression_wide_trans)
gene_expression_wide_knn <- t(gene_expression_wide_trans_knn$data)

## Normalization - Autoscaling
gene_expression_wide_norm <- apply(gene_expression_wide_knn, 2, function(x) (x - mean(x, na.rm = FALSE))/sd(x, na.rm = FALSE))

## Dimension reduction (PCA/UMAP)
if(pca) {
  exp_reduced <- prcomp(gene_expression_wide_norm, scale. = FALSE, center = FALSE)
  exp_reduced <- exp_reduced$x
  
} else {
  exp_reduced <- uwot::umap(gene_expression_wide_norm,
                            n_neighbors = 15, 
                            n_components = sqrt(nrow(gene_expression_wide_norm))
  )
}

## Transpose PCA/UMAP components/embeddings
exp_reduced_t <- exp_reduced %>% 
  t() %>% 
  as.data.frame() %>% 
  janitor::clean_names()

# Define cell line names
cell_lines <- colnames(exp_reduced_t)

cell_line_combs <- data.frame(comb = apply(combn(cell_lines, 2), 2, paste, collapse = ";")) %>% 
  separate(comb, into = c("cell1", "cell2"), sep = ";")

# ONLY FOR TESTING
# cell_line_combs <- cell_line_combs %>%
#   dplyr::filter(cell1 == "ach_000739" | cell2 == "ach_000739" ) # HEPG2

# Run models
cell_model_p <- list()
for(i in 1:nrow(cell_line_combs)) {
  
  cell1 <- cell_line_combs %>%
    dplyr::slice(i) %>% 
    dplyr::pull(cell1)
  
  cell2 <- cell_line_combs %>% 
    dplyr::slice(i) %>% 
    dplyr::pull(cell2)
  
  # LM
  part1 <- noquote(cell1)
  part2 <- " ~ "
  frag1 <- paste(part1, part2)
  part3 <- noquote(cell2)
  frag2 <- paste0(frag1, part3)
  fmla <- as.formula(frag2)
  
  cell_lm <- lm(formula = fmla, data = exp_reduced_t)
  cell_lm_sum <- summary(cell_lm)
  cell_lm_p <- cell_lm_sum$coefficients[2, 4]
  cell_lm_coef <- cell_lm_sum$coefficients[2, 1]
  
  # FINAL RESULT
  cell_model_p[[i]] <- data.frame(cell1 = cell1,
                                  cell2 = cell2,
                                  coef = cell_lm_coef,
                                  pval = cell_lm_p)
  
  # TRACK LOOP
  print(paste0("Done ", i, " out of ", nrow(cell_line_combs)))
  
}

# Extract p-values
cell_line_exp_sim <- cell_model_p %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(fdr = p.adjust(pval, method = "fdr"),
                bonferroni = p.adjust(pval, method = "bonferroni"),
                cell1 = gsub("ach_", "ACH-", cell1),
                cell2 = gsub("ach_", "ACH-", cell2)
  ) %>% 
  dplyr::left_join(cell_metadata, by = c("cell1" = "X1")) %>% 
  dplyr::left_join(cell_metadata, by = c("cell2" = "X1")) %>% 
  dplyr::rename(cell1_name = cell_line.x,
                cell2_name = cell_line.y) %>% 
  dplyr::relocate(cell1_name, .after = cell1) %>% 
  dplyr::relocate(cell2_name, .after = cell2)

saveRDS(cell_line_exp_sim, file = here::here(app_data_dir, paste0(release, "_cell_line_exp_sim.Rds")))

