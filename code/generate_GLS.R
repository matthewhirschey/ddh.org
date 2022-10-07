
# load libraries
library(tidyverse)
library(here)
library(impute)
library(reticulate)
library(biomaRt)

# read current release information to set parameters for download
source(here::here("code", "current_release.R"))

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
achilles_raw <- read_csv(achilles_url, col_names = TRUE)

sample_info <- readRDS(file = here::here("data", paste0(release, "_expression_meta.Rds"))) %>%
  dplyr::filter(X1 %in% achilles_raw$DepMap_ID) %>%
  dplyr::rename(DepMap_ID = X1,
                CCLE_name = ccle_name)

# KNN imputation
DepMap_ID <- achilles_raw$DepMap_ID

achilles_raw <- t(achilles_raw[,-1])
achilles_raw_knn <- impute::impute.knn(achilles_raw)
achilles_raw <- t(achilles_raw_knn$data)

# add names again
achilles_raw <- achilles_raw %>%
  as.data.frame() %>%
  dplyr::mutate(DepMap_ID = DepMap_ID) %>%
  dplyr::select(DepMap_ID, dplyr::everything())

# get GO-MF "olfactory receptor activity" genes
go_term <- "GO:0004984"
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

olfactory_genes <- getBM(attributes = "hgnc_symbol",
                         filters = "go",
                         values = go_term, 
                         mart = ensembl)

# save files
readr::write_csv(achilles_raw, file = here::here("data/gene_effect.csv"))
readr::write_csv(sample_info, file = here::here("data/sample_info.csv"))
write.table(olfactory_genes, file = here::here("data/olfactory_genes.txt"), row.names = FALSE)

# Run Python scripts
reticulate::use_virtualenv("gls_regression")
reticulate::py_run_file(here::here("code/gene_pairs.py"))

# Acsess objects created in Python from R
GLS_p <- py$GLS_p
genes <- read.table(here::here("data/genes.txt"))

# Format data
rownames(GLS_p) <- genes$V1
colnames(GLS_p) <- genes$V1

GLS_p_long <- GLS_p %>%
  as.table() %>%
  as.data.frame() %>%
  dplyr::filter(Var1 != Var2) %>%
  drop_na() %>%
  dplyr::rename(gene1 = Var1,
                gene2 = Var2,
                pval = Freq) %>%
  dplyr::mutate(fdr = p.adjust(pval))

saveRDS(GLS_p_long, file = here::here("data", paste0(release, "_GLS_p_long.Rds")))

