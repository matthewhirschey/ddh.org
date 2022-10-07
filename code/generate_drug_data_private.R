#rm(list=ls()) 
time_begin_data <- Sys.time()

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))
source(here::here("code", "generate_drug_data.R")) #script to generate drug names (for search)

#load data
prism_meta <- readRDS(file=here::here("data", paste0(release, "_prism_meta.Rds")))
prism_names <- readRDS(file=here::here("data", paste0(release, "_prism_names.Rds")))
gene_summary <- readRDS(file=here::here("data", paste0(release, "_gene_summary.Rds")))

#calc like depmap
#These are log-fold change collapsed replicates with outliers and controls removed of PRISM
prism_raw <- 
  readr::read_csv(prism_url, col_names = TRUE) %>% 
  janitor::clean_names() %>% 
  dplyr::select(-any_of(censor_ids)) #censor drugs that have gene names #"lta", "nppb", "prima1" defined in generate_drug_data.R

prism_long <- 
  prism_raw %>% #need this for joining below
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "drug", values_to = "log2fc") %>% 
  dplyr::left_join(prism_names, by = c("drug" = "clean_drug")) %>% 
  dplyr::filter(!is.na(name)) %>% #remove drugs which don't have a searchable name
  dplyr::select(x1, name, log2fc)

prism <-
  prism_long %>% 
  tidyr::pivot_wider(names_from = "name", values_from = "log2fc")

prism_cor_long <- 
  prism %>%
  dplyr::select(-1) %>% 
  corrr::correlate() %>% 
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "name", values_to = "r2")

#Permutation tests
virtual_prism_cor <- 
  prism_cor_long %>% 
  dplyr::filter(!is.na(r2)) %>%   
  moderndive::rep_sample_n(size = 4000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r2), max = max(r2), min = min(r2), sd = sd(r2))

mean_virtual_prism_cor <- mean(virtual_prism_cor$mean)
sd_virtual_prism_cor <- mean(virtual_prism_cor$sd)

drug_cor_sd_threshold <- 3

prism_cor_upper <- mean_virtual_prism_cor + drug_cor_sd_threshold*sd_virtual_prism_cor
prism_cor_lower <- mean_virtual_prism_cor - drug_cor_sd_threshold*sd_virtual_prism_cor

prism_cor_nest <- 
  prism_cor_long %>%
  dplyr::rename(fav_drug = 1) %>% 
  dplyr::group_by(fav_drug) %>% 
  tidyr::nest()

#save
saveRDS(prism_long, file = here::here("data", paste0(release, "_prism_long.Rds")))
saveRDS(prism_cor_nest, file = here::here("data", paste0(release, "_prism_cor_nest.Rds")))
saveRDS(drug_cor_sd_threshold, file = here::here("data", paste0(release, "_drug_cor_sd_threshold.Rds")))
saveRDS(prism_cor_lower, file = here::here("data", paste0(release, "_prism_cor_lower.Rds")))
saveRDS(prism_cor_upper, file = here::here("data", paste0(release, "_prism_cor_upper.Rds")))
saveRDS(mean_virtual_prism_cor, file = here::here("data", paste0(release, "_mean_virtual_prism_cor.Rds")))
saveRDS(sd_virtual_prism_cor, file = here::here("data", paste0(release, "_sd_virtual_prism_cor.Rds")))

#-----
#recalc using log-fold
#first, get log2FC values of ACHILLES for integration
achilles_log2fc_raw <- 
  readr::read_csv(achilles_log_url, col_names = TRUE) %>% 
  dplyr::rename("X1" = 1)
achilles_guide_map <- read_csv(achilles_guide_map_url, col_names = TRUE)
achilles_rep_map <- read_csv(achilles_rep_map_url, col_names = TRUE)

#clean
achilles_guide_map$gene <- str_remove_all(achilles_guide_map$gene, "\\s\\(\\d+\\)")

achilles_log2fc <- 
  achilles_guide_map %>% 
  dplyr::left_join(achilles_log2fc_raw, by = c("sgrna" = "X1")) #%>% select(1:100)

achilles_log2fc_long <- 
  achilles_log2fc %>% 
  tidyr::pivot_longer(cols = "143B-311Cas9_RepA_p6_batch3":"BT549-311cas9 Rep A p5_batch2", names_to = "cell_line", values_to = "log2fc")

achilles_log2fc_long <- 
  achilles_log2fc_long %>% 
  dplyr::left_join(achilles_rep_map, by = c( "cell_line" = "replicate_ID")) %>% 
  dplyr::select(DepMap_ID, gene, log2fc)

achilles_log2fc_long_mean <- 
  achilles_log2fc_long %>% 
  dplyr::group_by(DepMap_ID, gene) %>% 
  dplyr::summarize(meanlog2fc = mean(log2fc), 
            sdlog2fc = sd(log2fc)) %>% #add QC for SD that is too high by filter()?
  dplyr::ungroup()

#combine, and pivot_wider for matrix generation
combined <- 
  achilles_log2fc_long_mean %>% 
  dplyr::select(x1 = DepMap_ID, name = gene, log2fc = meanlog2fc) %>% 
  dplyr::bind_rows(prism_long) %>% 
  dplyr::rename(DepMap_ID = x1) %>% 
  tidyr::pivot_wider(names_from = name, values_from = log2fc)

#Combined CORRELATION MATRIX
combined_cor <- 
  combined %>%
  dplyr::select(-1) %>% 
  corrr::correlate() %>% 
  dplyr::rename(rowname = 1)

#test gene corrs ##comment out
# fav_gene <- "TSC1"
# combined_cor %>%
#   focus(!!fav_gene) %>%
#   arrange(desc(.[[2]]))

#resample for stats
#make some long files
combined_cor_long <- 
  combined_cor %>% 
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "r")

#Permutation tests
virtual_prism <- 
  combined_cor_long %>% 
  dplyr::filter(!is.na(r)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_prism <- mean(virtual_prism$mean)
sd_virtual_prism <- mean(virtual_prism$sd)

drug_sd_threshold <- 2

prism_upper <- mean_virtual_prism + drug_sd_threshold*sd_virtual_prism
prism_lower <- mean_virtual_prism - drug_sd_threshold*sd_virtual_prism

#save
saveRDS(drug_sd_threshold, file = here::here("data", paste0(release, "_drug_sd_threshold.Rds")))
saveRDS(prism_lower, file = here::here("data", paste0(release, "_prism_lower.Rds")))
saveRDS(prism_upper, file = here::here("data", paste0(release, "_prism_upper.Rds")))
saveRDS(mean_virtual_prism, file = here::here("data", paste0(release, "_mean_virtual_prism.Rds")))
saveRDS(sd_virtual_prism, file = here::here("data", paste0(release, "_sd_virtual_prism.Rds")))

#cutoff and make tables
gene_drugs_cor_table <- tibble(
  fav_gene = character(), 
  data = list()
)
drug_genes_cor_table <- tibble(
  fav_drug = character(), 
  data = list()
)

#define list
#genes <- sample(names(select(combined_cor, A1BG:ZZZ3)), size = 10) #comment this out
#drugs <- sample(names(select(combined_cor, !rowname:ZZZ3)), size = 10) #comment this out
genes <- names(dplyr::select(combined_cor, A1BG:ZZZ3))
drugs <- names(dplyr::select(combined_cor, !1:ZZZ3))
#18524+4514 == 23038 (same as combined_cor)

#drug table for a gene
for (i in genes) {
  message("Drug tables for ", i)
  gene_top <- 
    combined_cor %>% 
    dplyr::select(1, all_of(i)) %>% 
    dplyr::arrange(desc(.[[2]])) %>% #use column index
    dplyr::filter(rowname %in% drugs, #remove genes
           .[[2]] > prism_upper) %>% #mean +/- 2sd
    dplyr::rename(drug = 1, 
         r2 = 2) %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_prism)/sd_virtual_prism, 1))

  gene_table <- 
    gene_top %>% 
    dplyr::left_join(prism_meta, by = c("drug" = "name")) %>% 
    dplyr::select(drug, moa, z_score, r2) %>% 
    dplyr::mutate(fav_gene = i) %>% 
    dplyr::group_by(fav_gene) %>% 
    tidyr::nest()
  
  gene_drugs_cor_table <- gene_drugs_cor_table %>% 
    dplyr::bind_rows(gene_table)
}

#gene table for a drug query
for (i in drugs) {
  message("Gene tables for ", i)
  drug_top <- 
    combined_cor %>% 
    dplyr::select(1, all_of(i)) %>% 
    dplyr::arrange(desc(.[[2]])) %>% #use column index
    dplyr::filter(rowname %in% genes, #remove drugs
           .[[2]] > prism_upper) %>% #mean +/- 2sd
    dplyr::rename(gene = 1, 
           r2 = 2) %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_prism)/sd_virtual_prism, 1))
  
  drug_table <- drug_top %>% 
    dplyr::left_join(gene_summary, by = c("gene" = "approved_symbol")) %>% 
    dplyr::select(gene, approved_name, z_score, r2) %>% 
    dplyr::mutate(fav_drug = i) %>% 
    dplyr::group_by(fav_drug) %>% 
    tidyr::nest()
  
  drug_genes_cor_table <- 
    drug_genes_cor_table %>% 
    dplyr::bind_rows(drug_table)
}
  
top_100_drug_correlations <-
  combined_cor %>% 
  dplyr::select(1, all_of(drugs)) %>% 
  tidyr::pivot_longer(-rowname, names_to = "drugs", values_to = "r2") %>% 
  dplyr::filter(!rowname %in% drugs) %>% 
  dplyr::arrange(desc(r2)) %>% 
  dplyr::select(drug = drugs, gene = rowname, r2) %>% 
  dplyr::slice(1:100)

#TEST get data out
# make_drug_table <- function(gene_data = gene_drugs_table, gene_symbol) {
#   gene_data %>%
#     dplyr::filter(fav_gene %in% gene_symbol) %>%
#     tidyr::unnest(data) %>%
#     dplyr::arrange(desc(r2)) %>%
#     dplyr::rename("Query" = "fav_gene", "Drug" = "drug", "R^2" = "r2", "Z Score" = "z_score")
# }
# make_gene_table <- function(drug_data = drug_genes_table, drug_name) {
#   drug_data %>%
#     dplyr::filter(fav_drug %in% drug_name) %>%
#     tidyr::unnest(data) %>%
#     dplyr::arrange(desc(r2)) %>%
#     dplyr::rename("Query" = "fav_drug", "Gene" = "gene", "R^2" = "r2", "Z Score" = "z_score")
# }
#combined_cor_long %>% arrange(desc(r)) %>% filter(x %in% drugs) %>% filter(y %in% genes)

#get drugs associated with genes by querying gene
compound_table <- 
  prism_meta %>% 
  dplyr::select(fav_drug = name, fav_gene = target, moa) %>% 
  tidyr::separate_rows(fav_gene, sep = ",") %>% 
  dplyr::mutate(fav_gene = stringr::str_trim(fav_gene, side = "both")) %>% 
  drop_na()

gene_drugs_table <-
  compound_table %>% 
  dplyr::group_by(fav_gene) %>% 
  nest()

#get genes associated with drugs by querying drugs
drug_genes_table <-
  compound_table %>% 
  dplyr::left_join(gene_summary, by = c("fav_gene" = "approved_symbol")) %>% 
  dplyr::select(fav_drug, fav_gene, approved_name) %>% 
  dplyr::group_by(fav_drug) %>% 
  nest()

#save files
saveRDS(prism, file = here::here("data", paste0(release, "_prism.Rds")))
saveRDS(gene_drugs_cor_table, file=here::here("data", paste0(release, "_gene_drugs_cor_table.Rds")))
saveRDS(drug_genes_cor_table, file=here::here("data", paste0(release, "_drug_genes_cor_table.Rds")))
saveRDS(gene_drugs_table, file=here::here("data", paste0(release, "_gene_drugs_table.Rds")))
saveRDS(drug_genes_table, file=here::here("data", paste0(release, "_drug_genes_table.Rds")))
saveRDS(top_100_drug_correlations, file=here::here("data", paste0(release, "_top_100_drug_correlations.Rds")))

#how long
time_end_data <- Sys.time()
