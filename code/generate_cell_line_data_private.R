#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#generate public data
source(here::here("code", "generate_cell_line_data.R"))

achilles_long <- readRDS(file=here::here("data", paste0(release, "_achilles_long.Rds")))
expression_names <- readRDS(file=here::here("data", paste0(release, "_expression_names.Rds")))
prism_long <- readRDS(file=here::here("data", paste0(release, "_prism_long.Rds")))

#Cell_Line correlation matrix
achilles_cell_line <- 
  achilles_long %>%
  dplyr::left_join(expression_names, by = "X1") %>%
  dplyr::select(-X1, -lineage, -lineage_subtype) 

achilles_cell_line_flip <- 
  tidyr::pivot_wider(achilles_cell_line, names_from = cell_line, values_from = dep_score)

achilles_cell_line_cor <- 
  achilles_cell_line_flip %>%
  dplyr::select(-1) %>% 
  corrr::correlate()

#Stats
achilles_cell_line_cor_long <- 
  achilles_cell_line_cor %>%
  dplyr::rename(fav_cell = 1) %>% 
  tidyr::pivot_longer(cols = -fav_cell, names_to = "cell", values_to = "r2") 

#generate nested df for later
achilles_cell_line_cor_nest <- 
  achilles_cell_line_cor_long %>% 
  dplyr::group_by(fav_cell) %>% 
  tidyr::nest()

virtual_achilles_cell_line <- 
  achilles_cell_line_cor_long %>% 
  dplyr::filter(!is.na(r2)) %>%   
  moderndive::rep_sample_n(size = 1000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r2), max = max(r2), min = min(r2), sd = sd(r2))

mean_virtual_achilles_cell_line <- mean(virtual_achilles_cell_line$mean)
sd_virtual_achilles_cell_line <- mean(virtual_achilles_cell_line$sd)

sd_threshold_cell <- 2

achilles_cell_line_upper <- mean_virtual_achilles_cell_line + sd_threshold_cell*sd_virtual_achilles_cell_line
achilles_cell_line_lower <- mean_virtual_achilles_cell_line - sd_threshold_cell*sd_virtual_achilles_cell_line

#Tables
master_top_table_cell_line <- tibble(
  fav_cell = character(), 
  data = list()
)
master_bottom_table_cell_line <- tibble(
  fav_cell = character(), 
  data = list()
)

make_dep_table_cell_line <- function(dep_table = achilles_cell_line_cor_nest, 
                                     summary_table = expression_names, 
                                     query_cell, 
                                     upper = achilles_cell_line_upper, 
                                     lower = achilles_cell_line_lower, 
                                     top = TRUE) {
  dep <- 
    dep_table %>% 
    dplyr::filter(fav_cell %in% query_cell) %>%
    tidyr::unnest(data) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-fav_cell) %>% 
    {if (top == TRUE) dplyr::filter(., .[[2]] > upper) else dplyr::filter(., .[[2]] < lower)} %>% #mean +/- 3sd
    dplyr::arrange(desc(.[[2]])) %>% #use column index
    dplyr::left_join(summary_table, by = c("cell" = "cell_line")) %>% 
    dplyr::select(-X1) 
  return(dep)
}

#define list
full_cell <- achilles_cell_line_cor_nest %>% pull(fav_cell)
cell_group <- full_cell #(~60' on a laptop); change to sample for testing, or methods

for (fav_cell in cell_group) {
  
  message(" Dep tables for ", fav_cell)
  dep_top_cell <- make_dep_table_cell_line(query_cell = fav_cell, 
                                           top = TRUE)
  dep_top_cell <- 
    dep_top_cell %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles_cell_line)/sd_virtual_achilles_cell_line, 1)) %>% 
    dplyr::select(cell, lineage, lineage_subtype, z_score, r2)
  
  top_table_cell <- dep_top_cell %>% 
    dplyr::mutate(fav_cell = fav_cell) %>% 
    dplyr::group_by(fav_cell) %>% 
    tidyr::nest()
  
  master_top_table_cell_line <- 
    master_top_table_cell_line %>% 
    dplyr::bind_rows(top_table_cell)
  
  dep_bottom_cell <- make_dep_table_cell_line(query_cell = fav_cell, 
                                              top = FALSE)
  
  dep_bottom_cell <-
    dep_bottom_cell %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles_cell_line)/sd_virtual_achilles_cell_line, 1)) %>% 
    dplyr::select(cell, lineage, lineage_subtype, z_score, r2)
  
  bottom_table_cell <- 
    dep_bottom_cell %>% 
    dplyr::mutate(fav_cell = fav_cell) %>% 
    dplyr::group_by(fav_cell) %>% 
    tidyr::nest()
  
  master_bottom_table_cell_line <- 
    master_bottom_table_cell_line %>% 
    dplyr::bind_rows(bottom_table_cell)
}

#Save
#achilles_long (long raw data) is already saved
saveRDS(achilles_cell_line_cor_nest, file = here::here("data", paste0(release, "_achilles_cell_line_cor_nest.Rds")))
#saveRDS(achilles_cell_line_flip, file = here::here("data", paste0(release, "_achilles_cell_line.Rds")))
saveRDS(sd_threshold_cell, file = here::here("data", paste0(release, "_sd_threshold_cell.Rds")))
saveRDS(achilles_cell_line_lower, file = here::here("data", paste0(release, "_achilles_cell_line_lower.Rds")))
saveRDS(achilles_cell_line_upper, file = here::here("data", paste0(release, "_achilles_cell_line_upper.Rds")))
saveRDS(mean_virtual_achilles_cell_line, file = here::here("data", paste0(release, "_mean_virtual_achilles_cell_line.Rds")))
saveRDS(sd_virtual_achilles_cell_line, file = here::here("data", paste0(release, "_sd_virtual_achilles_cell_line.Rds")))

saveRDS(master_top_table_cell_line, file=here::here("data", paste0(release, "_master_top_table_cell_line.Rds")))
saveRDS(master_bottom_table_cell_line, file=here::here("data", paste0(release, "_master_bottom_table_cell_line.Rds")))

##Common Essential Genes----------------------------------------------------------
common_essentials <- read_csv(common_essentials_url, col_names = TRUE) %>%
  mutate(gene = word(gene, 1, -2))

saveRDS(common_essentials, file = here::here("data", paste0(release, "_common_essentials.Rds")))

##Unique Essential Genes-------------------------------------------------------------------

var_unique_essential_genes_cutoff <- n_distinct(achilles_long$X1) * 0.05

unique_essential_genes <- achilles_long %>%
  dplyr::filter(dep_score < -1) %>%
  count(gene) %>%
  filter(n < var_unique_essential_genes_cutoff)

saveRDS(unique_essential_genes, file = here::here("data", paste0(release, "_unique_essential_genes.Rds")))

##Prism------------------------------------------------------------------

var_prism_cell_cutoff <- n_distinct(prism_long$x1) * 0.05

prism_unique_toxic <- prism_long %>% 
  dplyr::filter(log2fc < -1) %>% 
  dplyr::count(name) %>%
  dplyr::filter(n < var_prism_cell_cutoff)

saveRDS(prism_unique_toxic, file = here::here("data", paste0(release, "_prism_unique_toxic.Rds")))

##Metabolites---------------------------------------------------------

temp <- tempfile()
download.file(metabolites_cell_url, temp)
cell_metabolites <- read.csv(temp)
saveRDS(cell_metabolites, file = here::here("data", paste0(release, "_cell_metabolites.Rds")))

