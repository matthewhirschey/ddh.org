library(tidyverse)
source(here::here("code", "current_release.R"))

pathway_go <-  "1902965"
tests_data_dir <- "tests/data"
no_data_gene <- "WASH7P"

# create tests/data directory if it doesn't exist
dir.create(here::here(tests_data_dir), showWarnings = FALSE, recursive = TRUE)

#READ DATA FOR IDS------
message("Generating data")
pathways_filename <- paste0(release, "_pathways.Rds")
pathways_orig <- readRDS(here::here("data", pathways_filename))
pathways <- pathways_orig %>% filter(go == pathway_go)
pathway_genes <- pathways %>% pull(data) %>% first() %>% pull(gene)

master_bottom_table_filename <- paste0(release, "_master_bottom_table.Rds")
master_bottom_table_orig <- readRDS(file=here::here("data", master_bottom_table_filename))

master_top_table_filename <- paste0(release, "_master_top_table.Rds")
master_top_table_orig <- readRDS(file=here::here("data", master_top_table_filename))

#GET TEST DATA IDS------
mbt_genes <- master_bottom_table_orig %>%
  filter(fav_gene %in% pathway_genes) %>%
  pull(data) %>%
  map("gene") %>%
  unlist()

mtt_genes <- master_top_table_orig %>%
  filter(fav_gene %in% pathway_genes) %>%
  pull(data) %>%
  map("gene") %>%
  unlist()

all_genes <- pathway_genes %>%
  append(mbt_genes) %>%
  append(mtt_genes) %>%
  append(c(no_data_gene))

all_genes_and_x1 <- append(c("X1"), all_genes)

#FILTER DATA -----
gene_summary_filename <- paste0(release, "_gene_summary.Rds")
gene_summary <- readRDS(here::here("data", gene_summary_filename)) %>% 
  filter(approved_symbol %in% all_genes)

achilles_filename <- paste0(release, "_achilles.Rds")
achilles <- readRDS(file=here::here("data", achilles_filename)) %>% 
  select(any_of(all_genes_and_x1))

# achilles_cor_filename <- paste0(release, "_achilles_cor.Rds")
# achilles_cor <- readRDS(file=here::here("data", achilles_cor_filename)) %>% 
#   select(any_of(all_genes_and_x1))

achilles_cor_nest_filename <- paste0(release, "_achilles_cor_nest.Rds")
achilles_cor_nest <- readRDS(file=here::here("data", achilles_cor_nest_filename)) %>% 
  filter(fav_gene %in% all_genes)

expression_filename <- paste0(release, "_expression.Rds")
expression <- readRDS(file=here::here("data", expression_filename)) %>% 
  select(any_of(all_genes_and_x1))

expression_meta_filename <- paste0(release, "_expression_meta.Rds")
expression_meta <- readRDS(file=here::here("data", expression_meta_filename)) #no select or filter

expression_names_filename <- paste0(release, "_expression_names.Rds")
expression_names <- readRDS(file=here::here("data", expression_names_filename)) #no select or filter

proteins_filename <- paste0(release, "_proteins.Rds")
proteins <- readRDS(file=here::here("data", proteins_filename)) %>%
  filter(gene_name %in% all_genes)

master_bottom_table <- master_bottom_table_orig %>%
  filter(fav_gene %in% all_genes)

master_top_table <- master_top_table_orig %>%
  filter(fav_gene %in% all_genes)

# Filter the sub-tables of the master tables to only include genes within the all_genes list
temp_bottom_table = tibble(fav_gene = as.character(), data = c())
for(geneName in master_bottom_table$fav_gene){
  filteredData <- master_bottom_table %>% 
    filter(fav_gene == geneName) %>% 
    select(-fav_gene) %>%
    unnest(cols = c(data)) %>% 
    filter(gene %in% all_genes) %>% 
    nest(data=everything())
  newRow <- tibble(fav_gene=geneName, filteredData)
  temp_bottom_table <- temp_bottom_table %>% 
    bind_rows(newRow)
}

temp_top_table = tibble(fav_gene = as.character(), data = c())
for(geneName in master_top_table$fav_gene){
  filteredData <- master_top_table %>% 
    filter(fav_gene == geneName) %>% 
    select(-fav_gene) %>%
    unnest(cols = c(data)) %>% 
    filter(gene %in% all_genes) %>% 
    nest(data=everything())
  newRow <- tibble(fav_gene=geneName, filteredData)
  temp_top_table <- temp_top_table %>% 
    bind_rows(newRow)
}
master_bottom_table <- temp_bottom_table
master_top_table <- temp_top_table

master_positive_filename <- paste0(release, "_master_positive.Rds")
master_positive <- readRDS(file=here::here("data", master_positive_filename)) %>%
  filter(fav_gene %in% all_genes)

master_negative_filename <- paste0(release, "_master_negative.Rds")
master_negative <- readRDS(file=here::here("data", master_negative_filename)) %>%
  filter(fav_gene %in% all_genes)

surprise_genes_filename <- paste0(release, "_surprise_genes.Rds")
surprise_genes <- all_genes %>% head()

censor_genes_filename <- paste0(release, "_censor_genes.Rds")
censor_genes <- readRDS(file=here::here("data", censor_genes_filename)) %>%
  filter(genes %in% all_genes)

subcell_filename <- paste0(release, "_subcell.Rds")
subcell <- readRDS(file=here::here("data", subcell_filename)) %>%
  filter(gene_name %in% all_genes)

#SAVE DATA------
message("Saving pathways.Rds for go ", pathway_go)
saveRDS(pathways, here::here(tests_data_dir, pathways_filename))

message("Saving gene_summary")
saveRDS(gene_summary, here::here(tests_data_dir, gene_summary_filename))

message("Saving achilles")
saveRDS(achilles, here::here(tests_data_dir, achilles_filename))

# message("Saving achilles cor")
# saveRDS(achilles_cor, here::here(tests_data_dir, achilles_cor_filename))

message("Saving achilles cor nest")
saveRDS(achilles_cor_nest, here::here(tests_data_dir, achilles_cor_nest_filename))

message("Saving expression")
saveRDS(expression, here::here(tests_data_dir, expression_filename))

message("Saving expression_meta")
saveRDS(expression_meta, here::here(tests_data_dir, expression_meta_filename))

message("Saving expression_names")
saveRDS(expression_names, here::here(tests_data_dir, expression_names_filename))

message("Saving proteins")
saveRDS(proteins, here::here(tests_data_dir, proteins_filename))

#read data from generate_depmap_stats.R
file_suffixes_to_copy <- c("_sd_threshold.Rds",
                           "_achilles_lower.Rds",
                           "_achilles_upper.Rds",
                           "_mean_virtual_achilles.Rds",
                           "_sd_virtual_achilles.Rds", 
                           "_expression_upper.Rds",
                           "_expression_lower.Rds",
                           "_mean_virtual_expression.Rds",
                           "_sd_virtual_expression.Rds")
message("Saving value Rds files directly")
for (file_suffix in file_suffixes_to_copy)
  file.copy(
    here::here("data", paste0(release, file_suffix)),
    here::here(tests_data_dir, paste0(release, file_suffix)))

message("Saving master_bottom_table")
saveRDS(master_bottom_table, here::here(tests_data_dir, master_bottom_table_filename))

message("Saving master_top_table")
saveRDS(master_top_table, here::here(tests_data_dir, master_top_table_filename))

message("Saving master_positive")
saveRDS(master_positive, here::here(tests_data_dir, master_positive_filename))

message("Saving master_negative")
saveRDS(master_negative, here::here(tests_data_dir, master_negative_filename))

message("Saving surprise_genes")
saveRDS(surprise_genes, here::here(tests_data_dir, surprise_genes_filename))

message("Saving censor_genes")
saveRDS(censor_genes, here::here(tests_data_dir, censor_genes_filename))

message("Saving subcell")
saveRDS(subcell, here::here(tests_data_dir, subcell_filename))
