library(tidyverse)

#SET TEST DATA LOCATION & CLEAR -----
tests_data_dir <- "tests/data"

# remove tests/data direction, to nuke old data
unlink(here::here(tests_data_dir), recursive = TRUE)

# create tests/data directory if it doesn't exist
dir.create(here::here(tests_data_dir), showWarnings = FALSE, recursive = TRUE)

#SET TEST DATA -----
pathway_go <-  "1902965"
no_data_gene <- "WASH7P"
extra_gene <- "BRCA1"
cell_line_subtype <- "Hepatocellular Carcinoma"
moa_target <- "cyclooxygenase inhibitor"

#READ DATA FOR IDS------
message("Generating data")
pathways_filename <- "pathways.Rds"
pathways <- readRDS(here::here("data", pathways_filename))
pathways_sub <- pathways %>% 
  dplyr::filter(go == pathway_go)
pathway_genes <- pathways_sub %>% 
  dplyr::pull(data) %>% #gene_list
  dplyr::first() %>% 
  dplyr::pull(gene)

master_bottom_table_filename <- "master_bottom_table.Rds"
master_bottom_table_orig <- readRDS(file=here::here("data", master_bottom_table_filename))

master_top_table_filename <- "master_top_table.Rds"
master_top_table_orig <- readRDS(file=here::here("data", master_top_table_filename))

#GET TEST DATA IDS------
mbt_genes <- 
  master_bottom_table_orig %>%
  dplyr::filter(fav_gene %in% pathway_genes) %>%
  pull(data) %>%
  map(slice_head, n = 20) %>% 
  map("gene") %>%
  unlist()

mtt_genes <- 
  master_top_table_orig %>%
  dplyr::filter(fav_gene %in% pathway_genes) %>%
  pull(data) %>%
  map(slice_head, n = 20) %>% 
  map("gene") %>%
  unlist()

all_genes <- 
  pathway_genes %>%
  append(mbt_genes) %>%
  append(mtt_genes) %>%
  append(c(extra_gene)) %>% 
  append(c(no_data_gene))

all_genes_and_x1 <- append(c("X1"), all_genes)

#GET TEST DATA DRUG IDS------
all_drugs <- readRDS(file=here::here("data", paste0(release, "_prism_meta.Rds"))) %>% 
  dplyr::filter(moa %in% moa_target) %>% 
  dplyr::pull(name)

all_genes_and_drugs <- c(all_drugs, all_genes)

#GET TEST DATA CELL NAMES------
all_cells <- readRDS(file=here::here("data", paste0(release, "_expression_names.Rds"))) %>% 
  dplyr::filter(lineage_subtype %in% cell_line_subtype) %>% 
  dplyr::pull(cell_line)

#GET TEST DATA METABOLITE NAMES------
all_metabolites <- c("Adenosine triphosphate", "NAD", "Calcium")

#FILTER DATA -----
gene_summary_filename <- paste0(release, "_gene_summary.Rds")
gene_summary <- readRDS(here::here("data", gene_summary_filename)) %>% 
  dplyr::filter(approved_symbol %in% all_genes)

gene_location_filename <- paste0(release, "_gene_location.Rds")
gene_location <- readRDS(here::here("data", gene_location_filename)) %>% 
  dplyr::filter(approved_symbol %in% all_genes)

pubmed_filename <- paste0(release, "_pubmed.Rds")
pubmed <- readRDS(here::here("data", pubmed_filename)) %>% 
  dplyr::filter(name %in% all_genes_and_drugs)

cellosaurus_filename <- paste0(release, "_cellosaurus.Rds")
cellosaurus <- readRDS(here::here("data", cellosaurus_filename)) %>% 
  dplyr::filter(name %in% all_cells)
cellosaurus_key_filename <- paste0(release, "_cellosaurus_key.Rds")
cellosaurus_key <- readRDS(here::here("data", cellosaurus_key_filename))

chromosome_filename <- paste0(release, "_chromosome.Rds")
chromosome <- readRDS(here::here("data", chromosome_filename))

achilles_long_filename <- paste0(release, "_achilles_long.Rds")
achilles_long <- readRDS(file=here::here("data", achilles_long_filename)) %>% 
  dplyr::filter(gene %in% all_genes)

achilles_cor_nest_filename <- paste0(release, "_achilles_cor_nest.Rds")
achilles_cor_nest <- readRDS(file=here::here("data", achilles_cor_nest_filename)) %>% 
  dplyr::filter(fav_gene %in% all_genes)

cell_line_mat_filename <- paste0(release, "_cell_line_mat.Rds")
cell_line_mat <- readRDS(file=here::here("data", cell_line_mat_filename)) %>%
  dplyr::select_at(vars(matches(all_cells)))

cell_line_dep_sim_filename <- paste0(release, "_cell_line_dep_sim.Rds")
cell_line_dep_sim <- readRDS(file=here::here("data", cell_line_dep_sim_filename)) %>%
  dplyr::filter(cell1_name %in% all_cells & cell2_name %in% all_cells)

cell_line_exp_sim_filename <- paste0(release, "_cell_line_exp_sim.Rds")
cell_line_exp_sim <- readRDS(file=here::here("data", cell_line_exp_sim_filename)) %>%
  dplyr::filter(cell1_name %in% all_cells & cell2_name %in% all_cells)

# gls_regression_filename <- paste0(release, "_GLS_p_long.Rds")
# gls_regression <- readRDS(file=here::here("data", gls_regression_filename)) %>% 
  # dplyr::filter(gene1 %in% all_genes | gene2 %in% all_genes)

achilles_cell_line_cor_nest_filename <- paste0(release, "_achilles_cell_line_cor_nest.Rds")
achilles_cell_line_cor_nest <- readRDS(file=here::here("data", achilles_cell_line_cor_nest_filename)) %>% 
  dplyr::filter(fav_cell %in% all_cells)

expression_long_filename <- paste0(release, "_expression_long.Rds")
expression_long <- readRDS(file=here::here("data", expression_long_filename)) %>%
  dplyr::filter(gene %in% all_genes)

expression_meta_filename <- paste0(release, "_expression_meta.Rds")
expression_meta <- readRDS(file=here::here("data", expression_meta_filename)) #no select or dplyr::filter

expression_names_filename <- paste0(release, "_expression_names.Rds")
expression_names <- readRDS(file=here::here("data", expression_names_filename)) #no select or dplyr::filter

proteins_filename <- paste0(release, "_proteins.Rds")
proteins <- readRDS(file=here::here("data", proteins_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

signatures_filename <- paste0(release, "_protein_signatures.Rds")
signatures <- readRDS(file=here::here("data", signatures_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

sequence_cluster_filename <- paste0(release, "_signature_clusters.Rds")
sequence_clusters <- readRDS(file=here::here("data", sequence_cluster_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

protein_domains_filename <- paste0(release, "_protein_domains.Rds")
protein_domains <- readRDS(file=here::here("data", protein_domains_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

protein_cluster_names_filename <- paste0(release, "_protein_cluster_names.Rds")
protein_cluster_names <- readRDS(file=here::here("data", protein_cluster_names_filename))

protein_cluster_enrichment_filename <- paste0(release, "_protein_cluster_enrichment.Rds")
protein_cluster_enrichment <- readRDS(file=here::here("data", protein_cluster_enrichment_filename))

prism_names_filename <- paste0(release, "_prism_names.Rds")
prism_names <- readRDS(file=here::here("data", prism_names_filename)) %>%
  dplyr::filter(name %in% all_drugs)

prism_meta_filename <- paste0(release, "_prism_meta.Rds")
prism_meta <- readRDS(file=here::here("data", prism_meta_filename)) %>%
  dplyr::filter(name %in% all_drugs)

prism_long_filename <- paste0(release, "_prism_long.Rds")
prism_long <- readRDS(file = here::here("data", prism_long_filename)) %>%
  dplyr::filter(name %in% all_drugs)

prism_cor_nest_filename <- paste0(release, "_prism_cor_nest.Rds")
prism_cor_nest <- readRDS(file = here::here("data", prism_cor_nest_filename)) %>%
  dplyr::filter(fav_drug %in% all_drugs)

gene_drugs_cor_table_filename <- paste0(release, "_gene_drugs_cor_table.Rds")
gene_drugs_cor_table <- readRDS(file=here::here("data", gene_drugs_cor_table_filename)) %>%
  dplyr::filter(fav_gene %in% all_genes)

drug_genes_cor_table_filename <- paste0(release, "_drug_genes_cor_table.Rds")
drug_genes_cor_table <- readRDS(file=here::here("data", drug_genes_cor_table_filename)) %>%
  dplyr::filter(fav_drug %in% all_drugs)

gene_drugs_table_filename <- paste0(release, "_gene_drugs_table.Rds")
gene_drugs_table <- readRDS(file=here::here("data", gene_drugs_table_filename)) %>%
  dplyr::filter(fav_gene %in% all_genes)

drug_genes_table_filename <- paste0(release, "_drug_genes_table.Rds")
drug_genes_table <- readRDS(file=here::here("data", drug_genes_table_filename)) %>%
  dplyr::filter(fav_drug %in% all_drugs)

hmdb_names_filename <- paste0(release, "_hmdb_names.Rds")
hmdb_names <- readRDS(file=here::here("data", hmdb_names_filename)) %>%
  dplyr::filter(name %in% all_metabolites)

hmdb_meta_filename <- paste0(release, "_hmdb_meta.Rds")
hmdb_meta <- readRDS(file=here::here("data", hmdb_meta_filename)) %>%
  dplyr::filter(name %in% all_metabolites)

hmdb_proteins_filename <- paste0(release, "_hmdb_proteins.Rds")
hmdb_proteins <- readRDS(file=here::here("data", hmdb_proteins_filename)) %>%
  dplyr::filter(fav_gene %in% all_genes)

hmdb_metabolites_filename <- paste0(release, "_hmdb_metabolites.Rds")
hmdb_metabolites <- readRDS(file=here::here("data", hmdb_metabolites_filename)) %>%
  dplyr::filter(fav_metabolite %in% all_metabolites)

cell_metabolites_filename <- paste0(release, "_cell_metabolites.Rds")
cell_metabolites <- readRDS(file=here::here("data", cell_metabolites_filename)) %>%
  dplyr::left_join(expression_names, by = c("DepMap_ID" = "X1")) %>%
  dplyr::filter(cell_line %in% all_cells) %>%
  dplyr::select(-cell_line, -lineage, -lineage_subtype)

common_essentials_filename <- paste0(release, "_common_essentials.Rds")
common_essentials <- readRDS(file=here::here("data", common_essentials_filename)) %>%
  dplyr::filter(gene %in% all_genes)

unique_essential_genes_filename <- paste0(release, "_unique_essential_genes.Rds")
unique_essential_genes <- readRDS(file=here::here("data", unique_essential_genes_filename)) %>%
  dplyr::filter(gene %in% all_genes)

prism_unique_toxic_filename <- paste0(release, "_prism_unique_toxic.Rds")
prism_unique_toxic <- readRDS(file=here::here("data", prism_unique_toxic_filename)) %>%
  dplyr::filter(name %in% all_drugs)

#master achilles tables
master_bottom_table <- master_bottom_table_orig %>%
  dplyr::filter(fav_gene %in% all_genes)

master_top_table <- master_top_table_orig %>%
  dplyr::filter(fav_gene %in% all_genes)

# dplyr::filter the sub-tables of the master tables to only include genes within the all_genes list
temp_bottom_table = tibble::tibble(fav_gene = as.character(), data = list())
for(geneName in master_bottom_table$fav_gene){
  filteredData <- master_bottom_table %>% 
    dplyr::filter(fav_gene == geneName) %>% 
    dplyr::select(-fav_gene) %>%
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::filter(gene %in% all_genes) %>% 
    tidyr::nest(data=everything())
  newRow <- tibble(fav_gene=geneName, filteredData)
  temp_bottom_table <- temp_bottom_table %>% 
    bind_rows(newRow)
}

temp_top_table = tibble(fav_gene = as.character(), data = list())
for(geneName in master_top_table$fav_gene){
  filteredData <- master_top_table %>% 
    dplyr::filter(fav_gene == geneName) %>% 
    dplyr::select(-fav_gene) %>%
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::filter(gene %in% all_genes) %>% 
    tidyr::nest(data=everything())
  newRow <- tibble(fav_gene=geneName, filteredData)
  temp_top_table <- temp_top_table %>% 
    bind_rows(newRow)
}
master_bottom_table <- temp_bottom_table
master_top_table <- temp_top_table

#master cell line tables
master_bottom_table_cell_line_filename <- paste0(release, "_master_bottom_table_cell_line.Rds")
master_bottom_table_cell_line_orig <- readRDS(file=here::here("data", master_bottom_table_cell_line_filename))

master_top_table_cell_line_filename <- paste0(release, "_master_top_table_cell_line.Rds")
master_top_table_cell_line_orig <- readRDS(file=here::here("data", master_top_table_cell_line_filename))

master_bottom_table_cell_line <- master_bottom_table_cell_line_orig %>%
  dplyr::filter(fav_cell %in% all_cells)

master_top_table_cell_line <- master_top_table_cell_line_orig %>%
  dplyr::filter(fav_cell %in% all_cells)

# dplyr::filter the sub-tables of the master tables to only include genes within the all_genes list
temp_bottom_table_cell_line = tibble::tibble(fav_cell = as.character(), data = list())
for(cellName in master_bottom_table_cell_line$fav_cell){
  filteredData <- master_bottom_table_cell_line %>% 
    dplyr::filter(fav_cell == cellName) %>% 
    dplyr::select(-fav_cell) %>%
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::filter(cell %in% all_cells) %>% 
    tidyr::nest(data=everything())
  newRow <- tibble(fav_cell=cellName, filteredData)
  temp_bottom_table_cell_line <- temp_bottom_table_cell_line %>% 
    bind_rows(newRow)
}

temp_top_table_cell_line = tibble::tibble(fav_cell = as.character(), data = list())
for(cellName in master_top_table_cell_line$fav_cell){
  filteredData <- master_top_table_cell_line %>% 
    dplyr::filter(fav_cell == cellName) %>% 
    dplyr::select(-fav_cell) %>%
    tidyr::unnest(cols = c(data)) %>% 
    dplyr::filter(cell %in% all_cells) %>% 
    tidyr::nest(data=everything())
  newRow <- tibble(fav_cell=cellName, filteredData)
  temp_top_table_cell_line <- temp_top_table_cell_line %>% 
    bind_rows(newRow)
}

master_bottom_table_cell_line <- temp_bottom_table_cell_line
master_top_table_cell_line <- temp_top_table_cell_line

master_positive_filename <- paste0(release, "_master_positive.Rds")
master_positive <- readRDS(file=here::here("data", master_positive_filename)) %>%
  dplyr::filter(fav_gene %in% all_genes)

master_negative_filename <- paste0(release, "_master_negative.Rds")
master_negative <- readRDS(file=here::here("data", master_negative_filename)) %>%
  dplyr::filter(fav_gene %in% all_genes)

surprise_genes_filename <- paste0(release, "_surprise_genes.Rds")
surprise_genes <- all_genes %>% head()

censor_genes_filename <- paste0(release, "_censor_genes.Rds")
censor_genes <- readRDS(file=here::here("data", censor_genes_filename)) %>%
  dplyr::filter(genes %in% all_genes)

subcell_filename <- paste0(release, "_subcell.Rds")
subcell <- readRDS(file=here::here("data", subcell_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

male_tissue_filename <- paste0(release, "_male_tissue.Rds")
male_tissue <- readRDS(file=here::here("data", male_tissue_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

female_tissue_filename <- paste0(release, "_female_tissue.Rds")
female_tissue <- readRDS(file=here::here("data", female_tissue_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

tissue_filename <- paste0(release, "_tissue.Rds")
tissue <- readRDS(file=here::here("data", tissue_filename)) %>%
  dplyr::filter(gene_name %in% all_genes)

#COPY SUMMARY PAGE IMAGES-----
#assumes all data are already generated, and will simply copy some. 

# set folders
current_dir <- here::here("data", "images")
new_dir <- here::here(tests_data_dir, "images")

#set vars
#all_genes <- c("BRCA1", "C10orf114") #for testing
#all_genes defined above

#loop through copying
copy_images <- function(type, 
                        vec) {
  # copy the images
  dir.create(here::here(new_dir, type), showWarnings = FALSE, recursive = TRUE)
  for (i in vec) {
    if(file.exists(glue::glue('{current_dir}/{type}/{i}'))){
      file.copy(from = glue::glue('{current_dir}/{type}/{i}'),
                to = here::here(new_dir, type), 
                recursive = TRUE)
      message(glue::glue('copied {i} to test dir'))
    } else {
      #prevent errors if one dir is not found
      message(glue::glue('could not find {i} in data dir'))
    }
  }
}

#copy genes
#test_genes <- c("ROCK1", "ROCK2")
copy_images(type = "gene", 
            vec = all_genes) # test_genes
#copy cells
#test_cells <- c("HEPG2", "HUH7")
copy_images(type = "cell", 
            vec = all_cells)

#copy compounds
#test_compounds <- c("aspirin", "ibuprofen")
copy_images(type = "compound", 
            vec =  all_drugs) #test_compounds)

#SAVE DATA------
message("Saving pathways.Rds for go ", pathway_go)
saveRDS(pathways, here::here(tests_data_dir, pathways_filename))

message("Saving gene_summary")
saveRDS(gene_summary, here::here(tests_data_dir, gene_summary_filename))

message("Saving gene_location")
saveRDS(gene_location, here::here(tests_data_dir, gene_location_filename))

message("Saving pubmed")
saveRDS(pubmed, here::here(tests_data_dir, pubmed_filename))

message("Saving cellosaurus")
saveRDS(cellosaurus, here::here(tests_data_dir, cellosaurus_filename))
saveRDS(cellosaurus_key, here::here(tests_data_dir, cellosaurus_key_filename))

message("Saving chromosome")
saveRDS(chromosome, here::here(tests_data_dir, chromosome_filename))

message("Saving achilles long")
saveRDS(achilles_long, here::here(tests_data_dir, achilles_long_filename))

message("Saving achilles cor nest")
saveRDS(achilles_cor_nest, here::here(tests_data_dir, achilles_cor_nest_filename))

message("Saving achilles cell line cor nest")
saveRDS(achilles_cell_line_cor_nest, here::here(tests_data_dir, achilles_cell_line_cor_nest_filename))

message("Saving cell line similarity table")
saveRDS(cell_line_mat, here::here(tests_data_dir, cell_line_mat_filename))
saveRDS(cell_line_dep_sim, here::here(tests_data_dir, cell_line_dep_sim_filename))
saveRDS(cell_line_exp_sim, here::here(tests_data_dir, cell_line_exp_sim_filename))

# message("Saving achilles GLS pvalues")
# saveRDS(gls_regression, here::here(tests_data_dir, gls_regression_filename))

message("Saving expression_meta")
saveRDS(expression_meta, here::here(tests_data_dir, expression_meta_filename))

message("Saving expression_names")
saveRDS(expression_names, here::here(tests_data_dir, expression_names_filename))

message("Saving proteins")
saveRDS(proteins, here::here(tests_data_dir, proteins_filename))

message("Saving sequence clusters")
saveRDS(signatures, here::here(tests_data_dir, signatures_filename))
saveRDS(sequence_clusters, here::here(tests_data_dir, sequence_cluster_filename))
saveRDS(protein_cluster_names, here::here(tests_data_dir, protein_cluster_names_filename))
saveRDS(protein_cluster_enrichment, here::here(tests_data_dir, protein_cluster_enrichment_filename))
saveRDS(protein_domains, here::here(tests_data_dir, protein_domains_filename))

message("Saving updated expression_long")
saveRDS(expression_long, here::here(tests_data_dir, expression_long_filename))

#PRISM
message("Saving prism names")
saveRDS(prism_names, here::here(tests_data_dir, prism_names_filename))

message("Saving prism meta")
saveRDS(prism_meta, here::here(tests_data_dir, prism_meta_filename))

message("Saving prism long")
saveRDS(prism_long, here::here(tests_data_dir, prism_long_filename))

message("Saving prism cor nest")
saveRDS(prism_cor_nest, here::here(tests_data_dir, prism_cor_nest_filename))

message("Saving drug genes table")
saveRDS(drug_genes_table, here::here(tests_data_dir, drug_genes_table_filename))

message("Saving gene drugs table")
saveRDS(gene_drugs_table, here::here(tests_data_dir, gene_drugs_table_filename))

message("Saving drug genes COR table")
saveRDS(drug_genes_cor_table, here::here(tests_data_dir, drug_genes_cor_table_filename))

message("Saving gene drugs COR table")
saveRDS(gene_drugs_cor_table, here::here(tests_data_dir, gene_drugs_cor_table_filename))

message("Saving metabolite names")
saveRDS(hmdb_names, here::here(tests_data_dir, hmdb_names_filename))

message("Saving metabolite meta")
saveRDS(hmdb_meta, here::here(tests_data_dir, hmdb_meta_filename))

message("Saving metabolite genes table")
saveRDS(hmdb_metabolites, here::here(tests_data_dir, hmdb_metabolites_filename))

message("Saving gene metabolites table")
saveRDS(hmdb_proteins, here::here(tests_data_dir, hmdb_proteins_filename))

message("Saving cell metabolites table")
saveRDS(cell_metabolites, here::here(tests_data_dir, cell_metabolites_filename))

message("Saving common essentials table")
saveRDS(common_essentials, here::here(tests_data_dir, common_essentials_filename))

message("Saving unique essential genes table")
saveRDS(unique_essential_genes, here::here(tests_data_dir, unique_essential_genes_filename))

message("Saving prism unique toxic table")
saveRDS(prism_unique_toxic, here::here(tests_data_dir, prism_unique_toxic_filename))

#read data from stats
file_suffixes_to_copy <- c("_sd_threshold.Rds",
                           "_achilles_lower.Rds",
                           "_achilles_upper.Rds",
                           "_mean_virtual_achilles.Rds",
                           "_sd_virtual_achilles.Rds", 
                           "_sd_threshold_cell.Rds",
                           "_achilles_cell_line_lower.Rds",
                           "_achilles_cell_line_upper.Rds",
                           "_mean_virtual_achilles_cell_line.Rds",
                           "_sd_virtual_achilles_cell_line.Rds",
                           "_gene_expression_upper.Rds",
                           "_gene_expression_lower.Rds",
                           "_mean_virtual_gene_expression.Rds",
                           "_sd_virtual_gene_expression.Rds",
                           "_protein_expression_upper.Rds",
                           "_protein_expression_lower.Rds",
                           "_mean_virtual_protein_expression.Rds",
                           "_sd_virtual_protein_expression.Rds", 
                           "_drug_cor_sd_threshold.Rds",
                           "_prism_cor_lower.Rds",
                           "_prism_cor_upper.Rds",
                           "_mean_virtual_prism_cor.Rds",
                           "_sd_virtual_prism_cor.Rds"
                           )
message("Saving value Rds files directly")
for (file_suffix in file_suffixes_to_copy)
  file.copy(
    here::here("data", paste0(release, file_suffix)),
    here::here(tests_data_dir, paste0(release, file_suffix)))

message("Saving master_bottom_table")
saveRDS(master_bottom_table, here::here(tests_data_dir, master_bottom_table_filename))

message("Saving master_top_table")
saveRDS(master_top_table, here::here(tests_data_dir, master_top_table_filename))

message("Saving master_top_table_cell_line")
saveRDS(master_top_table_cell_line, here::here(tests_data_dir, master_top_table_cell_line_filename))

message("Saving master_bottom_table_cell_line")
saveRDS(master_bottom_table_cell_line, here::here(tests_data_dir, master_bottom_table_cell_line_filename))

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

message("Saving male_tissue")
saveRDS(male_tissue, here::here(tests_data_dir, male_tissue_filename))

message("Saving female_tissue")
saveRDS(female_tissue, here::here(tests_data_dir, female_tissue_filename))

message("Saving tissue")
saveRDS(tissue, here::here(tests_data_dir, tissue_filename))

message("Saving all_genes") #this is for generate images
saveRDS(all_genes, here::here(tests_data_dir, "all_genes.Rds"))
