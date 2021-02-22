#upon each release, run this file locally

#read current release information 
source(here::here("code", "current_release.R"))

source(here::here("code", "create_gene_summary.R"))
create_gene_summary(gene_names_url, entrez_key, here::here("data", paste0(release, "_", gene_summary_output_filename)), gene_symbol = "TP53")
#copy file (no need to regenerate)
#file.copy(from = here::here("data","20Q3_gene_summary.Rds"), to = here::here("data", paste0(release, "_gene_summary.Rds")), overwrite = TRUE)

source(here::here("code", "find_threshold.R"))
source(here::here("code", "generate_depmap_data.R"))
source(here::here("code", "generate_pubmed_data.R"))
source(here::here("code", "generate_depmap_stats.R"))
source(here::here("code", "generate_subcell_data.R")) #script to generate subcell data

#generate table data
#go to generate_depmap_tables.R
#set methods = TRUE in header, and then source
if (!file.exists(here::here("data", paste0(release, "_master_top_table.Rds")))) {
  stop("you forgot to generate tables")
}

#generate pathway data
source(here::here("code", "generate_depmap_pathways.R"))
master_positive <- generate_positive_data(gene_group = c("TP53"), achilles_cor = achilles_cor, achilles_upper = achilles_upper, gene_summary = gene_summary)
saveRDS(master_positive, file=here::here("data", paste0(release, "_", master_positive_filename)))
master_negative <- generate_negative_data(gene_group = c("TP53"), achilles_cor = achilles_cor, achilles_lower = achilles_lower, gene_summary = gene_summary)
saveRDS(master_negative, file=here::here("data", paste0(release, "_", master_negative_filename)))

#then rerun the methods Rmd document to generate source for methods.md
rmarkdown::render(input = here::here("code", "methods.Rmd"))

