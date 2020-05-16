library(rmarkdown)
library(here)

##run one time
source(here::here("code", "create_gene_summary.R"))

##run each time pubmed and pubtator are updated; no sense yet of how often this will happen
source(here::here("code", "update_gene_summary.R"))

##run sequentially upon each release
source(here::here("code", "find_threshold.R")) #run this first to generate the na_cutoff
source(here::here("code", "generate_depmap_data.R")) #script to generate ddh correlation matrix
source(here::here("code", "generate_depmap_stats.R")) #script to generate ddh stats
source(here::here("code", "generate_pubmed_data.R")) #script to generate pubtator relationships
source(here::here("code", "generate_subcell_data.R")) #script to generate subcell data

#after running above sequential scripts, then run these two in either order
source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
source(here::here("code", "generate_depmap_pathways.R")) #fourth script to generate ddh pathways; needs ||

#then rerun the methods Rmd document to generate source for methods.md
rmarkdown::render(input = here::here("code", "methods.Rmd"))
