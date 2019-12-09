library(rmarkdown)
library(here)

##run one time
source(here::here("code", "create_gene_summary.R"))

##run upon each release
source(here::here("code", "generate_depmap_data.R")) #first script to generate ddh correlation matrix
source(here::here("code", "generate_depmap_stats.R")) #second script to generate ddh stats

#after running above two lines, then run these two in either order
source(here::here("code", "generate_depmap_tables.R")) #third script to generate ddh tables
source(here::here("code", "generate_depmap_pathways.R")) #fourth script to generate ddh pathways; needs ||

#then rerun the methods Rmd document to generate source for methods.md
rmarkdown::render(input = here::here("code", "methods.Rmd"))
