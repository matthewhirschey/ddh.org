library(rmarkdown)
library(here)

##run one time
source(here::here("code", "create_gene_summary.R"))

##run upon each release
#script to generate ddh correlation matrix
rmarkdown::render(input = here::here("code", "depmap_generate.Rmd"), output_dir = here::here("data"))
#file.remove() could be used to programtically remove the file generated

#script to generate ddh stats
rmarkdown::render(input = here::here("code", "depmap_generate_stats.Rmd"), output_dir = here::here("data"))

#script to generate ddh pathways
rmarkdown::render(input = here::here("code", "depmap_generate_pathways.Rmd"), output_dir = here::here("data"))