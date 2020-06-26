## THIS CODE GENERATES censor_genes

#Load Libraries

library(tidyverse)

#Load Data
master_top_table <- readRDS(file=here::here("data", paste0(release, "_master_top_table.Rds")))

num_genes <- nrow(master_top_table)

genes <- character(num_genes)
num_sim <- numeric(num_genes)

for (i in seq_along(genes)) {
  genes[i] <- master_top_table$fav_gene[i]
  num_sim[i] <- nrow(master_top_table[[2]][[i]])
}

censor_genes <- tibble(genes, num_sim)

#save
saveRDS(censor_genes, here::here("data", paste0(release, "_censor_genes.Rds")))
