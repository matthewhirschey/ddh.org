library(tidyverse)

achilles_long <- readRDS(file=here::here("data/achilles_long.Rds"))

total_cells <- 
  achilles_long %>% 
  distinct(X1) %>% 
  nrow()

essentail_genes_vec <- 
  achilles_long %>% 
  ungroup() %>% 
  filter(dep_score < -1) %>% #is this the right threshold?
  group_by(gene) %>% 
  summarise(num = n()) %>% 
  arrange(desc(num)) %>% 
  filter(num == total_cells) %>% 
  #nrow()
  pull(gene)

#what do about genes like AARS1?
#AARS1 is nearly essentail (less than -1 in 807 cell lines, and the 808th is -0.817)
achilles_long %>% 
  ungroup() %>% 
  filter(gene == "AARS1") %>% 
  arrange(desc(dep_score))


#1. get essentail gene list
#2. get pathways
#viz ideas: 
#network of function using dep scores (using funs in ddh)
#pathway analyzer? 
#mean dep score for these genes (what's the range? tight and low or broad?)

