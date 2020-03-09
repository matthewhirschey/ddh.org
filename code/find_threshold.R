#load libraries
library(tidyverse)
library(janitor)
library(corrr)
library(moderndive)

#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

achilles <- read_csv(achilles_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))

expression <- read_csv(ccle_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))

achilles_long <- achilles %>% 
  pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

#filter achilles to remove no expression dep scores(special sauce)
expression_long <- expression %>% 
  filter(expression$X1 %in% achilles$X1 == TRUE) %>% #matches cells
  gather("gene", "gene_expression", -X1) %>% 
  arrange(desc(gene_expression))

no_expression <- expression_long %>% 
  filter(gene_expression == 0) %>% 
  unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>% 
  pull(match)

achilles_no0 <- achilles_long %>% 
  unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>% 
  filter(match %in% no_expression == FALSE) %>% 
  select(-match) %>%
  spread(gene, dep_score)

achilles_no0_plot <- achilles_no0 %>% 
  summarise_all(list(~sum(is.na(.)))) %>% 
  gather(gene, NAs) %>% 
  arrange(desc(NAs)) %>% 
  mutate(pos = sum(achilles$X1 %in% expression$X1)-NAs)

#skip toomanyNAs setps

#clean Achilles correlation matrix
achilles_cor <- achilles_no0 %>%
  select(-X1) %>% 
  correlate() #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA

#make some long files
achilles_cor_long <- achilles_cor %>% 
  stretch()

#Permutation tests
virtual_achilles <- achilles_cor_long %>% 
  dplyr::filter(!is.na(r)) %>%   
  rep_sample_n(size = 20000, reps = 1000) %>%
  group_by(replicate) %>% 
  summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_achilles <- mean(virtual_achilles$mean)
sd_virtual_achilles <- mean(virtual_achilles$sd)

sd_threshold <- 3

achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

na_cutoff <- achilles_cor_long %>% 
  dplyr::filter(r > achilles_upper) %>%  #| achilles_correlation_raw < achilles_lower) %>% 
  dplyr::group_by(x) %>% 
  dplyr::summarize(count = n()) %>% 
  dplyr::left_join(achilles_no0_plot, by = c("x" = "gene")) %>% 
  dplyr::arrange(pos) %>% 
  dplyr::top_frac(-.05, pos) %>% 
  dplyr::arrange(pos) %>% 
  dplyr::pull(NAs) %>% 
  min(.)

print(na_cutoff)

#save
saveRDS(na_cutoff, file = here::here("data", paste0(release, "_na_cutoff.Rds")))

