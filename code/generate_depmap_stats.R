#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

#ACHILLES STATS-----
#LOAD data 
achilles_cor_nest <- readRDS(file = here::here("data", paste0(release, "_achilles_cor_nest.Rds")))

#make some long files & Permutation tests
virtual_achilles <- 
  achilles_cor_nest %>%
  tidyr::unnest(cols = c("data")) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(r2)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r2), max = max(r2), min = min(r2), sd = sd(r2))

mean_virtual_achilles <- mean(virtual_achilles$mean)
sd_virtual_achilles <- mean(virtual_achilles$sd)

sd_threshold <- 2

achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

#save
saveRDS(sd_threshold, file = here::here("data", paste0(release, "_sd_threshold.Rds")))
saveRDS(achilles_lower, file = here::here("data", paste0(release, "_achilles_lower.Rds")))
saveRDS(achilles_upper, file = here::here("data", paste0(release, "_achilles_upper.Rds")))
saveRDS(mean_virtual_achilles, file = here::here("data", paste0(release, "_mean_virtual_achilles.Rds")))
saveRDS(sd_virtual_achilles, file = here::here("data", paste0(release, "_sd_virtual_achilles.Rds")))

#GENE EXPRESSION STATS-----
#LOAD expression data 
expression_long <- readRDS(file=here::here("data", paste0(release, "_expression_long.Rds")))

#Permutation tests for gene expression
virtual_expression <- 
  expression_long %>% 
  dplyr::select(X1, gene, gene_expression) %>% 
  dplyr::filter(!is.na(gene_expression)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(gene_expression), 
                   max = max(gene_expression), 
                   min = min(gene_expression), 
                   sd = sd(gene_expression))

mean_virtual_gene_expression <- mean(virtual_expression$mean)
sd_virtual_gene_expression <- mean(virtual_expression$sd)

sd_threshold <- 3

gene_expression_upper <- mean_virtual_gene_expression + sd_threshold*sd_virtual_gene_expression
gene_expression_lower <- mean_virtual_gene_expression - sd_threshold*sd_virtual_gene_expression

#save
saveRDS(gene_expression_upper, file = here::here("data", paste0(release, "_gene_expression_upper.Rds")))
saveRDS(gene_expression_lower, file = here::here("data", paste0(release, "_gene_expression_lower.Rds")))
saveRDS(mean_virtual_gene_expression, file = here::here("data", paste0(release, "_mean_virtual_gene_expression.Rds")))
saveRDS(sd_virtual_gene_expression, file = here::here("data", paste0(release, "_sd_virtual_gene_expression.Rds")))

#PROTEIN EXPRESSION STATS-----
#LOAD expression data 
#expression_long already loaded

#Permutation tests for protein expression
virtual_protein_expression <- expression_long %>% 
  dplyr::select(X1, gene, protein_expression) %>% 
  dplyr::filter(!is.na(protein_expression)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(protein_expression), 
                   max = max(protein_expression), 
                   min = min(protein_expression), 
                   sd = sd(protein_expression))

mean_virtual_protein_expression <- mean(virtual_protein_expression$mean)
sd_virtual_protein_expression <- mean(virtual_protein_expression$sd)

#sd_threshold already set

protein_expression_upper <- mean_virtual_protein_expression + sd_threshold*sd_virtual_protein_expression
protein_expression_lower <- mean_virtual_protein_expression - sd_threshold*sd_virtual_protein_expression

#save
saveRDS(protein_expression_upper, file = here::here("data", paste0(release, "_protein_expression_upper.Rds")))
saveRDS(protein_expression_lower, file = here::here("data", paste0(release, "_protein_expression_lower.Rds")))
saveRDS(mean_virtual_protein_expression, file = here::here("data", paste0(release, "_mean_virtual_protein_expression.Rds")))
saveRDS(sd_virtual_protein_expression, file = here::here("data", paste0(release, "_sd_virtual_protein_expression.Rds")))
