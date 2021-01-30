#load libraries
library(tidyverse)
library(here)
library(corrr)
library(moderndive)

#rm(list=ls()) 
time_begin_stats <- Sys.time()

#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

#LOAD data 
achilles_cor <- readRDS(file = here::here("data", paste0(release, "_achilles_cor.Rds")))

#make some long files
achilles_cor_long <- achilles_cor %>% 
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "gene", values_to = "r")

#Permutation tests
virtual_achilles <- achilles_cor_long %>% 
  dplyr::filter(!is.na(r)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_achilles <- mean(virtual_achilles$mean)
sd_virtual_achilles <- mean(virtual_achilles$sd)

sd_threshold <- 3

achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

#save
saveRDS(sd_threshold, file = here::here("data", paste0(release, "_sd_threshold.Rds")))
saveRDS(achilles_lower, file = here::here("data", paste0(release, "_achilles_lower.Rds")))
saveRDS(achilles_upper, file = here::here("data", paste0(release, "_achilles_upper.Rds")))
saveRDS(mean_virtual_achilles, file = here::here("data", paste0(release, "_mean_virtual_achilles.Rds")))
saveRDS(sd_virtual_achilles, file = here::here("data", paste0(release, "_sd_virtual_achilles.Rds")))

#LOAD expression data 
expression <- readRDS(file=here::here("data", paste0(release, "_expression.Rds")))

#make some long files
expression_long <- expression %>% 
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "cell_exp")

#Permutation tests
virtual_expression <- expression_long %>% 
  dplyr::filter(!is.na(cell_exp)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(cell_exp), max = max(cell_exp), min = min(cell_exp), sd = sd(cell_exp))

mean_virtual_expression <- mean(virtual_expression$mean)
sd_virtual_expression <- mean(virtual_expression$sd)

sd_threshold <- 3

expression_upper <- mean_virtual_expression + sd_threshold*sd_virtual_expression
expression_lower <- mean_virtual_expression - sd_threshold*sd_virtual_expression

#save
saveRDS(expression_upper, file = here::here("data", paste0(release, "_expression_upper.Rds")))
saveRDS(expression_lower, file = here::here("data", paste0(release, "_expression_lower.Rds")))
saveRDS(mean_virtual_expression, file = here::here("data", paste0(release, "_mean_virtual_expression.Rds")))
saveRDS(sd_virtual_expression, file = here::here("data", paste0(release, "_sd_virtual_expression.Rds")))

#how long
time_end_stats <- Sys.time()
