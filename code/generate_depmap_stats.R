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

#convert cor
class(achilles_cor) <- c("cor_df", "tbl_df", "tbl", "data.frame")

#make some long files
achilles_cor_long <- achilles_cor %>% 
  stretch()

#Permutation tests
virtual_achilles <- achilles_cor_long %>% 
  filter(!is.na(r)) %>%   
  rep_sample_n(size = 20000, reps = 1000) %>%
  group_by(replicate) %>% 
  summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_achilles <- mean(virtual_achilles$mean)
sd_virtual_achilles <- mean(virtual_achilles$sd)

sd_threshold <- 3

achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

#save
saveRDS(sd_threshold, file = here::here("data", "sd_threshold.Rds"))
saveRDS(achilles_lower, file = here::here("data", "achilles_lower.Rds"))
saveRDS(achilles_upper, file = here::here("data", "achilles_upper.Rds"))
saveRDS(mean_virtual_achilles, file = here::here("data", "mean_virtual_achilles.Rds"))
saveRDS(sd_virtual_achilles, file = here::here("data", "sd_virtual_achilles.Rds"))

#how long
time_end_stats <- Sys.time()
