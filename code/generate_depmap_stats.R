#load libraries
library(tidyverse)
library(here)
library(corrr)
library(moderndive)

#rm(list=ls()) 

#set params
release <- "19Q3"
start_time <- Sys.time()

#LOAD data 
load(file = here::here("data", paste0(release, "_achilles_cor.Rdata")))

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
saveRDS(sd_threshold, file = here::here("data", "sd_threshold.rds"))
saveRDS(achilles_lower, file = here::here("data", "achilles_lower.rds"))
saveRDS(achilles_upper, file = here::here("data", "achilles_upper.rds"))
saveRDS(mean_virtual_achilles, file = here::here("data", "mean_virtual_achilles.rds"))
saveRDS(sd_virtual_achilles, file = here::here("data", "sd_virtual_achilles.rds"))

#how long
end_time <- Sys.time()
