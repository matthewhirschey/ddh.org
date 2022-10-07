source(here::here("code", "current_release.R"))

#  tissue data-----
#import
zip_file <- tempfile(fileext = ".zip")
download.file(tissue_url, zip_file, mode = "wb")
tissue <- 
  readr::read_tsv(zip_file) %>% 
  janitor::clean_names()

#  male tissue data-----
data("hgMale_key")

#cleaning
hgMale_key$value <- NULL
colnames(tissue) <- c("gene", "gene_name", "organ", "value")
tissue$organ <- stringr::str_replace(tissue$organ, " ", "_") 

#join
male_tissue <- tissue
male_tissue <- 
  male_tissue %>% 
  dplyr::full_join(hgMale_key, by = "organ") 

#  female tissue data-----

#import
data("hgFemale_key")

#clean
hgFemale_key$value <- NULL

#join
female_tissue <- tissue
female_tissue <- 
  female_tissue %>% 
  dplyr::full_join(hgFemale_key, by = "organ") 

#save files
saveRDS(female_tissue, file = here::here("data", paste0(release, "_female_tissue.Rds")))
saveRDS(male_tissue, file = here::here("data", paste0(release, "_male_tissue.Rds")))
saveRDS(tissue, file = here::here("data", paste0(release, "_tissue.Rds")))
