#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#get proteins.Rds from zenodo for protein information
tmp_proteins <- tempfile()
download.file(proteins_url, tmp_proteins)
proteins <- readRDS(tmp_proteins)

proteins <- proteins %>% 
  group_by(gene_name) %>% 
  mutate(gene_name_uq = case_when(
    n() > 1 ~ paste0(gene_name, " (", gsub("\\_.*", "", entry_name), ")"),
    n() == 1 ~ gene_name)
    ) %>% 
  ungroup()

#save files
saveRDS(proteins, file = here::here("data", paste0(release, "_proteins.Rds")))
