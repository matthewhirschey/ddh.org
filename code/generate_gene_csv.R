source(here::here("code", "current_release.R"))
source(here::here("code", "app_params.R"))

gene_summary <- readRDS(here::here(app_data_dir, paste0(release, "_gene_summary.Rds")))
df = gene_summary[c("approved_symbol", "uni_prot_id_supplied_by_uni_prot")]
dir.create(here::here(app_cache_dir, "protein_pdbs"), showWarnings = TRUE)
write.csv(df, here::here(app_cache_dir, "protein_pdbs/symbol_and_uniprot.csv"), row.names = FALSE)
