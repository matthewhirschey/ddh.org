
library(org.Hs.eg.db)
library(clusterProfiler)

# Read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

# GENERATE SIGNATURE CLUSTER ENRICHMENT ----------------------------------------
## Load Data
signature_final_clusters <- readRDS(file = here::here("data", paste0(release, "_signature_final_clusters.Rds")))

## Enrichment analysis
### Define gene universe and get ENTREZ IDs
all_genes <- signature_final_clusters$gene_name %>%
  unique() %>% 
  AnnotationDbi::select(org.Hs.eg.db, 
                        keys = .,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")

gene_universe <- all_genes %>% 
  dplyr::pull(ENTREZID)

### Compute ORA
clust_num <- unique(signature_final_clusters$clust)

enriched_clusters <- list()
for (i in 1:length(clust_num)) {
  cluster <- clust_num[[i]]
  
  gene_list <- signature_final_clusters %>%
    dplyr::filter(clust == cluster) %>%
    dplyr::select(uniprot_id, gene_name) %>% 
    dplyr::left_join(all_genes, by = c("gene_name" = "SYMBOL")) %>% 
    dplyr::pull(ENTREZID)
  
  # Gene list size
  minGSSize_ddh <- 10
  maxGSSize_ddh <- 500
  
  if(length(gene_list) > maxGSSize_ddh) {
    maxGSSize_ddh <- length(gene_list)
  }
  
  ## MF
  mf_results <- enrichGO(gene = gene_list,
                         universe = gene_universe,
                         OrgDb = org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE,
                         minGSSize = minGSSize_ddh,
                         maxGSSize = maxGSSize_ddh
                         )
  
  if(!is.null(mf_results)) {
    mf_results_df <- mf_results@result %>%
      dplyr::mutate(ont = "MF")
  } else {
    mf_results_df <- NULL
  }
  
  ## BP
  bp_results <- enrichGO(gene = gene_list,
                         universe = gene_universe,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE,
                         minGSSize = minGSSize_ddh,
                         maxGSSize = maxGSSize_ddh
                         )
  
  if(!is.null(bp_results)) {
    bp_results_df <- bp_results@result %>%
      dplyr::mutate(ont = "BP")
  } else {
    bp_results_df <- NULL
  }
  
  ## CC
  cc_results <- enrichGO(gene = gene_list,
                         universe = gene_universe,
                         OrgDb = org.Hs.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable = TRUE,
                         minGSSize = minGSSize_ddh,
                         maxGSSize = maxGSSize_ddh
                         )
  
  if(!is.null(cc_results)) {
    cc_results_df <- cc_results@result %>%
      dplyr::mutate(ont = "CC")
  } else {
    cc_results_df <- NULL
  }
  
  ## MERGE
  enriched_clusters[[i]] <- dplyr::bind_rows(mf_results_df, bp_results_df, cc_results_df) %>% 
    mutate(cluster = as.factor(cluster))
  
  ## TRACK LOOP
  print(paste0("Done ", i, " out of ", length(clust_num)))
  
}

## MASTER TABLE
enriched_clusters <- enriched_clusters %>% 
  dplyr::bind_rows() %>% 
  remove_rownames()

## SAVE
saveRDS(enriched_clusters, file = here::here("data", paste0(release, "_protein_cluster_enrichment.Rds")))

