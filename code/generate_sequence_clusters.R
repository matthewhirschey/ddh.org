
library(tidyverse)

# GENERATE AA SIGNATURES ----------------------------------------
## Load Data
ddh::load_ddh_data(app_data_dir = app_data_dir, object_name = "proteins")

## Compute AA Signature
data_signature <- proteins %>%
  mutate(sequence_split = str_replace_all(sequence, "(.{1})", "\\1 "), 
         sequence_split = str_split(sequence_split, " ")
         )

list_freqs <- list()
for (i in 1:nrow(data_signature)) { # this will take ~1-2m
  seq <- data_signature$sequence_split[[i]][data_signature$sequence_split[[i]] != ""]
  list_freqs[[i]] <- as.data.frame(100*(table(seq)/length(seq))) %>%
    pivot_wider(names_from = "seq", values_from = Freq)
}

protein_signatures <- bind_rows(list_freqs) %>%
  mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>%
  mutate(uniprot_id = proteins$uniprot_id) %>%
  dplyr::inner_join(proteins, by = "uniprot_id") %>%
  dplyr::select(uniprot_id:ec, A:U)

# SAVE
saveRDS(protein_signatures, file = here::here(app_data_dir, "protein_signatures.Rds"))

# GENERATE AA CLUSTERS ----------------------------------------
## Load Data
ddh::load_ddh_data(app_data_dir = app_data_dir, object_name = "protein_signatures")

data_dr <- protein_signatures %>% # Subset AAs
  dplyr::select(A:U)

## Compute "boring"/control proteins
# signatures_clean <- signatures %>%
#   dplyr::select(uniprot_id, A:U)
# 
# aa_mean <- signatures_clean %>%
#   dplyr::summarise_if(is.numeric, list(mean = mean)) %>%
#   dplyr::rename_all(~ str_remove(., "_mean"))
# 
# aa_sd <- signatures_clean %>%
#   dplyr::summarise_if(is.numeric, list(sd = sd)) %>%
#   dplyr::rename_all(~ str_remove(., "_sd"))
# 
# aa_mean_sd <- bind_rows(aa_mean, aa_sd) %>%
#   t() %>%
#   as_tibble() %>%
#   dplyr::mutate(amino_acid = colnames(aa_mean),
#                 lower = V1 - V2,
#                 upper = V1 + V2) %>%
#   dplyr::select(amino_acid, lower, upper)
# 
# tmp_res <- list(1:nrow(signatures_clean))
# for (i in 1:nrow(signatures_clean)) {
# 
#   tmp_prot <- signatures_clean %>%
#     dplyr::slice(i) %>%
#     dplyr::select(-uniprot_id) %>%
#     t() %>%
#     as_tibble() %>%
#     mutate(lower = V1 < aa_mean_sd$lower,
#            upper = V1 > aa_mean_sd$upper)
# 
#   tmp_res[[i]] <- all(!tmp_prot$lower) & all(!tmp_prot$upper)
#   print(paste0("Loop ", i, " completed!"))
# }
# 
# boring_proteins <- data.frame(uniprot_id = signatures_clean$uniprot_id,
#                               boring = unlist(tmp_res)) # 989 boring proteins
# 
# ## Batch correction
# boring_signatures <- signatures %>%
#   dplyr::select(uniprot_id, A:U) %>%
#   dplyr::filter(uniprot_id %in% (boring_proteins %>%
#                                    dplyr::filter(boring) %>%
#                                    dplyr::pull(uniprot_id))
#   ) %>%
#   dplyr::select(A:U)
# 
# data_dr <- signatures %>% # Subset AAs
#   dplyr::select(A:U)
# 
# ## PCA on "boring"/control
# pca_boring_res <- prcomp(boring_signatures)
# # summary(pca_boring_res)$importance[3,15] # First 15 components explain 95% of total variance
# pca_boring_res_df <- pca_boring_res$rotation[, 1:15]
# 
# ## Substraction
# data_dr <- as.matrix(data_dr) - (as.matrix(data_dr) %*% tcrossprod(pca_boring_res_df, pca_boring_res_df))

## Set Parameters
umap_neighbors <- sqrt(nrow(data_dr)) # https://en.wikipedia.org/wiki/Random_coil

## Dimension Reduction UMAP
res_umap <- uwot::umap(data_dr,
                       n_neighbors = umap_neighbors, 
                       n_components = 10,
                       metric = "cosine",
                       min_dist = 0.01
                       )

res_umap_dr <- data.frame(res_umap)

## Visualization UMAP
res_umap2 <- uwot::umap(res_umap_dr,
                        n_neighbors = 10, 
                        n_components = 2
                        )

res_umap_dr2 <- data.frame(res_umap2)

## HDBSCAN Clustering
signature_clust <- dbscan::hdbscan(res_umap2,
                                   minPts = 40 # 25 == Optimized for the number of clusters ("elbow" method)
                                   )

signature_clusters <- data.frame(uniprot_id = protein_signatures$uniprot_id,
                                 gene_name = protein_signatures$gene_name,
                                 res_umap2,
                                 clust = as.factor(signature_clust$cluster),
                                 member_prob = signature_clust$membership_prob)

# ggplot(signature_clusters, aes(X1, X2, color = clust)) +
#   geom_point() +
#   theme_bw() +
#   theme(legend.position = "none") +
#   scale_colour_viridis_d()
# 
# ggplot(signature_clusters %>% 
#          filter(clust != 0), 
#        aes(X1, X2, color = clust)) +
#   geom_point() +
#   theme_bw() +
#   theme(legend.position = "none") +
#   scale_colour_viridis_d()


## SAVE
# saveRDS(boring_proteins, file = here::here("data", paste0(release, "_boring_proteins.Rds")))
saveRDS(signature_clusters, file = here::here(app_data_dir, "signature_clusters.Rds"))

# GENERATE CLUSTER NAMES ----------------------------------------
## LOAD DATA
text_data <- readRDS(file = here::here("data", paste0(release, "_gene_summary.Rds")))
ddh::load_ddh_data(app_data_dir = app_data_dir, object_name = "signature_clusters")

## PREPARE DATA
prep_data <- text_data %>%
  dplyr::select(approved_symbol, entrez_summary) %>%
  dplyr::rename(description = entrez_summary) %>%
  mutate(description = ifelse(str_detect(description, "[[:alpha:]]"), description, NA)) %>%
  dplyr::filter(!duplicated(approved_symbol)) %>%
  dplyr::filter(!is.na(description))

## PROCESS DATA (this step will take few minutes)
remove_words <- c("gene", "encodes", "kd", "typeand")

proc_text <- prep_data %>%
  mutate(summary = tolower(description), # convert to lower case
         summary = str_replace_all(summary, "\\[.*]", ""), # remove text within square brackets (references)
         summary = str_replace_all(summary, "[^-[:^punct:]]", " "), # remove punctuation marks except -
         summary = str_remove_all(summary, "\\s+(\\d+)\\s"), # remove numbers outside words
         summary = str_remove_all(summary, pattern = paste(paste0("\\b", tidytext::get_stopwords()$word, "\\b"), collapse = "|")), # remove stop words using tidytext
         summary = str_remove_all(summary, pattern = paste(paste0("\\b", lexicon::pos_action_verb, "\\b"), collapse = "|")), # remove verbs
         summary = str_remove_all(summary, pattern = paste(paste0("\\b", remove_words, "\\b"), collapse = "|")), # remove text-specific words
         summary = str_squish(summary), # remove whitespaces from start and end of string and repeated whitespaces inside a string
         n_words = str_count(summary, "\\w+") # count words in each string
  )

# key_words <- proc_text %>%
#   tidytext::unnest_tokens(word, summary) %>%
#   # mutate(word = pluralize::singularize(word)) %>% # make singulars
#   dplyr::filter(nchar(word) != 1) %>% # filter strings of length 1
#   group_by(approved_symbol) %>%
#   count(word) %>%
#   dplyr::arrange(desc(n)) %>%
#   dplyr::slice(1:20) %>% # top 20 words
#   dplyr::mutate(key_words = paste(word, collapse = " ")) %>%
#   dplyr::slice(1) %>%
#   dplyr::select(approved_symbol, key_words)

# proc_data_final <- proc_text %>%
#   left_join(key_words, by = "approved_symbol") %>%
#   dplyr::relocate(summary, .after = description) %>%
#   dplyr::relocate(key_words, .after = summary) %>%
#   dplyr::relocate(n_words, .after = key_words)

## CLUSTER KEYWORDS
remove_words2 <- c("protein", "proteins", "encoded", "family",
                   "region", "member", "involved", "associated",
                   "transcript", "variants", "genes", "transcription",
                   "encoding", "also", "role")

cluster_names <- signature_clusters %>%
  dplyr::left_join(proc_text, by = c("gene_name" = "approved_symbol")) %>%
  dplyr::select(uniprot_id, gene_name, clust, summary) %>%
  tidytext::unnest_tokens(word, summary) %>%
  group_by(clust) %>%
  count(word) %>%
  dplyr::filter(!is.na(word)) %>%
  dplyr::filter(nchar(word) != 1) %>% # filter strings of length 1
  dplyr::filter(!word %in% remove_words2) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::slice(1:20) %>% # top 20 words
  dplyr::mutate(cluster_name = paste(word, collapse = "; ")) %>%
  dplyr::slice(1) %>%
  dplyr::select(clust, cluster_name)

## SAVE
saveRDS(proc_text, file = here::here("data", paste0(release, "_gene_summary_processed.Rds")))
saveRDS(cluster_names, file = here::here("data", paste0(release, "_protein_cluster_names.Rds")))

