
library(tidyverse)

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

# Load data
proteins <- readRDS(file = here::here("data", paste0(release, "_proteins.Rds")))

baseurl = "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
url_end = "" # in case you want the query to be more specific

# Get protein IDs
protein_names <- proteins %>% 
  dplyr::pull(uniprot_id) %>% 
  unique()

# API QUERY LOOP --------------------
features_all <- list()
for (i in 1:length(protein_names)) {
  
  url <- paste0(baseurl, protein_names[i], url_end)
  prots_feat <- httr::GET(url, httr::accept_json())
  prots_feat_red <- httr::content(prots_feat)
  
  if(length(prots_feat_red) == 0) {
    next
  }
  
  features <- tribble(~ type, ~ description, ~ category, 
                      ~ begin, ~ end, ~ url)
  
  for (j in 1:length(prots_feat_red[[1]]$features)) {
    
    if (is.null(prots_feat_red[[1]]$features[[j]]$description)) {
      features_temp <- c(prots_feat_red[[1]]$features[[j]]$type, 
                         "NONE", as.character(prots_feat_red[[1]]$features[[j]]$category),
                         as.numeric(prots_feat_red[[1]]$features[[j]]$begin), 
                         as.numeric(prots_feat_red[[1]]$features[[j]]$end))
    } else {
      features_temp <- c(prots_feat_red[[1]]$features[[j]]$type,
                         as.character(prots_feat_red[[1]]$features[[j]]$description),
                         as.character(prots_feat_red[[1]]$features[[j]]$category),
                         as.numeric(prots_feat_red[[1]]$features[[j]]$begin), 
                         as.numeric(prots_feat_red[[1]]$features[[j]]$end))
    }
    if (is.null(prots_feat_red[[1]]$features[[j]]$evidences[[1]]$source$url)) {
      features_temp[length(features)] <- "NONE"
    }
    else {
      features_temp[length(features)] <- prots_feat_red[[1]]$features[[j]]$evidences[[1]]$source$url
    }
    
    names(features_temp) <- c("type", "description", "category", "begin", "end", "url")
    features <- bind_rows(features, features_temp)
  }
  
  features_all[[i]] <- features %>% 
    mutate(accession = protein_names[i],
           entryName = prots_feat_red[[1]]$entryName,
           taxid = prots_feat_red[[1]]$taxid,
           seq_len = nchar(prots_feat_red[[1]]$sequence)) %>% 
    dplyr::select(accession, entryName, everything())
  
  print(paste0("Done ", i, " out of ", length(protein_names), " - Uniprot ID: ", protein_names[i]))
  
}

# Combine data
protein_domains <- bind_rows(features_all) %>% 
  mutate(begin = as.numeric(begin),
         end = as.numeric(end),
         length = end - begin) %>% 
  dplyr::rename(uniprot_id = accession) %>% 
  dplyr::left_join(proteins %>% 
                     dplyr::select(uniprot_id, gene_name),
                   by = "uniprot_id") %>% 
  dplyr::relocate(gene_name, .after = uniprot_id)

# SAVE --------------------
saveRDS(protein_domains, file = here::here("data", paste0(release, "_protein_domains.Rds")))

