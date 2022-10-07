
library(tidyverse)
library(rvest)
library(RCurl)
library(jsonlite)
library(magrittr)
library(glue)

# Read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

# GENERATE AA SIGNATURES ----------------------------------------
## Load Data
proteins <- readRDS(file = here::here("data", paste0(release, "_proteins.Rds")))

## Pull UniProt IDs
uniprot_list <- proteins %>% 
  pull(uniprot_id)

uniprot_pdb_table_ids <- tibble(uniprot = character(),
                                pdb = list())

num <- 0
for (i in uniprot_list) {
  
  num <- num + 1
  
  from <- "UniProtKB_AC-ID"
  to <- "PDB"
  ids <- i 
  
  fetch_pdb_id <- function(data, id) {
    ## QUERY API
    api_query <- postForm("https://rest.uniprot.org/idmapping/run", 
                          from = from,
                          to = to,
                          ids = id)
    
    api_query <- gsub('\\{\"jobId\":\"', '', noquote(api_query[1]))
    api_query <- as.character(gsub('\\"}', '', api_query))
    
    uniprot_link_results <- glue::glue("https://rest.uniprot.org/idmapping/results/{api_query}")
    
    ## READ JSON FROM URL 
    tmp_pdb <- uniprot_link_results %>%
      jsonlite::read_json()
    
    if(length(tmp_pdb$results) > 0) {
      tmp_pdb <- tmp_pdb$results %>% 
        bind_rows() %>% 
        pull(to) 
    } else {
      tmp_pdb <- NULL
    }
    
    ## SAVE RESULTS
    uniprot_pdb_table_tmp <- tibble(uniprot = id,
                                    pdb = list(tmp_pdb))
    
    data %<>% 
      bind_rows(uniprot_pdb_table_tmp)
    
    return(data)
  }
  
  tryCatch({
    uniprot_pdb_table_ids %<>%
      fetch_pdb_id(id = ids)},
    error = function(y){})
  
  print(paste0("Done ", num, " out of ", length(uniprot_list)))

}

# SAVE
saveRDS(uniprot_pdb_table_ids, file = here::here("data", paste0(release, "_uniprot_pdb_table_ids.Rds")))

# Web scraping to obtain structure titles (https://www.rcsb.org/structure/) -----------------
## Load Data
uniprot_pdb_table_ids <- readRDS(file = here::here("data", paste0(release, "_uniprot_pdb_table_ids.Rds")))

uniprot_pdb_table_annotated <- tibble(uniprot = character(),
                                      data = list())

## Web scraping loop
for (i in 1:nrow(uniprot_pdb_table_ids)) {
  uniprot_pdb_table_tmp <- uniprot_pdb_table_ids %>% 
    dplyr::slice(i) %>% 
    tidyr::unnest(pdb)
  
  if(nrow(uniprot_pdb_table_tmp) > 0) {
    pdb_list <- list()
    
    for (j in 1:nrow(uniprot_pdb_table_tmp)) {
      pdb_tmp <- uniprot_pdb_table_tmp %>% 
        dplyr::slice(j)
      
      fetch_pdb_anno <- function(uniprot, pdb) {
        rcsb_link <- rvest::read_html(glue::glue("https://www.rcsb.org/structure/{pdb}")) %>% 
          rvest::html_elements("body")
        
        title <- rcsb_link %>% 
          rvest::html_element("#structureTitle") %>%
          rvest::html_text2()
        
        doi <- rcsb_link %>% 
          rvest::html_element("#header_doi") %>%
          rvest::html_text2() %>% 
          stringr::word(2, sep = "&nbsp")
        
        organism <- rcsb_link %>% 
          rvest::html_element("#header_organism") %>%
          rvest::html_text2() %>% 
          stringr::word(2, sep = "&nbsp")
        
        expression_system <- rcsb_link %>% 
          rvest::html_element("#header_expression-system") %>%
          rvest::html_text2() %>% 
          stringr::word(2, sep = "&nbsp")
        
        pdb_tmp_anno <- tibble::tibble(uniprot, pdb, title, doi, organism, expression_system)
        
        return(pdb_tmp_anno)
      }
      
      tryCatch({
        pdb_list[[j]] <- fetch_pdb_anno(uniprot = pdb_tmp$uniprot,
                                        pdb = pdb_tmp$pdb)},
        error = function(y){})
    }
    
    # if list returns a NULL element, it is silently dropped in the bind_rows() call
    pdb_list %<>%
      bind_rows() %>% 
      tidyr::nest(data = -uniprot)
    
    uniprot_pdb_table_annotated %<>% 
      dplyr::bind_rows(pdb_list)
  } 
  
  print(paste0("Done ", i, " out of ", nrow(uniprot_pdb_table_ids)))

}

# SAVE
saveRDS(uniprot_pdb_table_annotated, file = here::here("data", paste0(release, "_uniprot_pdb_table.Rds")))

  