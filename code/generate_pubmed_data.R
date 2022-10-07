time_begin_pubmed <- Sys.time()

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))
regenerate <- FALSE

#LOAD DATA -----
gene_summary <- readRDS(file = here::here("data", paste0(release, "_gene_summary.Rds")))
achilles_cor_nest <- readRDS(file = here::here("data", paste0(release, "_achilles_cor_nest.Rds")))
prism_meta <- readRDS(file = here::here("data", paste0(release, "_prism_meta.Rds")))
expression_names <- readRDS(file=here::here("data", paste0(release, "_expression_names.Rds")))

# r <- "rowname" #need to drop "rowname"
# full <- (names(achilles_cor))[!(names(achilles_cor)) %in% r] #f[!f %in% r]
full <-
  achilles_cor_nest %>% 
  pull(fav_gene)

#DOWNLOAD pubtator data, which has all of the gene "concepts" within a paper
bioconcepts2pubtator <- readr::read_tsv(bioconcepts2pubtator_url, 
                                        col_names = c("pmid", "type", "concept_id", "mentions", "resource")) 
gene2pubtator <- 
  bioconcepts2pubtator %>% 
  dplyr::filter(type == "Gene") %>% 
  dplyr::select(pmid, concept_id, mentions) %>% 
  tidyr::separate_rows(concept_id, sep = ";") %>% 
  dplyr::filter(concept_id %in% gene_summary$ncbi_gene_id == TRUE)

gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)

#DOWNLOAD pmids + metadata PMC-ids.csv.gz file via FTP Service, which maps an article’s standard IDs to each other and to other article metadata elements.
pubmed_meta <- 
  read_csv(pubmed_url, #spec(ids_raw) grabs cols that were guessed, so I can manually set the errors
           col_types = cols( #manually set some of these because parser incorrectly guess ~1M
             `Journal Title` = col_character(),
             ISSN = col_character(),
             eISSN = col_character(),
             Year = col_double(),
             Volume = col_character(),
             Issue = col_character(),
             Page = col_character(),
             DOI = col_character(),
             PMCID = col_character(),
             PMID = col_double(),
             `Manuscript Id` = col_character(),
             `Release Date` = col_character()
           )) %>% 
  clean_names() %>% 
  arrange(pmid)

#CO-OCCURENCE -----
#co-occurrence of concepts: count which two genes co-occur in a paper, across all papers
pubmed_concept_pairs <- gene2pubtator %>%
  widyr::pairwise_count(concept_id, pmid, sort = TRUE) %>% 
  left_join(gene_summary, by = c("item1" = "ncbi_gene_id")) %>% 
  left_join(gene_summary, by = c("item2" = "ncbi_gene_id")) %>% 
  transmute(target_gene = approved_symbol.x, target_gene_pair = approved_symbol.y, n)

#nest for faster searching
pubmed_concept_pairs <- pubmed_concept_pairs %>% 
  nest(nested = c(target_gene_pair, n))

#prevent join errors in generate_depmap_tables.R by adding missing genes from names(achilles_cor) to this dataset
missing <- full[which(!full %in% pubmed_concept_pairs$target_gene)]

build_missing_df <- function(missing_gene_vec) {
  missing_dataframe <- tibble(target_gene = character(), 
                              nested = list())
  for (i in missing_gene_vec) {
    mt_df <- tibble(target_gene = as.character(i), 
                    target_gene_pair = as.character(NA), 
                    n = as.numeric(NA))
    
    mt_nest <- mt_df %>% 
      nest(nested = c(target_gene_pair, n))
    
    missing_dataframe <- missing_dataframe %>% 
      bind_rows(mt_nest)
  }
  return(missing_dataframe)
}
missing_df <- build_missing_df(missing)

pubmed_concept_pairs <- pubmed_concept_pairs %>% 
  bind_rows(missing_df)

#PUBMED PLOTS -----
pubmed_gene <- 
  gene2pubtator %>% 
  dplyr::left_join(gene_summary, by = c("concept_id" = "ncbi_gene_id")) %>% #to get gene ids
  dplyr::left_join(pubmed_meta, by = c("pmid" = "pmid")) %>% #to get pub years
  dplyr::select(pmid, name = approved_symbol, year, pmcid)

pubmed_compound <- 
  bioconcepts2pubtator %>% 
  dplyr::filter(type == "Chemical") %>% 
  dplyr::filter(mentions %in% prism_meta$name) %>% 
  #dplyr::mutate(concept_id = stringr::str_replace(pubmed_compound$concept_id, "MESH:", "")) %>% 
  dplyr::left_join(pubmed_meta, by = c("pmid" = "pmid")) %>% #to get pub years
  dplyr::select(pmid, name = mentions, year, pmcid)

pubmed_cell_line <- 
  bioconcepts2pubtator %>% 
  dplyr::filter(type == "CellLine") %>% 
  dplyr::mutate(mentions = str_to_upper(mentions)) %>%
  dplyr::filter(str_detect(mentions, paste(expression_names$cell_line, collapse = "|"))) %>%
  #dplyr::mutate(concept_id = stringr::str_replace(pubmed_compound$concept_id, "MESH:", "")) %>% 
  dplyr::left_join(pubmed_meta, by = c("pmid" = "pmid")) %>% #to get pub years
  dplyr::select(pmid, name = mentions, year, pmcid)

pubmed <-
  pubmed_gene %>% 
  dplyr::bind_rows(pubmed_compound) %>%
  dplyr::bind_rows(pubmed_cell_line)

#fix missing years
if (regenerate == TRUE) {
fetch_year <- function(pmid) {
  record <- RISmed::EUtilsGet(pmid)
  Sys.sleep(0.1)
  year <- RISmed::YearPubmed(record)
  return(year)
}

pubmed_years <- 
  pubmed %>% 
  dplyr::filter(!is.na(year))

pubmed_no_years <- 
  pubmed %>% 
  dplyr::filter(is.na(year))

pubmed_no_years <-
  pubmed_no_years %>% 
  dplyr::mutate(year = furrr::future_map(pmid, fetch_year)) %>%  #error w/ map_dbl: must be a single double, not a double vector of length 3
  unnest(year) #use map and unnest

pubmed <-
  pubmed_years %>% 
  dplyr::bind_rows(pubmed_no_years) %>% 
  dplyr::arrange(pmid)
  
}

#remove dupes
pubmed <- 
  pubmed %>% 
  dplyr::distinct(pmid, .keep_all = TRUE)

pubmed$pmid <- as.character(pubmed$pmid) #make char for ahref in table later

#nest
pubmed <-
  pubmed %>% 
  dplyr::group_by(name) %>% 
  tidyr::nest()

#save files
saveRDS(pubmed_concept_pairs, file = here::here("data", paste0(release, "_pubmed_concept_pairs.Rds")))
saveRDS(pubmed, file = here::here("data", paste0(release, "_pubmed.Rds")))

#how long
time_end_pubmed <- Sys.time()
paste0("this took ", round(time_end_pubmed - time_begin_pubmed, digits = 1)/60, " min")
