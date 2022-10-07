#regenerate from raw (TRUE), or load from saved (FALSE)
regenerate <- FALSE

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

if(regenerate == TRUE) {
#LOAD BIG data 
temp <- tempfile()
download.file(metabolites_url, temp)
hmdb_full <- xmlToDataFrame(unzip(temp)) %>% #extracting XML took a long time for a 4GB XML file!
  janitor::clean_names()
unlink(temp)
unlink("hmdb_metabolites.xml")

#hmdb_full <- xmlToDataFrame(here::here("tmp", "hmdb_metabolites.xml")) %>% janitor::clean_names()

#for getting metabolite classes
class_finder <- function(string) {
  if(stringr::str_detect(string, "class\\b")) {
    class_compounds <- function(string) {
      tmp <- purrr::map(string, ~ unlist(str_split(., pattern = "\\."))) %>% 
        purrr::map(., ~keep(.x, ~ str_detect(.x, "class"))) %>% flatten_chr(.)
      class <- stringr::str_split(tmp, "class of organic compounds known as ")
      return(stringr::str_to_sentence(class[[1]][[2]]))}
    tryCatch(class_compounds(string), 
             error = function(x){"Organic compound"})
  } else if(stringr::str_detect(string, "classified\\b")) {
    classified_compounds <- function(string) {
      tmp <- purrr::map(string, ~ unlist(str_split(., pattern = "\\."))) %>% 
        purrr::map(., ~keep(.x, ~ str_detect(.x, "classified"))) %>% flatten_chr(.)
      class <- stringr::str_split(tmp, "classified as a member of the ")
      return(stringr::str_to_sentence(class[[1]][[2]]))}
    tryCatch(classified_compounds(string), 
             error = function(x){"Organic compound"})
  } else if(stringr::str_detect(string, "steroid")) {
    return("Steroid hormone")
  } else if(stringr::str_detect(string, "nucleoside")) {
    return("Nucleoside")
  } else if(stringr::str_detect(string, "acid")) {
    return("Organic acid")
  } else {
    return("Organic compound")
  }
}

#purrr::map(hmdb_class$description, ~ class_finder(.))

hmdb_class <- 
  hmdb_full %>% 
  #slice(1:10) %>% #for testing
  dplyr::select(name, description) %>% 
  dplyr::mutate(class = purrr::map(description, ~ class_finder(.))) %>% 
  dplyr::select(-description)

hmdb_names <-
  hmdb_full %>%
  #slice(1:10) %>% #for testing
  dplyr::select(name, synonyms, cid = pubchem_compound_id) %>% 
  dplyr::left_join(hmdb_class, by = "name")

hmdb_meta <-
  hmdb_full %>%
  dplyr::select(name, synonyms, cid = pubchem_compound_id, accession, description, chemical_formula, average_molecular_weight, wikipedia_id) %>% 
  dplyr::left_join(hmdb_class, by = "name")

temp <- tempfile()
download.file(metabolite_proteins_url, temp)
hmdb_proteins_raw <- 
  XML::xmlToDataFrame(unzip(temp)) %>% 
  janitor::clean_names()
unlink(temp)
unlink("hmdb_proteins.xml")

censor_proteins <- c("")
censor_metabolites <- c("","Water", "Hydrogen Ion")

hmdb_proteins_long <- 
  hmdb_proteins_raw %>%
  #dplyr::slice(1:20) %>% #for testing
  dplyr::mutate(metabolite_associations = stringr::str_replace_all(metabolite_associations, "HMDB", "\\.HMDB")) %>% 
  tidyr::separate_rows(metabolite_associations, sep = "\\.") %>% 
  dplyr::filter(metabolite_associations != "") %>% 
  tidyr::separate(metabolite_associations, 
                  into = c("metabolite_accession", "metabolite_name"), 
                  sep = 11) %>% 
  dplyr::select(gene_name, metabolite_name, gene_accession = accession, metabolite_accession) %>% 
  dplyr::filter(!gene_name %in% censor_proteins, 
                !metabolite_name %in% censor_metabolites, 
                stringr::str_detect(metabolite_accession, "HMDB")) 

# tmp <- "Water"
# tmp <- "Glucose-6-Phosphate"
# tmp <- "3-Carbamoyl-2-phenylpropionaldehyde"
# tmp <- "-1(11),7,9-trien-11-ol"
# metabolite_string <- "R-95913"


#collapse metabolite name if it's long
collapse_metabolites <- function(metabolite_string) {
  # #skip if missing
  # if(stringr::str_length(metabolite_string) == 0 | is.na(metabolite_string)) {
  #   return(metabolite_string)
  # }
  #skip if too short
  if(stringr::str_length(metabolite_string) < 7) {
    return(metabolite_string)
  }
  #get first word
  new_metabolite_string <- 
    stringr::str_extract(metabolite_string, "[[:alpha:]]\\w+") #alpha omits numbers and punct, w gets word character, + gets one or more
  #if new_string is NA b/c it has a bizarre name (like a drug name), code breaks; so return original string
  if(is.na(new_metabolite_string)) {
    return(metabolite_string)
  }
  #add plus to indicate that it got collapsed
  if(metabolite_string != new_metabolite_string) {
    new_metabolite_string <- glue::glue("{new_metabolite_string} +")}
  return(new_metabolite_string)
}
# collapse_metabolites(tmp)

# simplify_metabolites <- function(df) {
#   # if (nrow(df) < 20) {
#   #   return(df)
#   # }
#   new_df <- 
#     df %>% 
#     mutate(metabolite_name_simple = map_chr(.x = df[[2]], .f = collapse_metabolites))
#   return(new_df)
# }

hmdb_proteins_full <-
  hmdb_proteins_long %>% 
  dplyr::group_by(fav_gene = gene_name) %>% 
  tidyr::nest() %>% 
  dplyr::filter(!fav_gene %in% censor_proteins) %>%
  dplyr::mutate(original_num_rows = map_int(data, nrow)) 

hmdb_proteins <-
  hmdb_proteins_long %>% 
  dplyr::mutate(metabolite_name_simple = map_chr(metabolite_name, collapse_metabolites)) %>% 
  dplyr::distinct(gene_name, metabolite_name_simple, .keep_all = TRUE) %>% 
  dplyr::group_by(fav_gene = gene_name) %>% 
  tidyr::nest() %>% 
  dplyr::filter(!fav_gene %in% censor_proteins) %>% 
  dplyr::mutate(num_rows = map_int(data, nrow)) %>% 
  dplyr::left_join(hmdb_proteins_full, by = "fav_gene", suffix = c("_collapsed", "_original")) %>% 
  dplyr::arrange(desc(num_rows))

#colnames(hmdb_proteins) are fav_gene, data_collapsed, num_rows, data_original, original_num_rows

hmdb_metabolites <-
  hmdb_proteins_long %>% 
  #dplyr::mutate(metabolite_name = stringr::str_to_lower(metabolite_name)) %>% 
  dplyr::group_by(fav_metabolite = metabolite_name) %>% 
  tidyr::nest() %>% 
  dplyr::filter(!fav_metabolite %in% censor_metabolites) %>% 
  dplyr::mutate(num_rows = map_int(data, nrow)) %>% 
  dplyr::arrange(desc(num_rows))

} else {
  #get Rds from zenodo for metabolite information
  tmp_hmdb_meta <- tempfile()
  download.file(hmdb_meta_url, tmp_hmdb_meta)
  hmdb_meta <- readRDS(tmp_hmdb_meta)
  
  tmp_hmdb_metabolites <- tempfile()
  download.file(hmdb_metabolites_url, tmp_hmdb_metabolites)
  hmdb_metabolites <- readRDS(tmp_hmdb_metabolites)
  
  tmp_hmdb_names <- tempfile()
  download.file(hmdb_names_url, tmp_hmdb_names)
  hmdb_names <- readRDS(tmp_hmdb_names)
  
  tmp_hmdb_proteins <- tempfile()
  download.file(hmdb_proteins_url, tmp_hmdb_proteins)
  hmdb_proteins <- readRDS(tmp_hmdb_proteins)
}

#save files
#saveRDS(hmdb_full, file = here::here("data", paste0(release, "_hmdb_full.Rds")))
saveRDS(hmdb_names, file = here::here("data", paste0(release, "_hmdb_names.Rds")))
saveRDS(hmdb_meta, file = here::here("data", paste0(release, "_hmdb_meta.Rds")))
saveRDS(hmdb_proteins, file = here::here("data", paste0(release, "_hmdb_proteins.Rds")))
saveRDS(hmdb_metabolites, file = here::here("data", paste0(release, "_hmdb_metabolites.Rds")))

