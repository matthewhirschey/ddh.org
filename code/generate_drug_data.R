#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

regenerate = FALSE

#NAMES-----
#This makes the names for drug search
time_begin_data <- Sys.time()

#censor drugs that have gene names #"lta", "nppb", "prima1"
censor_names <- c("LTA", "NPPB", "PRIMA1")
censor_ids <- c("brd_k52914072_001_01_5_2_5_mts004", "brd_k89272762_001_12_7_2_5_hts", "brd_k15318909_001_10_5_2_5_hts")

#get meta file -----
prism_meta <- readr::read_csv(prismmeta_url, col_names = TRUE) %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(clean_drug = janitor::make_clean_names(column_name)) %>% 
  dplyr::distinct(name, .keep_all = TRUE) %>%  #drop rows that have duplicate names
  dplyr::filter(!name %in% censor_names) %>% #censor 3 drugs from meta
  dplyr::distinct(name, .keep_all = TRUE)

#get CIDS -----
if (regenerate == TRUE) {
  # cids consists of "query" column and "cid" column
  cids <- get_cid(prism_meta$name)
  saveRDS(cids, file = here::here("data", paste0(release, "_cids.Rds")))
  
} else {
  # or grab from
  # cids <- readRDS(file=here::here("data", paste0(release, "_cids.Rds")))
  #get Rds from zenodo for backup
  tmp_cids <- tempfile()
  download.file(cids_url, tmp_cids)
  cids <- readRDS(tmp_cids)
}

cids <- 
  cids %>% 
  dplyr::distinct(query, .keep_all = TRUE)

prism_meta <- 
  prism_meta %>% 
  dplyr::left_join(cids, by = c("name" = "query")) %>%  #"smiles" = "CanonicalSMILES"
  dplyr::filter(!is.na(name))

#get some descriptions -----
if (regenerate == TRUE) {
  get_compound_descriptions <- function(cid) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/description/JSON")
    json <- httr::GET(url)
    if(json$status_code != 200) {return("Reeeeeally interesting compound, but not so interesting that I could get a description of it")}
    
    df <- jsonlite::fromJSON(rawToChar(json$content)) 
    compound_desc <- df[["InformationList"]][["Information"]][["Description"]][[2]]
    return(compound_desc)
  }
  
  compound_descriptions <- tibble::tibble(
    cid = character(), 
    description = character()
  )
  
  for (i in prism_meta$cid) {
    j <- 
      get_compound_descriptions(i)
    tmp <- 
      c(cid = i, 
        description = j) #name cols here
    compound_descriptions <-
      compound_descriptions %>% 
      dplyr::bind_rows(tmp)
    Sys.sleep(0.15) #per https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access
    
  }
  
  compound_descriptions <-
    compound_descriptions %>% 
    dplyr::distinct(cid, .keep_all = TRUE) %>% 
    tidyr::drop_na(cid) #otherwise I add WAY too many descriptions
  
} else {
  # or grab from
  #compound_descriptions <- readRDS(file=here::here("data", paste0(release, "_compound_descriptions.Rds")))
  #get Rds from zenodo for backup
  tmp_compound_descriptions <- tempfile()
  download.file(compound_descriptions_url, tmp_compound_descriptions)
  compound_descriptions <- readRDS(tmp_compound_descriptions)
}
saveRDS(compound_descriptions, file = here::here("data", paste0(release, "_compound_descriptions.Rds")))

#join descriptions
prism_meta <-
  prism_meta %>% 
  dplyr::left_join(compound_descriptions, by = "cid")


#get MESH terms -----
# if (regenerate == TRUE) {
#   get_compound_mesh <- function(name) {
#     #name = "A-33903" #error
#     #name = "valsartan" #404
#     #name = "gepefrine" #supplemental
#     #name = "RS-0481" #direct 4th node
#     #name = "prolylleucylglycinamide" #direct 9th node
#     #name = "cis-urocanic acid" #404
#     if(stringr::str_detect(name, "[:space:]")) {
#       url_name <- stringr::str_replace_all(name, " ", "%20")
#     } else {
#       url_name <- name
#     }
#     
#     mesh_url <- paste0("https://www.ncbi.nlm.nih.gov/mesh/?term=", url_name)
#     
#     mesh <- read_html(mesh_url)
#     
#     if(sum(stringr::str_detect(mesh %>% 
#                                html_elements("a") %>% #"a"
#                                html_text(), "Supplementary")) >= 1) {
#       mesh <- 
#         session(mesh_url) %>%
#         session_follow_link("Supplementary") %>%
#         read_html()
#     }        
#     
#     elements <- 
#       mesh %>% 
#       html_elements(css = "p , #maincontent") 
#     
#     text <- 
#       elements %>%  #found prolylleucylglycinamide with child(10)
#       html_text()
#     
#     id <- stringr::str_extract(text, "MeSH Unique ID: [:upper:]{1}[:digit:]{6,8}") 
#     id <- stringr::str_extract(id, "[:upper:]{1}[:digit:]{6,8}") 
#     id <- id[!is.na(id)]
#     id <- unique(id)
#     
#     if (identical(id, character(0))) {
#       id <- NA
#     }
#     return(id)
#   }
#   
#   compound_mesh <- tibble::tibble(
#     name = character(), 
#     meshid = character()
#   )
#   
#   for (i in prism_meta$name) { #a few http errors req'd restarts [4537:4540]
#     j <- 
#       stringr::str_to_sentence(i) %>% 
#       get_compound_mesh()
#     tmp <- 
#       c(name = i, meshid = j) #name cols here
#     compound_mesh <-
#       compound_mesh %>% 
#       dplyr::bind_rows(tmp)
#     Sys.sleep(1) #per https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access
#     
#     #manual cleaning steps
#     compound_mesh <- 
#       compound_mesh %>% distinct(name, .keep_all = TRUE)
#     
#     saveRDS(compound_mesh, file = here::here("data", paste0(release, "_compound_mesh.Rds")))
#   }
# } else {
#   #load data
#   #compound_mesh <- readRDS(file=here::here("data", paste0(release, "_compound_mesh.Rds")))
#   #get Rds from zenodo for backup
#   tmp_compound_mesh <- tempfile()
#   download.file(compound_mesh_url, tmp_compound_mesh)
#   compound_mesh <- readRDS(tmp_compound_mesh)
#   
# }
# 
# #join MESH
# prism_meta <-
#   prism_meta %>% 
#   dplyr::left_join(compound_mesh, by = "name")

prism_meta <- 
  prism_meta %>% 
  dplyr::distinct(name, .keep_all = TRUE)


#make name/join/search df -----
prism_names <- 
  prism_meta %>% 
  dplyr::select("name", "moa", "cid", "clean_drug") #, "clean_drug"...reqd for drug private

#save files
saveRDS(prism_meta, file = here::here("data", paste0(release, "_prism_meta.Rds")))
saveRDS(prism_names, file = here::here("data", paste0(release, "_prism_names.Rds")))

#how long
time_end_data <- Sys.time()
