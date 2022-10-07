#rm(list = ls())

#read current release information to set url for download
source(here::here("code", "current_release.R"))

#GO pathways
go_bp <- vroom::vroom(go_bp_url, col_names = FALSE, delim = "\t") %>%  #, col_names = FALSE, .name_repair = "universal") %>% 
  janitor::clean_names() %>% 
  dplyr::select(-x2) %>% 
  tidyr::unite(col = "gene", x3:x15, sep = "\t") %>% 
  tidyr::separate_rows(gene, sep = "\t") %>% 
  dplyr::rename("pathway" = "x1") %>% 
  dplyr::filter(!is.na(gene), 
         gene != "NA", 
         gene != "") %>% 
  dplyr::group_by(pathway)

go_bp_count <- go_bp %>%  
  dplyr::summarize(count = n()) %>% 
  dplyr::ungroup()

pathways <- go_bp %>% 
  tidyr::nest() %>% #gene_list = c(gene)
  dplyr::left_join(go_bp_count, by = "pathway") %>% 
  dplyr::filter(count < 41) %>%  #5103 total; <50 removes 839; <30 removes 1462; <20 removes 2118
  tidyr::separate(col = "pathway", into = c("pathway", "go"), sep = "\\(GO\\:") %>% 
  tidyr::separate(col = "go", into = "go", sep = "\\)", extra = "drop") %>% 
  dplyr::mutate(pathway = str_trim(pathway, side = "right"), 
         pathway = str_to_title(pathway), 
         pathway = str_replace_all(pathway, "Ii", "II"),
         pathway = str_replace_all(pathway, "ii", "II"))

#GO definitions
go_def <- readr::read_delim(go_def_url, delim = "\n", col_names = "X1") %>% 
  tidyr::separate(X1, into = c("X1", "X2"), sep = ":", extra = "merge") %>% 
  dplyr::filter(X1 == "id" | X1 == "name" | X1 == "def") %>% 
  dplyr::mutate(id = dplyr::case_when(
      X1 == "id" ~ stringr::str_extract(.[[2]], "[:digit:]+"),
      TRUE ~ NA_character_
    )
  ) %>% 
  tidyr::fill(id) %>% #populate ids across all observations, so that spread can work
  tidyr::separate(X2, into = "X2", sep = "\\[") %>%  #drop anything after "[", could save it if you want PMIDs, etc.
  dplyr::select(X1, X2, id) %>% 
  dplyr::filter(X1 != "id") %>% 
  tidyr::pivot_wider(names_from = X1, values_from = X2) %>% 
  dplyr::mutate_at(vars(def), unlist) %>% # unlist def column for rocker/tidyverse:3.6.1 compatibility
  dplyr::mutate(def = str_remove_all(def, '\\"'))

#join
pathways <- pathways %>% 
  dplyr::left_join(go_def, by = c("go" = "id"))
#fix empties (only 28!)
pathways <- pathways %>% 
  dplyr:: mutate(def = dplyr::case_when(
    is.na(def) ~ "No pathway definition", 
    TRUE ~ def), 
  def = stringr::str_trim(def, side = "both")) %>% 
  dplyr::select(-name) %>% 
  dplyr::filter(str_detect(pathways$def, "OBSOLETE", negate = TRUE))

saveRDS(pathways, file = here::here("data", paste0(release, "_pathways.Rds")))
