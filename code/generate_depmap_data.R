
#read current release information to set parameters for download
source(here::here("code", "current_links.R"))

time_begin_data <- Sys.time()
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 3)

##BROAD
achilles_raw <- vroom::vroom(achilles_url, delim = ",")
colnames(achilles_raw) <- gsub("\\s\\(\\d+\\)", "", colnames(achilles_raw))
colnames(achilles_raw)[1] <- "X1"

#add name cleaning step
# ddh::load_ddh_data("tests/data/", object_name = "gene_summary")
# achilles_raw <- ddh::clean_colnames(achilles_raw)

achilles_long <- achilles_raw %>% 
  tidyr::pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

#EXPRESSION(BROAD)
expression_raw <- vroom::vroom(ccle_url, delim = ",")
colnames(expression_raw) <- gsub("\\s\\(\\d+\\)", "", colnames(expression_raw))
colnames(expression_raw)[1] <- "X1"

#repeat cleaning step for expression
# expression_raw <- clean_colnames(expression_raw)

expression_meta <- readr::read_csv(cclemeta_url, col_names = TRUE) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(X1 = dep_map_id, cell_line = stripped_cell_line_name) %>% 
  dplyr::mutate_at("lineage", function(str) {
    str <- stringr::str_replace_all(str, "\\_", " ")
    str <- stringr::str_to_title(str)
    return(str)
  }) %>% 
  dplyr::mutate_at("lineage_subtype", function(str) {
    str <- stringr::str_replace_all(str, "\\_", " ")
    str <- dplyr::if_else(stringr::str_detect(str, "^[:lower:]"), stringr::str_to_title(str), str)
    return(str)
  })

expression_names <- expression_meta %>%
  dplyr::select(X1, cell_line, lineage, lineage_subtype)

expression_names <- expression_meta %>% 
  select(X1, cell_line, lineage, lineage_subtype)
  
#filter achilles to remove no expression dep scores(special sauce)
expression_long <- expression_raw %>% 
  dplyr::filter(expression_raw$X1 %in% achilles_raw$X1) %>% #matches cells
  tidyr::gather("gene", "gene_expression", -X1) %>% 
  dplyr::arrange(desc(gene_expression))

no_expression <- expression_long %>% 
  dplyr::filter(gene_expression == 0) %>% 
  tidyr::unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>% 
  dplyr::pull(match)

achilles_no0 <- achilles_long %>% 
  tidyr::unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>% 
  dplyr::filter(!match %in% no_expression) %>% 
  dplyr::select(-match) %>%
  tidyr::spread(gene, dep_score)

toomanyNAs <- achilles_no0 %>% 
  dplyr::summarise_all(list(~sum(is.na(.)))) %>% 
  tidyr::gather(gene, NAs) %>% 
  dplyr::arrange(desc(NAs)) %>% 
  dplyr::filter(NAs > na_cutoff) %>% #set in current_release.R
  dplyr::pull(gene)

achilles <- achilles_no0 %>% 
  dplyr::select(-one_of(toomanyNAs)) #check to see if achilles has fewer variables than achilles_raw

#clean Achilles correlation matrix
achilles_cor <- achilles %>%
  dplyr::select(-X1) %>% 
  corrr::correlate() %>% #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA
  dplyr::rename(rowname = 1) #update to Corrr 0.4.3 changes to term by default, so changing it back to rowname here using column index

achilles_cor_nest <- achilles_cor %>% 
  dplyr::rename(fav_gene = 1) %>% 
  tidyr::pivot_longer(cols = -fav_gene, names_to = "gene", values_to = "r2") %>% 
  dplyr::group_by(fav_gene) %>% 
  tidyr::nest()

#get it back out
# achilles_cor_nest %>%
#   #filter(fav_gene == "A2M") %>%
#   unnest(cols = c("data")) %>%
#   ungroup()

#save files
#saveRDS(achilles, file = here::here("data", paste0(release, "_achilles.Rds")))
saveRDS(achilles_long, file = here::here("data", paste0(release, "_achilles_long.Rds")))
saveRDS(expression_long, file = here::here("data", paste0(release, "_expression_long.Rds")))
saveRDS(expression_meta, file = here::here("data", paste0(release, "_expression_meta.Rds")))
saveRDS(expression_names, file = here::here("data", paste0(release, "_expression_names.Rds")))
saveRDS(achilles_cor, file = here::here("data", paste0(release, "_achilles_cor.Rds"))) #required for step 4
saveRDS(achilles_cor_nest, file = here::here("data", paste0(release, "_achilles_cor_nest.Rds")))

#how long
time_end_data <- Sys.time()

