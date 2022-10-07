#read current release information to set url for download
source(here::here("code", "current_release.R"))

master_top_table <-readRDS(file=here::here("data", paste0(release, "_master_top_table.Rds")))
achilles_cor_nest <- readRDS(file=here::here("data", paste0(release, "_achilles_cor_nest.Rds")))
achilles_upper <- readRDS(file = here::here("data", paste0(release, "_achilles_upper.Rds")))
censor_genes <- readRDS(file = here::here("data", paste0(release, "_censor_genes.Rds")))

#censor threshold at 300
censor_threshold <- 300
top_threshold <- 10

censor_n <- 
  censor_genes %>% 
  filter(num_sim > censor_threshold) %>% 
  pull(genes)

#started with a for-loop on master_top_table, but instead grabbed all the data, took the top 10, filtered for r2; WAY faster
connected_table <- 
  achilles_cor_nest %>% 
  tidyr::unnest(data) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(r2 > achilles_upper, 
                !fav_gene %in% censor_n, 
                !gene %in% censor_n) %>% 
  dplyr::group_by(fav_gene) %>% 
  dplyr::slice_max(r2, n = top_threshold) %>% 
  dplyr::ungroup()
   

#make network
graph_network <- 
  tidygraph::as_tbl_graph(connected_table)

degree_table <-  
  tibble::as_tibble(graph_network) %>%
  dplyr::mutate(degree = igraph::degree(graph_network), 
                degree = dplyr::case_when(
                  name %in% censor_n ~ 0, #censor genes with too many correlations by artifically setting to 0
                  TRUE ~ degree
                )) %>% 
  dplyr::arrange(desc(degree))

#this step fills in dummy data for master_table fav_genes missing from degree_table, which throws an error
missing <- master_top_table %>% 
  dplyr::filter(!fav_gene %in% degree_table$name) %>% 
  dplyr::pull(fav_gene)

missing_table <- tibble::tibble(
  name = c(missing), 
  degree = 0
)

degree_table <-
  degree_table %>% 
  dplyr::bind_rows(missing_table)

#to determine good cutoff for connections > n
decile_filter <- degree_table %>% mutate(decile = ntile(degree, 10)) %>% filter(decile == 10) %>% tail(1) %>% pull(degree)
#degree_table %>% mutate(decile = ntile(degree, 10)) %>% ggplot() + geom_histogram(aes(x = degree, fill = factor(decile)))

#make surprise gene list
find_good_candidate <- function(gene_symbol) { #table = master_top_table
  #this gets the top 10 correlation values
  top_10 <- 
    master_top_table %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(desc(r2)) %>% 
    dplyr::slice(1:10)
  #this looks for 'positive controls' in the top 10...smoking guns
  above <- 
    top_10 %>% 
    dplyr::filter(concept_index > 90) %>% 
    pull(fav_gene)
  #this looks for genes within the top 10 that are under-studied
  below <- 
    top_10 %>% 
    dplyr::filter(concept_index < 10) %>% 
    pull(fav_gene)
  #checks degree table
  connections <-
    degree_table %>% 
    filter(gene_symbol %in% name) %>% 
    pull(degree)
    
  #if both are true, then it's a good candidate for further study
  if(length(above) > 0 && length(below) > 0 && connections > decile_filter){ #connections in top decile; only removes ~100 genes
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#TESTING
#find_good_candidate("SSNA1")
#genes <- c("TP53", "TP53BP1")
#map_lgl(genes, ~ find_good_candidate(.))

#GENERATE SURPRISE
surprise_genes <- master_top_table %>% 
  mutate(good = purrr::map_lgl(fav_gene, ~ find_good_candidate(.))) %>% 
  dplyr::filter(good == TRUE) %>% 
  pull(fav_gene)

saveRDS(degree_table, here::here("data", paste0(release, "_degree_table.Rds")))
saveRDS(surprise_genes, here::here("data", paste0(release, "_surprise_genes.Rds")))
