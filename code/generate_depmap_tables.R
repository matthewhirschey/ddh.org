## THIS CODE GENERATES master_top AND master_bottom TABLES

#methods
methods = FALSE

#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))

#LOAD data 
gene_summary <- readRDS(file = here::here("data", paste0(release, "_gene_summary.Rds")))
achilles_cor_nest <- readRDS(file = here::here("data", paste0(release, "_achilles_cor_nest.Rds")))
achilles_lower <- readRDS(file = here::here("data", paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here("data", paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_sd_virtual_achilles.Rds")))
pubmed_concept_pairs <- readRDS(file = here::here("data", paste0(release, "_pubmed_concept_pairs.Rds")))
# gls_regression_pvals <- readRDS(file = here::here("data", paste0(release, "_GLS_p_long.Rds")))

#setup containers
master_top_table <- tibble::tibble(
  fav_gene = character(), 
  data = list()
)
master_bottom_table <- tibble::tibble(
  fav_gene = character(), 
  data = list()
)

#table funs
get_concept_table <- function(concept_table = pubmed_concept_pairs, query_gene) {
  concept_tmp <- 
    concept_table %>% 
    dplyr::filter(target_gene %in% query_gene) %>% 
    tidyr::unnest(nested)
  return(concept_tmp)
}

make_dep_table <- function(dep_table = achilles_cor_nest, 
                           summary_table = gene_summary, 
                           # gls_p = gls_regression_pvals,
                           query_gene, 
                           upper = achilles_upper, 
                           lower = achilles_lower, 
                           top = TRUE) {
  dep <- 
    dep_table %>% 
    dplyr::filter(fav_gene %in% query_gene) %>%
    tidyr::unnest(data) %>% 
    dplyr::ungroup() %>% 
    # dplyr::inner_join(gls_p, by = c("fav_gene" = "gene1", "gene" = "gene2")) %>%
    dplyr::select(-fav_gene)
  
  if (top){
    dep_cor <-
      dep %>%
      filter(r2 > upper) # mean + 3sd
    
    # dep_gls <-
    #   dep %>%
    #   filter(pval < 0.05 & r2 > 0) # GLS top
    
    # dep <- bind_rows(dep_cor, dep_gls)
    dep <- dep_cor
  }
  else {
    dep_cor <-
      dep %>%
      filter(r2 < lower) # mean - 3sd
    
    # dep_gls <-
    #   dep %>%
    #   filter(pval < 0.05 & r2 < 0) # GLS bottom
    
    # dep <- bind_rows(dep_cor, dep_gls)
    dep <- dep_cor
  }
  
  dep <- 
    dep %>%
    dplyr::arrange(desc(.[[2]])) %>% #use column index
    dplyr::left_join(summary_table, by = c("gene" = "approved_symbol")) %>% 
    dplyr::rename(rowname = 1) %>% 
    dplyr::select(1:5)
  return(dep)
}

#define list
#sample <- achilles_cor_nest %>% dplyr::ungroup() %>% dplyr::slice_sample(n = 100) %>% pull(fav_gene) #comment this out
full <- achilles_cor_nest %>% pull(fav_gene)
gene_group <- full #(~60' on a laptop); change to sample for testing, or methods

#methods
make_graph_group <- function(graph_gene) {
  top10 <- make_dep_table(query_gene = graph_gene, top = TRUE) %>% 
    dplyr::top_n(10, wt = .[[2]]) %>% 
    dplyr::pull(var = 1)
  bottom10 <- make_dep_table(query_gene = graph_gene, top = FALSE) %>% 
    dplyr::top_n(-10, wt = .[[2]]) %>% 
    dplyr::pull(var = 1)
  graph_list <- c(graph_gene, top10, bottom10)
  return(graph_list)
}

if(methods == TRUE) {
  methods_gene_query <- make_graph_group(graph_gene = "TP53")
  methods_gene_group <- methods_gene_query
  for (i in methods_gene_query) {
    methods_list <- make_graph_group(graph_gene = i)
    methods_gene_group <- c(methods_gene_group, methods_list)
  }
  gene_group <- unique(methods_gene_group) #this overwrites gene_group with the methods subset
}


#make master_tables
for (fav_gene in gene_group) {
  concept_tmp <- get_concept_table(query_gene = fav_gene)
  
  message(" Dep tables for ", fav_gene)
  dep_top <- make_dep_table(query_gene = fav_gene, top = TRUE)
  dep_top <- 
    dep_top %>% 
    dplyr::left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    dplyr::rename(gene = rowname, 
           name = approved_name, 
           concept_count = n#,
           # GLSpvalue = pval
           ) %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0), 
           concept_index = round((concept_count/max(concept_count))*100), 0) %>% 
    dplyr::select(gene, name, z_score, r2, concept_count, concept_index) # GLSpvalue
  
  top_table <- dep_top %>% 
    dplyr::mutate(fav_gene = fav_gene) %>% 
    dplyr::group_by(fav_gene) %>% 
    tidyr::nest()
  
  master_top_table <- master_top_table %>% 
    dplyr::bind_rows(top_table)
  
  dep_bottom <- make_dep_table(query_gene = fav_gene, top = FALSE)
  
  dep_bottom <-
    dep_bottom %>% 
    dplyr::left_join(concept_tmp, by = c("rowname" = "target_gene_pair")) %>% 
    dplyr::rename(gene = rowname, 
           name = approved_name, 
           concept_count = n#,
           # GLSpvalue = pval
           ) %>% 
    dplyr::mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_achilles)/sd_virtual_achilles, 1), 
           concept_count = replace_na(concept_count, 0), 
           concept_index = round((concept_count/max(concept_count))*100), 0) %>% 
    dplyr::select(gene, name, z_score, r2, concept_count, concept_index) # GLSpvalue
  
  bottom_table <- dep_bottom %>% 
    dplyr::mutate(fav_gene = fav_gene) %>% 
    dplyr::group_by(fav_gene) %>% 
    tidyr::nest()
  
  master_bottom_table <- master_bottom_table %>% 
    dplyr::bind_rows(bottom_table)
}

#save
saveRDS(master_top_table, file=here::here("data", paste0(release, "_master_top_table.Rds")))
saveRDS(master_bottom_table, file=here::here("data", paste0(release, "_master_bottom_table.Rds")))

#Censor
num_genes <- nrow(master_top_table)

genes <- character(num_genes)
num_sim <- numeric(num_genes)

for (i in seq_along(genes)) {
  genes[i] <- master_top_table$fav_gene[i]
  num_sim[i] <- nrow(master_top_table[[2]][[i]])
}

censor_genes <- tibble::tibble(genes, num_sim)

#save
saveRDS(censor_genes, here::here("data", paste0(release, "_censor_genes.Rds")))
