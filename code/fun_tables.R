library(tidyverse)

make_top_table <- function(table, gene_symbol) {
  table %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(desc(r2)) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_bottom_table <- function(table, gene_symbol) {
  table %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(r2) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_enrichment_table <- function(table, gene_symbol) { #master_positive, master_negative
  table %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(table, expression_table, gene_symbol) {
  target_achilles <- table %>%
    dplyr::select(X1, any_of(gene_symbol)) %>%
    dplyr::left_join(expression_table, by = "X1") %>%
    dplyr::select(-X1) %>%
    dplyr::select(cell_line, lineage, everything()) %>% 
    dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    dplyr::mutate_at("lineage", function(str) {
      str <- str_replace_all(str, "\\_", " ")
      str <- str_to_title(str)
      return(str)
    }) %>% 
    dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage") %>% 
    dplyr::arrange(.[[3]])
  return(target_achilles)
}

make_query_results_table <- function(gene_summary, pathways, query_str, limit_pathways=10, limit_genes=10) {
  # Searches for query_str in gene_summary and pathways and limits the number of rows returned from each.
  # Returns df with columns:
  # - key - pathways$go or gene_summary$approved_symbol
  # - title - pathways$pathway or gene_summary$approved_name
  # - contents - 'pathway' or 'gene'
  # - data - variables from a single row of pathways or gene_summary

  # create regex that finds query at the start of the string(^) or after a space(\\s) and ignores case
  find_word_start_regex <- regex(paste0('(^', query_str, ')|(\\s', query_str, ')'), ignore_case=TRUE)

  # find pathway data and nest it underneath generic key, title, and contents columns
  pathways_data <- pathways %>%
    filter(str_detect(pathway, find_word_start_regex)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_pathways) %>%
    select(-length) %>% 
    mutate(key = go) %>%
    mutate(title = pathway) %>%
    add_column(contents='pathway') %>%
    group_by(key, title, contents) %>%
    nest()

  # find genes most specific
  genes_data_symbol <- gene_summary %>%
    filter(str_detect(approved_symbol, find_word_start_regex)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_genes) %>% 
    select(-length)

  # find genes most likely alternative
  genes_data_aka <- gene_summary %>%
    filter(str_detect(aka, find_word_start_regex)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_genes) %>% 
    select(-length)

  # find genes most generic
  genes_data_name <- gene_summary %>%
    filter(str_detect(approved_name, find_word_start_regex)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_genes) %>% 
    select(-length)

  # nest gene data underneath generic key, title, and contents columns
  genes_data <- bind_rows(genes_data_symbol, genes_data_aka, genes_data_name) %>%
    head(limit_genes) %>%
    mutate(key = approved_symbol,
           title = approved_name) %>%
    add_column(contents='gene') %>%
    group_by(key, title, contents) %>%
    nest()

  bind_rows(genes_data, pathways_data)
}