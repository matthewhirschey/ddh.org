library(tidyverse)

make_top_table <- function(toptable_data = master_top_table, gene_symbol) {
  toptable_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(desc(r2)) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_bottom_table <- function(bottomtable_data = master_bottom_table, gene_symbol) {
  bottomtable_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(r2) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_enrichment_table <- function(enrichmenttable_data, gene_symbol) { #master_positive, master_negative
  enrichmenttable_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(achilles_data, expression_data, gene_symbol) { #achilles, expression_join, gene_symbol
  target_achilles <- achilles_data %>%
    dplyr::select(X1, any_of(gene_symbol)) %>%
    dplyr::left_join(expression_data, by = "X1") %>%
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

make_pathway_table <- function(table_name = pathways) { #this makes a summary table of pathways for browsing
  pathway_table <- table_name %>% 
    mutate(genes = map_chr(table_name$data, function(.x) {
      y <- unlist(.x, use.names = FALSE)
      str_c(y, collapse = ', ')
    })) %>% 
    select(pathway, go, genes)
  return(pathway_table)
}

query_symbol_in_gene_summary <- function(gene_symbol, gene_summary) {
  # Searches for exact match on approved_symbol in gene_summary, when not found returns row with only approved_symbol filled in
  matches_gene_symbol_ignore_case <- regex(paste0('^', gene_symbol, '$'), ignore_case = TRUE)
  gene_summary_row <- gene_summary %>%
    filter(str_detect(approved_symbol, matches_gene_symbol_ignore_case))
  if (nrow(gene_summary_row) > 0) {
    gene_summary_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(approved_symbol=gene_symbol, known=FALSE)
  }
}

gene_list_query_results_table <- function(gene_summary, query_str) {
  # Splits query_str by ',' into a list of gene symbols.
  # Search in gene_summary for matching approved_symbols
  # Returns df with columns:
  # - key - query_str passed in
  # - contents - 'gene_list'
  # - data - gene_summary row for the specified key each gene symbol, when not found only the approved_symbol will be populated
  query_gene_symbols <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  
  # create a df containing valid gene summary rows and just the approved_symbol filled in for unknown gene symbols
  query_gene_symbols %>%
    map_dfr(query_symbol_in_gene_summary, gene_summary=gene_summary) %>%
    add_column(key=query_str) %>%
    add_column(contents='gene_list') %>%
    group_by(key, contents) %>%
    nest()
}

gene_or_pathway_query_results_table <- function(gene_summary, pathways, query_str, limit_pathways=10, limit_genes=10) {
  # Searches for query_str in gene_summary and pathways and limits the number of rows returned from each.
  # Returns df with columns:
  # - key - pathways$go or gene_summary$approved_symbol
  # - contents - 'pathway' or 'gene'
  # - data - variables from a single row of pathways or gene_summary
  find_word_start_regex <- regex(paste0('(^', query_str, ')|(\\s', query_str, ')'), ignore_case=TRUE)
  
  if (str_detect(query_str, ",")) {
    return(make_gene_list_query_results_table(gene_summary, query_str)) 
  }

  # create regex that finds query at the start of the string(^) or after a space(\\s) and ignores case
  find_word_start_regex <- regex(paste0('(^', query_str, ')|(\\s', query_str, ')'), ignore_case=TRUE)

  # find pathway data
  pathways_data_title <- pathways %>%
    filter(str_detect(pathway, find_word_start_regex)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_pathways) %>%
    select(-length)
  
  # find go data
  pathways_data_go <- pathways %>%
    filter(str_detect(go, find_word_start_regex)) %>%
    head(limit_pathways)
  
  # nest pathways data underneath generic key, title, and contents columns
  pathways_data <- unique(bind_rows(pathways_data_title, pathways_data_go)) %>%
    head(limit_pathways) %>%
    mutate(key = go,
           title = pathway) %>%
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
  genes_data <- unique(bind_rows(genes_data_symbol, genes_data_aka, genes_data_name)) %>%
    head(limit_genes) %>%
    mutate(key = approved_symbol) %>%
    add_column(contents='gene') %>%
    group_by(key, contents) %>%
    nest()

  bind_rows(genes_data, pathways_data)
}

make_cellanatogram_table <- function(cellanatogram_data = subcell, gene_symbol) {
  cellanatogram_data %>% 
    filter_all(any_vars(gene_name %in% gene_symbol)) %>% 
    filter(!is.na(type)) %>%
    add_count(main_location) %>% 
    transmute(Gene = gene_name, Reliability = reliability, Location = main_location, Count = as_factor(n)) %>% 
    arrange(desc(Count))
}
