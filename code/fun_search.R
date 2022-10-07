# These functions are for searching and retrieving data that will fill the queryString and reactives()

# SEARCH -----
word_starts_with_regex <- function(query_str) {
  regex(paste0('\\b', query_str), ignore_case=TRUE)
}

calc_rank <- function(query_str, value) {
  # query_str must be within value (though case may be different)
  # if both items are the same length it is a perfect match so returns 1
  # otherwise returns a value comparing the lengths of the strings
  return (str_length(query_str) / str_length(value))
}

search_tables <- function(gene_summary, pathways, expression_names, prism_names, hmdb_names, query_str) {
  if (grepl(',', query_str)) {
    custom_list_search_tables(gene_summary, expression_names, prism_names, query_str)
  } else {
    regular_search_tables(gene_summary, pathways, expression_names, prism_names, hmdb_names, query_str)
  }
}

custom_list_search_tables <- function(gene_summary, expression_names, prism_names, query_str) {
  bind_rows(
    custom_gene_list_search_tables(gene_summary, query_str),
    custom_cell_line_list_search_tables(expression_names, query_str),
    custom_compound_list_search_tables(prism_names, query_str)
  )
}

regular_search_tables <- function(gene_summary, pathways, expression_names, prism_names, hmdb_names, query_str, limit_rows=10) {
  gene_data_search_result <- search_gene_data(gene_summary, pathways, query_str, limit_rows)
  cell_line_search_result <- search_cell_line_data(expression_names, query_str, limit_rows)
  drug_search_result <- search_drug_data(prism_names, hmdb_names, query_str, limit_rows)
  bind_rows(gene_data_search_result, cell_line_search_result, drug_search_result) %>%
    arrange(desc(rank))
}

sort_dedup_and_limit <- function(df, limit_rows) {
  df %>%
    arrange(desc(rank)) %>%
    distinct(key, .keep_all = TRUE) %>%
    head(limit_rows)
}

search_cell_line_data <- function(expression_names, query_str, limit_rows) {
  word_starts_with_query_str <- word_starts_with_regex(query_str)
  
  # find cell lines that start with query_str
  cell_line_rows <- expression_names %>%
    filter(str_detect(cell_line, word_starts_with_query_str)) %>%
    mutate(key = cell_line,
           title = cell_line,
           rank = calc_rank(query_str, cell_line)) %>%
    sort_dedup_and_limit(limit_rows)
  
  # group cell lines into generic grouped format
  cell_line_data <- cell_line_rows %>%
    add_column(subtype='cell') %>%
    group_by(key, subtype, rank) %>%
    nest()
  
  # find lineage that start with query_str
  lineage_rows <- expression_names %>%
    filter(str_detect(lineage, word_starts_with_query_str)) %>%
    group_by(lineage) %>%
    group_nest() %>%
    mutate(key = lineage,
           title = lineage,
           rank = calc_rank(query_str, lineage)) %>%
    sort_dedup_and_limit(limit_rows)
  
  # group lineage data into generic grouped format
  lineage_data <- lineage_rows %>%
    add_column(subtype='lineage') %>%
    group_by(key, title, subtype, rank) %>%
    nest()
  
  # find lineage subtype that start with query_str
  lineage_subtype_rows <- expression_names %>%
    filter(str_detect(lineage_subtype, word_starts_with_query_str)) %>%
    group_by(lineage_subtype) %>%
    group_nest() %>%
    mutate(key = lineage_subtype,
           title = lineage_subtype,
           rank = calc_rank(query_str, lineage_subtype)) %>%
    sort_dedup_and_limit(limit_rows)
  
  # group lineage data into generic grouped format
  lineage_subtype_data <- lineage_subtype_rows %>%
    add_column(subtype='lineage_subtype') %>%
    group_by(key, title, subtype, rank) %>%
    nest()
  
  bind_rows(cell_line_data, lineage_data, lineage_subtype_data)
}

search_drug_data <- function(prism_names, hmdb_names, query_str, limit_rows) {
  word_starts_with_query_str <- word_starts_with_regex(query_str)
  
  prism_name_rows <- prism_names %>%
    filter(str_detect(name, word_starts_with_query_str)) %>%
    mutate(
      key = name,
      title = name,
      rank=calc_rank(query_str, name)) %>%
    sort_dedup_and_limit(limit_rows)
  
  prism_name_data <- prism_name_rows %>%
    add_column(subtype='compound') %>%
    group_by(key, subtype, rank) %>%
    nest()
  
  # find moa that start with query_str
  moa_rows <- prism_names %>%
    filter(str_detect(moa, word_starts_with_query_str)) %>%
    group_by(moa) %>%
    group_nest() %>%
    mutate(key = moa,
           title = moa,
           rank = calc_rank(query_str, moa)) %>%
    sort_dedup_and_limit(limit_rows)
  
  # group lineage data into generic grouped format
  moa_data <- moa_rows %>%
    add_column(subtype='moa') %>%
    group_by(key, title, subtype, rank) %>%
    nest()
  
  hmdb_name_rows <- hmdb_names %>%
    filter(str_detect(name, word_starts_with_query_str)) %>%
    mutate(
      key = cid,
      title = name,
      rank=calc_rank(query_str, name))
  
  hmdb_synonyms_rows <- hmdb_names %>%
    filter(str_detect(synonyms, word_starts_with_query_str)) %>%
    filter(!cid %in% hmdb_name_rows$cid) %>% # skip hmdb_names items already found by name
    mutate(
      key = cid,
      title = name,
      rank=calc_rank(query_str, name))
  
  hmdb_data <- bind_rows(hmdb_name_rows, hmdb_synonyms_rows) %>%
    add_column(subtype='metabolite') %>%
    sort_dedup_and_limit(limit_rows) %>%
    group_by(key, subtype, rank) %>%
    nest()
  
  bind_rows(prism_name_data, moa_data, hmdb_data)
}

search_gene_data <- function(gene_summary, pathways, query_str, limit_rows) {
  word_starts_with_query_str <- word_starts_with_regex(query_str)
  
  # find pathway data
  pathways_data_title <- pathways %>%
    filter(str_detect(pathway, word_starts_with_query_str)) %>%
    mutate(rank=calc_rank(query_str, pathway))
  
  # find go data
  pathways_data_go <- pathways %>%
    filter(str_detect(go, word_starts_with_query_str)) %>%
    mutate(rank=calc_rank(query_str, go))
  
  # nest pathways data underneath generic key, title, and subtype columns
  pathways_data <- unique(bind_rows(pathways_data_title, pathways_data_go)) %>%
    mutate(key = go, title = pathway) %>%
    sort_dedup_and_limit(limit_rows) %>%
    add_column(subtype='pathway') %>%
    group_by(key, title, subtype, rank) %>%
    nest()
  
  # find genes most specific
  genes_data_symbol <- gene_summary %>%
    filter(str_detect(approved_symbol, word_starts_with_query_str)) %>%
    mutate(rank=calc_rank(query_str, approved_symbol))
  
  # find genes most likely alternative
  genes_data_aka <- gene_summary %>%
    filter(str_detect(aka, word_starts_with_query_str)) %>%
    mutate(rank=calc_rank(query_str, aka))
  
  # find genes most generic
  genes_data_name <- gene_summary %>%
    filter(str_detect(approved_name, word_starts_with_query_str)) %>%
    mutate(rank=calc_rank(query_str, approved_name))
  
  # nest gene data underneath generic key, title, and subtype columns
  genes_data <- unique(bind_rows(genes_data_symbol, genes_data_aka, genes_data_name)) %>%
    mutate(key = approved_symbol, title = approved_name) %>%
    sort_dedup_and_limit(limit_rows) %>%
    add_column(subtype='gene') %>%
    group_by(key, subtype, rank) %>%
    nest()
  
  bind_rows(genes_data, pathways_data)
}

custom_gene_list_search_tables <- function(gene_summary, query_str) {
  # Search in gene_summary for matching approved_symbols
  # Returns df with columns:
  # - key - query_str passed in
  # - subtype - 'gene_list'
  # - data - gene_summary row for the specified key each gene symbol, when not found only the approved_symbol will be populated
  # create a df containing valid gene summary rows and just the approved_symbol filled in for unknown gene symbols
  query_gene_symbols <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  gene_symbol_with_known <- query_gene_symbols %>%
    map_dfr(query_symbol_in_gene_summary, gene_summary=gene_summary)
  if(any(gene_symbol_with_known$known)) {
    gene_symbol_with_known %>%
      add_column(key=query_str) %>%
      add_column(subtype='gene_list') %>%
      group_by(key, subtype) %>%
      nest()
  } else {
    tibble()
  }
}

query_symbol_in_gene_summary <- function(gene_symbol, gene_summary) {
  # Searches for exact match on approved_symbol in gene_summary, when not found returns row with only approved_symbol filled in
  matches_gene_symbol_ignore_case <- regex(paste0('^', gene_symbol, '$'), ignore_case = TRUE)
  gene_summary_row <- gene_summary %>%
    dplyr::filter(str_detect(approved_symbol, matches_gene_symbol_ignore_case))
  if (nrow(gene_summary_row) > 0) {
    gene_summary_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(approved_symbol=gene_symbol, known=FALSE)
  }
}

custom_cell_line_list_search_tables <- function(expression_names, query_str) {
  cell_line_symbols <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  cell_line_symbols_with_known <- cell_line_symbols %>%
    map_dfr(query_cell_line_in_expression_names, expression_names=expression_names)
  if(any(cell_line_symbols_with_known$known)) {
    cell_line_symbols_with_known %>%
      add_column(key=query_str) %>%
      add_column(subtype='cell_list') %>%
      group_by(key, subtype) %>%
      nest()
  } else {
    tibble()
  }
}

query_cell_line_in_expression_names <- function(cell_line, expression_names) {
  matches_cell_line_ignore_case <- regex(paste0('^', cell_line, '$'), ignore_case = TRUE)
  expression_names_row <- expression_names %>%
    dplyr::filter(str_detect(cell_line, matches_cell_line_ignore_case))
  if (nrow(expression_names_row) > 0) {
    expression_names_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(cell_line=cell_line, known=FALSE)
  }
}

custom_compound_list_search_tables <- function(prism_names, query_str) {
  compound_names <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  compound_names_with_known <- compound_names %>%
    map_dfr(query_compound_in_expression_names, prism_names=prism_names)
  if(any(compound_names_with_known$known)) {
    compound_names_with_known %>%
      add_column(key=query_str) %>%
      add_column(subtype='compound_list') %>%
      group_by(key, subtype) %>%
      nest()
  } else {
    tibble()
  }
}

query_compound_in_expression_names <- function(compound_name, prism_names) {
  matches_compound_name_ignore_case <- regex(paste0('^', compound_name, '$'), ignore_case = TRUE)
  prism_names_row <- prism_names %>%
    dplyr::filter(str_detect(name, matches_compound_name_ignore_case))
  if (nrow(prism_names_row) > 0) {
    prism_names_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(name=compound_name, known=FALSE)
  }
}

# PATHWAY TABLES -----
#this makes a summary table of ALL pathways for browsing
#it populates a conditionalPanel in shiny on the index 'search examples' page
make_pathway_table <- function(table_name = pathways) { 
  pathway_table <- 
    table_name %>% 
    dplyr::mutate(genes = map_chr(table_name$data, function(.x) {
      y <- unlist(.x, use.names = FALSE)
      str_c(y, collapse = ', ')
    })) %>% 
    dplyr::select(pathway, go, genes)
  return(pathway_table)
}
