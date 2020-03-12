make_top_table <- function(gene_symbol) {
  master_top_table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-fav_gene) %>%
    dplyr::mutate(r2 = round(r2, 3)) %>% 
    dplyr::arrange(desc(r2)) %>%
    dplyr::rename("Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_bottom_table <- function(gene_symbol) {
  master_bottom_table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-fav_gene) %>%
    dplyr::mutate(r2 = round(r2, 3)) %>% 
    dplyr::arrange(r2) %>%
    dplyr::rename("Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_enrichment_table <- function(table, gene_symbol) { #master_positive, master_negative
  table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(gene_symbol) {
  target_achilles <- achilles %>%
    dplyr::select(X1, gene_symbol) %>%
    dplyr::left_join(expression_join, by = "X1") %>%
    dplyr::rename(dep_score = gene_symbol) %>%
    dplyr::select(cell_line, lineage, dep_score) %>%
    dplyr::mutate(dep_score = round(dep_score, 3)) %>% 
    dplyr::arrange(dep_score) %>%
    dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Dependency Score" = "dep_score")
  return(target_achilles)
}