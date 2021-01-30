#this script contains functions for returning text, like summaries

# Gene Summary
summary_gene <- function(summary_table = gene_summary, input = list(), var = "approved_symbol") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$id)) {
    return (NULL)
  }
  gene_summary_var <- summary_table %>%
    dplyr::filter(approved_symbol == input$id) %>%
    dplyr::pull(var) #any column name
  return(gene_summary_var)
}
# Pathway Summary
summary_pathway <- function(summary_table = pathways, input = list(), var = "pathway") {
  if (is.null(input$id)) {
    return (NULL)
  }
  if (var == "data") {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == input$id) %>%
      unnest(data) %>% 
      dplyr::pull(gene) %>% 
      str_c(collapse = ", ")
  } else {
  pathway_summary_var <- summary_table %>%
    dplyr::filter(go == input$id) %>%
    dplyr::pull(var)
  }
  return(pathway_summary_var)
}

# Gene list Summary
summary_gene_list <- function(summary_table = gene_summary, input = list()) {
  if (is.null(input$id)) {
    return (NULL)
  }
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- summary_table %>%
    dplyr::filter(approved_symbol %in% input$gene_symbols) %>%
    dplyr::pull(approved_symbol) %>% 
    str_c(collapse = ", ")
  return(valid_gene_symbols)
}

#summary_gene(input = "SDHA", var = "entrez_summary")

# protein summary
#this is a master function to pull data out of the proteins df
summary_protein <- function(summary_table = proteins, input = list(), var = "gene_name") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$id)) {
    return (NULL)
  }
  protein_summary_var <- summary_table %>%
    dplyr::filter(gene_name %in% input$id) %>%
    dplyr::pull(var) #any column name
  return(protein_summary_var)
}
#summary_protein(input = list(id = "GSS"), var = "gene_name")
#summary_protein(input = list(id = "GSS"), var = "protein_name")
#summary_protein(summary_table = proteins, input = list(id = "GSS"), var = "sequence")

# cell summary
summary_cell <- function(summary_table = expression_names, input = list(), var = "cell_line") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$id)) {
    return (NULL)
  }
  cell_summary_var <- summary_table %>%
    dplyr::filter(cell_line == input$id) %>%
    dplyr::pull(var) #any column name
  return(cell_summary_var)
}

# cell lineage summary
summary_lineage <- function(summary_table = expression_names, input = list(), var = "cell_line") { #default so no error if empty, but this pulls the var out of the df
  if (is.null(input$id)) {
    return (NULL)
  }
  summary_table %>%
    dplyr::filter(lineage == input$id) %>%
    dplyr::pull(var) #any column name
}

summary_compound <- function(summary_table = prism_names, input = list(), var = "name") {
  if (is.null(input$id)) {
    return (NULL)
  }
  summary_table %>%
    dplyr::filter(name == input$id) %>%
    dplyr::pull(var) #any column name
}
