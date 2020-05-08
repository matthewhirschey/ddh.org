library(tidyverse)

make_summary <- function(gene_symbol, type) {
  if(type == "gene"){
    summary <- gene_summary %>% 
      filter(approved_symbol %in% gene_symbol) %>% 
      select(identifier = approved_symbol, 
             name = approved_name, 
             summary = entrez_summary)
    } else if (type == "pathway") {
      go_term <- pathways %>% 
        unnest(data) %>% 
        filter(gene %in% gene_symbol) %>% 
        count(go, sort = TRUE) %>% 
        slice(1) %>% 
        pull(go)
    summary <- pathways %>% 
      filter(go %in% go_term) %>% 
      select(identifier = go, 
             name = pathway, 
             summary = def)
    } else { #custom_gene_list
      #PLACEHOLDER
      summary <- tibble(
        identifier = character(), 
        name = character(), 
        summary = character())
  }
  return(summary)
}

make_summary(gene_symbol = "SDHA", type = "gene")

render_complete_report <- function (file, gene_symbol, type) {
  summary <- make_summary(gene_symbol, type)
  cellanatogram <- make_cellanatogram(achilles, gene_symbol)
  p1 <- make_celldeps(achilles, expression_join, gene_symbol, mean_virtual_achilles)
  p2 <- make_cellbins(achilles, expression_join, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles, expression_join, gene_symbol) %>% head(10)
  target_achilles_top <- make_achilles_table(achilles, expression_join, gene_symbol) %>% tail(10)
  dep_top <- make_top_table(master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_table(master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_table(master_negative, gene_symbol)
  graph_report <- make_graph_report(master_top_table, master_bottom_table, gene_symbol)
  rmarkdown::render("report_depmap_app.Rmd", output_file = file)
}

render_dummy_report <- function (file, gene_symbol, type) {
  summary <- make_summary(gene_symbol, type)
  rmarkdown::render("report_dummy_depmap.Rmd", output_file = file)
}

render_report_to_file <- function(file, gene_symbol, type) {
  if (gene_symbol %in% colnames(achilles)) { #length(gene_symbol) == 1 && 
    src <- normalizePath('report_depmap_app.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_depmap_app.Rmd', overwrite = TRUE)
    out <- render_complete_report(file, gene_symbol, type)
    file.rename(out, file)
  } else {
    src <- normalizePath('report_dummy_depmap.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_dummy_depmap.Rmd', overwrite = TRUE)
    out <- render_dummy_report(file, gene_symbol, type)
    file.rename(out, file)
  }
}
