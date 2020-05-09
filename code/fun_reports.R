library(tidyverse)

make_summary <- function(gene_symbol, 
                         type, 
                         gene_summary = gene_summary, 
                         pathways = pathways) {
  if(type == "gene"){
    summary_table <- gene_summary %>% 
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
    summary_table <- pathways %>% 
      filter(go %in% go_term) %>% 
      select(identifier = go, 
             name = pathway, 
             summary = def) %>% 
      mutate(identifier = str_c("GO:", identifier))
    } else { #custom_gene_list
      tmp <- gene_summary %>%
        filter(approved_symbol %in% gene_symbol) %>%
        pull(approved_symbol)
      id_vec <- c("custom gene list")
      name_vec <- str_c(tmp, collapse = ", ")
      summary_vec <- c("User defined gene input list")
      summary_table <- tibble(
        identifier = id_vec,
        name = name_vec,
        summary = summary_vec
      )
  }
  return(summary_table)
}

#make_summary(gene_symbol = "SDHA", type = "gene")
#make_summary(gene_symbol = c("GHRHR", "CDH3", "GHRH", "IGF1", "PHIP", "WNT1", "GH1"), type = "pathway")
#make_summary(gene_symbol = c("SDHA", "SDHB"), type = "custom_gene_list")

render_complete_report <- function (file, 
                                    gene_symbol, 
                                    type,
                                    gene_summary = gene_summary, 
                                    pathways = pathways,
                                    subcell = subcell, 
                                    achilles = achilles, 
                                    expression_join = expression_join, 
                                    mean_virtual_achilles = mean_virtual_achilles, 
                                    master_top_table = master_top_table, 
                                    master_bottom_table = master_bottom_table,
                                    master_positive = master_positive, 
                                    master_negative = master_negative
                                    ) {
  summary <- make_summary(gene_symbol, type)
  cellanatogram <- make_cellanatogram(subcell, gene_symbol)
  cellanatogram_table <- make_cellanatogram_table(subcell, gene_symbol)
  p1 <- make_celldeps(achilles, expression_join, gene_symbol, mean_virtual_achilles)
  p2 <- make_cellbins(achilles, expression_join, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles, expression_join, gene_symbol) %>% head(10)
  target_achilles_top <- make_achilles_table(achilles, expression_join, gene_symbol) %>% tail(10)
  dep_top <- make_top_table(master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_table(master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_table(master_negative, gene_symbol)
  graph_report <- make_graph_report(master_top_table, master_bottom_table, gene_symbol)
  rmarkdown::render(here::here("code", "report_app.Rmd"), output_file = file)
}
#render_complete_report(file = "tmp.pdf", gene_symbol = "SDHA", type = "gene")
#render_complete_report(file = "tmp.pdf", gene_symbol = c("GHRHR", "CDH3", "GHRH", "IGF1", "PHIP", "WNT1", "GH1"), type = "pathway")
#render_complete_report(file = "tmp.pdf", gene_symbol = c("SDHA", "SDHB"), type = "custom_gene_list")

render_dummy_report <- function (file, 
                                 gene_symbol, 
                                 type, 
                                 gene_summary = gene_summary, 
                                 pathways = pathways) {
  summary <- make_summary(gene_symbol, 
                          type, 
                          gene_summary = gene_summary, 
                          pathways = pathways)
  rmarkdown::render("report_dummy_app.Rmd", output_file = file)
}

render_report_to_file <- function(file, 
                                  gene_symbol, 
                                  type = type, 
                                  gene_summary = gene_summary, 
                                  pathways = pathways,
                                  subcell = subcell, 
                                  achilles = achilles, 
                                  expression_join = expression_join, 
                                  mean_virtual_achilles = mean_virtual_achilles, 
                                  master_top_table = master_top_table, 
                                  master_bottom_table = master_bottom_table,
                                  master_positive = master_positive, 
                                  master_negative = master_negative) {
  if (gene_symbol %in% colnames(achilles)) { #length(gene_symbol) == 1 && 
    src <- normalizePath('report_app.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_app.Rmd', overwrite = TRUE)
    out <- render_complete_report(file, 
                                  gene_symbol, 
                                  type = type, 
                                  gene_summary = gene_summary, 
                                  pathways = pathways,
                                  subcell = subcell, 
                                  achilles = achilles, 
                                  expression_join = expression_join, 
                                  mean_virtual_achilles = mean_virtual_achilles, 
                                  master_top_table = master_top_table, 
                                  master_bottom_table = master_bottom_table,
                                  master_positive = master_positive, 
                                  master_negative = master_negative)
    file.rename(out, file)
  } else {
    src <- normalizePath('report_dummy_app.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_dummy_app.Rmd', overwrite = TRUE)
    out <- render_dummy_report(file, 
                               gene_symbol, 
                               type, 
                               gene_summary = gene_summary, 
                               pathways = pathways)
    file.rename(out, file)
  }
}
