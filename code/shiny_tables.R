#LANDING PAGE TABLES----
# Browse Pathways
# module that displays a table of pathways when an link is clicked

browsePathwaysLink <- function (id) {
  ns <- NS(id)
  actionLink(inputId = ns("pathway_click"), h4("Browse the pathways"))
}

browsePathwaysLinkServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$pathway_click, {}) #event to store the 'click'
    }
  )
}

browsePathwaysPanel <- function (id) {
  ns <- NS(id)
  notZeroConditionalPanel(ns("pathway_click"),
                          tags$br(),
                          h4("GO Biological Processes"),
                          DT::dataTableOutput(outputId = ns("pathway_table")))
}

browsePathwaysPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_table <- DT::renderDataTable({
        DT::datatable(make_pathway_table(pathways) %>% 
                        dplyr::mutate(go = map_chr(go, internal_link)) %>% #from fun_helper.R
                        dplyr::rename(Pathway = pathway, GO = go, Genes = genes), 
                      escape = FALSE,
                      options = list(pageLength = 10))
      })      
    }
  )
}

#GENE-----
## GO pathways ---------------------
pathwayList <- function (id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_pathway_list")))),
    fluidRow(DT::dataTableOutput(outputId = ns("pathway_list")))
  )
}

pathwayListServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_pathway_list <- renderText({paste0("GO Biological Processes for ", str_c(data()$content, collapse = ", "))})
      output$pathway_list <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% unique(unlist(pathways$data, use.names = FALSE)), "Not found in any pathways")
        )
        DT::datatable(make_pathway_list(table_name = pathways, input = data()) %>% 
                        dplyr::mutate(go = map_chr(go, internal_link))  %>% #from fun_helper.R
                        dplyr::select(Pathway = pathway, GO = go),
                      escape = FALSE,
                      options = list(paging = FALSE, 
                                     searching = FALSE))
      })      
    }
  )
}

pathwayGeneList <- function (id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_pathway_gene_list")))),
    fluidRow(DT::dataTableOutput(outputId = ns("pathway_gene_list")))
  )
}

pathwayGeneListServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_pathway_gene_list <- renderText({paste0("Genes for GO Biological Processes ", data()$query)})
      output$pathway_gene_list <- DT::renderDataTable({
        shiny::validate(
          need(data()$query %in% pathways$go, "Not found in any pathways"))
        DT::datatable(make_pathway_genes(table_name = pathways, table_join = gene_summary, go_id = data()$query) %>% 
                        dplyr::mutate(gene = map_chr(gene, internal_link))  %>% #from fun_helper.R
                        dplyr::select(Gene = gene, Name = approved_name, AKA = aka),
                      escape = FALSE,
                      options = list(paging = FALSE, 
                                     searching = FALSE))
      })      
    }
  )
}

## Protein Cluster ---------------------
proteinClusterTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cluster_table")))),
    checkboxInput(inputId = ns("show_all_clusters_tab"),
                  "Show all proteins in the cluster",
                  value = FALSE),
    checkboxInput(inputId = ns("show_signatures"),
                  "Show amino acid frequencies",
                  value = FALSE),
    DT::dataTableOutput(outputId = ns("prot_clust_table"))
  )
}

proteinClusterTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cluster_table <- renderText({paste0("Amino Acid Signature Clusters for ", str_c(data()$content, collapse = ", "))})
      output$prot_clust_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        withProgress(message = 'Building a smart clustering table...', {
          DT::datatable(make_clustering_table(input = data(),
                                              cluster = input$show_all_clusters_tab,
                                              show_signature = input$show_signatures) %>%
                          dplyr::mutate(gene_name = map_chr(gene_name, internal_link),
                                        uniprot_id = map_chr(uniprot_id, uniprot_linkr)
                          ) %>%  #from fun_helper.R
                          dplyr::rename('Uniprot ID' = uniprot_id,
                                        'Gene Name' = gene_name, 
                                        'Protein Name' = protein_name, 
                                        'Cluster' = clust,
                                        'Cluster Keywords' = cluster_name),
                        escape = FALSE,
                        options = list(pageLength = 10))
        })
      })
    }
  )
}

## Cluster Enrichment ---------------------
proteinClusterEnrichmentTable <- function(id) {
  ns <- NS(id)
  uiOutput(outputId = ns("conditional_clusterenrichmenttable"))
}

proteinClusterEnrichmentTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_clusterenrichmenttable <- renderUI({
        
        if(!is.null(data()$content)) {
          sig_clust <- make_clustering_table(input = data()) %>%
            dplyr::pull(clust) %>%
            unique()
          
          sig_clust_len <- sig_clust %>% 
            length()
        } else {
          sig_clust <- 0
          sig_clust_len <- 0
        }
        
        if((sig_clust_len == 1) & (sig_clust != 0)) {
        tagList(
          tags$br(),
          fluidRow(actionLink(inputId = session$ns("cluster_tab_enrich_click"), "View cluster enrichment table")),
          tags$br(),
          conditionalPanel(condition = paste0("input['", session$ns("cluster_tab_enrich_click"), "'] != 0"), 
                           fluidRow(h4(textOutput(session$ns("text_cluster_enrichment_table")))),
                           fluidRow(checkboxInput(inputId = session$ns("prot_clust_enrich_table_filter_click"), 
                                                  label = "Filter table", value = FALSE)),
                           fluidRow(radioButtons(inputId = session$ns("pval_clust_enrich"),
                                                 label = "p-value type", 
                                                 choices = c("Raw p-value" = "FDR", # inverse because it's a subtraction 
                                                             "FDR" = "pvalue"), # idem
                                                 inline = TRUE)),
                           selectizeInput(session$ns("ontology_info"),
                                          "Ontology",
                                          choices = c("Biological Process (GO)" = "BP",
                                                      "Molecular Function (GO)" = "MF",
                                                      "Cellular Component (GO)" = "CC"), 
                                          selected = "BP"),
                           DT::dataTableOutput(outputId = session$ns("prot_clust_enrich_table")),
                           tags$br() 
                           )
        )
        } else {
          return(NULL)
        }
      })
      output$text_cluster_enrichment_table <- renderText({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        
        if(length(clust_num) == 1) {
          title_text <- paste0("Enrichment Analysis for Cluster ", 
                               clust_num)
        } else{
          title_text <- NULL
        }
        
        return(title_text)
      })
      output$prot_clust_enrich_table <- DT::renderDataTable({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"),
          need(length(clust_num) == 1, 
               "More than one cluster identified in the table above. Cluster enrichment analysis is only available for individual gene queries or multiple gene queries belonging to the same cluster."))
        withProgress(message = 'Building a smart enrichment table...', {
          DT::datatable(make_clustering_enrichment_table(input = data(),
                                                         ontology = input$ontology_info) %>%
                          dplyr::rename(FDR = p.adjust) %>%
                          dplyr::select(-qvalue, -geneID, -Count) %>%
                          dplyr::mutate(pvalue = signif(pvalue, digits = 3),
                                        FDR = signif(FDR, digits = 3)) %>% 
                          dplyr::select_at(vars(-matches(input$pval_clust_enrich))),
                        filter = if(input$prot_clust_enrich_table_filter_click == FALSE) {'none'} else {'top'},
                        escape = FALSE,
                        options = list(pageLength = 10))
        })
      })
    }
  )
}

##Pubmed----
# module that displays a table for pubmed

pubmedTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_pubmed_table")))),
    tags$br(),
    DT::dataTableOutput(outputId = ns("pubmed_table"))
  )
}

pubmedTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_pubmed_table <- renderText({paste0("Publication history table for ", str_c(data()$content, collapse = ", "))})
      output$pubmed_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% pubmed$name, "No data found for this query"))
        withProgress(message = 'Building a smart table...', {
          DT::datatable(make_pubmed_table(pubmed, input = data()) %>% 
                          dplyr::mutate(pmid = map_chr(pmid, pubmed_linkr, number_only = TRUE) #from fun_helper.R
                          ) %>% 
                          dplyr::mutate(pmcid = map_chr(pmcid, pmc_linkr) #from fun_helper.R
                          ) %>% 
                          dplyr::rename(Name = name, 'Pubmed ID' = pmid, Year = year, PMCID = pmcid),
                        escape = FALSE,
                        options = list(pageLength = 10))
        })
      })
    }
  )
}

##Expression----
# module that displays a table for cell anatogram

cellAnatogramTable <- function(id) {
  ns <- NS(id)
  DT::dataTableOutput(outputId = ns("cellanatogram_table"))
}

cellAnatogramTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogram_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% subcell$gene_name, ""))
        DT::datatable(make_cellanatogram_table(subcell, input = data()),
                      options = list(pageLength = 10))
      })
    }
  )
}

tissueTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(checkboxInput(inputId = ns("tissue_filter_click"), label = "Filter tissue table", value = FALSE)),
    fluidRow(DT::dataTableOutput(outputId = ns("tissueanatogram_table")))
  )
}

tissueTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tissueanatogram_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue Expression Data Found"))
        DT::datatable(make_humananatogram_table(tissue, input = data()),
                      filter = if(input$tissue_filter_click == FALSE) {'none'} else {'top'},
                      options = list(pageLength = 10))
      })
    }
  )
}

##Cell Expression Tables -----
cellGeneExpressionTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_gene_table")))),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_gene_table"))))
}

cellGeneExpressionTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_gene_table <- renderText({paste0("Gene expression table for ", str_c(data()$content, collapse = ", "))})
      output$cell_gene_table <- DT::renderDataTable({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, 
                 "No gene expression data found for this gene.")
          )
          DT::datatable(
            make_expression_table(input = data(), var = "gene") %>%
              dplyr::mutate(`Cell Line` = map_chr(`Cell Line`, cell_linkr, type = "cell")),
            escape = FALSE)
          
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0,
                 "No gene expression data found for this cell line.")
          )
          DT::datatable(
            make_expression_table(input = data(), var = "gene") %>%
              dplyr::mutate(Gene = map_chr(Gene, internal_link)),
            escape = FALSE)
        }
      })
    }
  )
}

cellProteinExpressionTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_protein_table")))),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_protein_table"))))
}

cellProteinExpressionTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_protein_table <- renderText({paste0("Protein expression table for ", str_c(data()$content, collapse = ", "))})
      output$cell_protein_table <- DT::renderDataTable({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter( gene %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
          DT::datatable(
            make_expression_table(input = data(), var = "protein") %>%
              dplyr::mutate(`Cell Line` = map_chr(`Cell Line`, cell_linkr, type = "cell")),
            escape = FALSE)
  
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
          DT::datatable(
            make_expression_table(input = data(), var = "protein") %>%
              dplyr::mutate(Gene = map_chr(Gene, internal_link)),
            escape = FALSE)
        }
      })
    }
  )
}

## Dep Table -----
cellDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_dep_table")))),
    fluidRow(checkboxInput(inputId = ns("dep_filter_click"), label = "Filter dependency table", value = FALSE)),
    fluidRow(DT::dataTableOutput(outputId = ns("target_achilles"))))
}

cellDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data()$content, collapse = ", "))})
      output$target_achilles <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_dep_table(input = data())) != 0, "No data found for this gene."))
        DT::datatable(make_dep_table(input = data()) %>%
                        dplyr::mutate(`Cell Line` = map_chr(`Cell Line`, cell_linkr, type = "cell")) #from fun_helper.R
                      , 
                      filter = if(input$dep_filter_click == FALSE) {'none'} else {'top'}, 
                      escape = FALSE,
                      options = list(pageLength = 10))
      })
    }
  )
}

compoundDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compound_dep_table")))),
    fluidRow(checkboxInput(inputId = ns("compound_dep_filter_click"), label = "Filter dependency table", value = FALSE)),
    fluidRow(DT::dataTableOutput(outputId = ns("compound_dep_table"))))
}

compoundDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_compound_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data()$content, collapse = ", "))})
      output$compound_dep_table <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_dep_table(input = data())) != 0, "No data found for this compound"))
        DT::datatable(make_dep_table(input = data()), 
                      filter = if(input$compound_dep_filter_click == FALSE) {'none'} else {'top'}, 
                      options = list(pageLength = 10))
      })
    }
  )
}
cellLineDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_line_dep_table")))),
    fluidRow(checkboxInput(inputId = ns("dep_cell_line_filter_click"), label = "Filter dependency table", value = FALSE)),
    fluidRow(checkboxGroupInput(inputId = ns("vars_essentials"), 
                                "Select columns:",
                                c("Unique Essential", "Pan Essential"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_line_achilles"))))
}

cellLineDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_line_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data()$content, collapse = ", "))})
      output$cell_line_achilles <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% expression_names$cell_line, "No data found for this cell line"))
        DT::datatable(make_dep_table(input = data()) %>%
                        dplyr::rename("Gene" = "gene", "Name" = "approved_name", "Unique Essential" = "unique_essential", "Pan Essential" = "common_essential") %>%
                        dplyr::mutate(Gene = map_chr(Gene, internal_link)) %>%  #from fun_helper.R
                        dplyr::select("Gene", "Name", contains(data()$content), input$vars_essentials),
                      escape = FALSE, 
                      filter = if(input$dep_cell_line_filter_click == FALSE) {'none'} else {'top'}, 
                      options = list(pageLength = 10))
      })
    }
  )
}

cellLineDrugDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_line_drug_dep_table")))),
    fluidRow(checkboxInput(inputId = ns("dep_cell_line_drug_filter_click"), label = "Filter dependency table", value = FALSE)),
    fluidRow(checkboxGroupInput(inputId = ns("vars_toxic_drugs"), 
                                "Select columns:",
                                c("Uniquely Toxic"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_line_drug_prism"))))
}

cellLineDrugDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_line_drug_dep_table <- renderText({paste0("Drug dependency table generated for ", str_c(data()$content, collapse = ", "))})
      output$cell_line_drug_prism <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% expression_names$cell_line, "No data found for this cell line"))
        DT::datatable(make_dep_table(input = data(), var = "drug") %>% 
                        dplyr::rename("Drug" = "name", "MOA" = "moa", "Uniquely Toxic" = "unique_toxic") %>%
                        dplyr::select("Drug", "MOA", "log2fc", input$vars_toxic_drugs), 
                      filter = if(input$dep_cell_line_drug_filter_click == FALSE) {'none'} else {'top'}, 
                      options = list(pageLength = 10))
      })
    }
  )
}
##Similar----
similarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_top")))),
    fluidRow(
      # shinyWidgets::prettySwitch(inputId = ns("gls_table_top"), 
      #                            "Show GLS table", 
      #                            value = FALSE),
      column(6, checkboxGroupInput(inputId = ns("vars_dep_top"), 
                                   "Select columns:",
                                   c("R^2", "Z-Score", "Co-publication Count", "Co-publication Index"), # "GLS p-value"
                                   selected = c("R^2", "Co-publication Count"), # "GLS p-value"
                                   inline = TRUE)),
      column(6, fluidRow(sliderInput(inputId = ns("num_sim_genes"),
                                     "Remove genes with <n associations:",
                                     min = 100,
                                     max = 1000,
                                     value = 1000,
                                     step = 100)), 
             fluidRow(column(3, actionButton(inputId = ns("censor"), "Submit")),
                      column(3, actionButton(inputId = ns("reset"), "Reset"))))
    ),
    hr(),
    fluidRow(DT::dataTableOutput(outputId = ns("dep_top"))), 
    tags$br()
  )
}

similarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      #censor reactive values
      censor_status <- reactiveValues(choice = FALSE, 
                                      num_sim_genes = 1000
      )
      
      observeEvent(input$censor, {
        censor_status$choice <- TRUE
        censor_status$num_sim_genes <- input$num_sim_genes
        censor_status$num <- nrow(make_top_table(input = data()#, 
                                                 # gls = input$gls_table_top
        ) %>% 
          censor(censor_genes, censor_status$choice, censor_status$num_sim_genes))
      })
      
      observeEvent(input$reset, {
        censor_status$choice <- FALSE
        censor_status$num_sim_genes <- 1000
        updateSliderInput(session, inputId = "num_sim_genes", value = 1000)
        censor_status$num <- nrow(make_top_table(input = data()#, 
                                                 # gls = input$gls_table_top
        )
        )
      })
      observeEvent(input$sim_pathway_click, { #event to store the 'click'
      })
      output$text_dep_top <- renderText({paste0(censor_status$num, " genes with similar dependencies as ", str_c(data()$content, collapse = ", "))})      
      output$dep_top <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% master_top_table$fav_gene, "No data found for this gene."))
        censor_status$num <- nrow(make_top_table(input = data()#, 
                                                 # gls = input$gls_table_top
        ) %>% 
          censor(censor_genes, censor_status$choice, censor_status$num_sim_genes))
        DT::datatable(
          make_top_table(input = data()#, 
                         # gls = input$gls_table_top
          ) %>% 
            censor(censor_genes, censor_status$choice, censor_status$num_sim_genes) %>% 
            dplyr::mutate(gene = map_chr(gene, internal_link, linkout_img = FALSE)) %>% #from fun_helper.R
            dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", 
                          # "GLS p-value" = "GLSpvalue", 
                          "R^2" = "r2", "Z-Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index") %>% 
            dplyr::select("Query", "Gene", "Name", input$vars_dep_top) %>%
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25)
        )
      })
    }
  )
}

similarPathwaysTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_pos_enrich")))),
    fluidRow(DT::dataTableOutput(outputId = ns("pos_enrich")))
  )
}

similarPathwaysTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_pos_enrich <- renderText({paste0("Pathways of genes with similar dependencies as ", str_c(data()$content, collapse = ", "))})
      output$pos_enrich <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% master_positive$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_top(master_positive, input = data()) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          options = list(pageLength = 25))
      })  
    }
  )
}

similarCellsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cells_dep_top")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_cells_dep_top"),
                                label = "Select columns:",
                                choices = c("Lineage", "Sub-lineage", "Estimate",
                                            "P.Value", "Bonferroni", "Sex", "Age", "Status"),
                                selected = c("Lineage", "Estimate", "Bonferroni"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cells_dep_top"))), 
    tags$br()
  )
}

similarCellsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_cells_dep_top <- renderText({paste0(nrow(make_cell_sim_table(input = data())$top_table),
                                                      " cells with similar dependency profiles as ", 
                                                      str_c(data()$content, collapse = ", "))})      
      output$cells_dep_top <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data())$top_table) > 0,
               "No data found for this cell line"))
        DT::datatable(
          make_cell_sim_table(input = data())$top_table %>%
            dplyr::mutate(cell2_name = map_chr(cell2_name, cell_linkr, type = "cell") #from fun_helper.R
            ) %>%
            dplyr::rename("Query" = "cell1_name", "Cell" = "cell2_name", "Lineage" = "lineage", 
                          "Sub-lineage" = "lineage_subtype", "Estimate" = "coef", 
                          "P.Value" = "pval", "Bonferroni" = "bonferroni",
                          "Sex" = "sex", "Age" = "age", "Status" = "status") %>% 
            dplyr::arrange(Bonferroni) %>% 
            dplyr::select(Query, Cell, input$vars_cells_dep_top) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

similarExpCellsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cells_exp_top")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_cells_exp_top"),
                                label = "Select columns:",
                                choices = c("Lineage", "Sub-lineage", "Estimate",
                                            "P.Value", "Bonferroni", "Sex", "Age", "Status"),
                                selected = c("Lineage", "Estimate", "Bonferroni"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cells_exp_top"))), 
    tags$br()
  )
}

similarExpCellsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_cells_exp_top <- renderText({paste0(nrow(make_cell_sim_table(input = data(),
                                                                               similarity = "expression")$top_table),
                                                      " cells with similar gene expression profiles as ", 
                                                      str_c(data()$content, collapse = ", "))})      
      output$cells_exp_top <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data(),
                                        similarity = "expression")$top_table) > 0,
               "No data found for this cell line"))
        DT::datatable(
          make_cell_sim_table(input = data(),
                              similarity = "expression")$top_table %>%
            dplyr::mutate(cell2_name = map_chr(cell2_name, cell_linkr, type = "cell") #from fun_helper.R
            ) %>%
            dplyr::rename("Query" = "cell1_name", "Cell" = "cell2_name", "Lineage" = "lineage", 
                          "Sub-lineage" = "lineage_subtype", "Estimate" = "coef", 
                          "P.Value" = "pval", "Bonferroni" = "bonferroni",
                          "Sex" = "sex", "Age" = "age", "Status" = "status") %>% 
            dplyr::arrange(Bonferroni) %>% 
            dplyr::select(Query, Cell, input$vars_cells_exp_top) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

similarCompoundsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compound_dep_top")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_compound_dep_top"), 
                                "Select columns:",
                                c("Mechanism", "R^2", "Z-Score"), #"Co-publication Count", "Co-publication Index"
                                selected = c("Z-Score"), #, "Co-publication Count"
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("compound_dep_top"))), 
    tags$br()
  )
}

similarCompoundsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_compound_dep_top <- renderText({paste0(nrow(make_compound_table(input = data(), top = TRUE)), " compounds with similar dependencies as ", str_c(data()$content, collapse = ", "))})      
      output$compound_dep_top <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_compound_table(input = data(), top = TRUE)) != 0, "No data found for this compound"))
        DT::datatable(
          make_compound_table(input = data(), top = TRUE) %>%
            dplyr::mutate(name = map_chr(name, drug_linkr), #from fun_helper.R
                          r2 = round(r2, 2), 
                          z_score = round((r2 - mean_virtual_prism_cor)/sd_virtual_prism_cor, 1) 
            ) %>% 
            dplyr::rename("R^2" = "r2", "Z-Score" = "z_score", "Mechanism" = "moa") %>% 
            dplyr::select("Query" = "fav_drug", "Drug" = "name", input$vars_compound_dep_top),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

##Dissimilar----
dissimilarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_bottom")))),
    # shinyWidgets::prettySwitch(inputId = ns("gls_table_bottom"), 
    #                            "Show GLS table", 
    #                            value = FALSE),
    fluidRow(checkboxGroupInput(inputId = ns("vars_dep_bottom"), 
                                "Select columns:",
                                c("R^2", "Z-Score", "Co-publication Count", "Co-publication Index"), # "GLS p-value"
                                selected = c("R^2", "Co-publication Count"), # "GLS p-value"
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("dep_bottom")))
  )
}

dissimilarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      observeEvent(input$dsim_pathway_click, { #event to store the 'click'
      })
      output$text_dep_bottom <- renderText({paste0(nrow(make_bottom_table(input = data()#, 
                                                                          # gls = input$gls_table_bottom
      )
      ), " genes with dissimilar dependencies as ", str_c(data()$content, collapse = ", "))})      
      output$dep_bottom <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% master_bottom_table$fav_gene, 
               "No data found for this gene."))
        DT::datatable(
          make_bottom_table(input = data()#, 
                            # gls = input$gls_table_bottom
          ) %>%
            dplyr::mutate(gene = map_chr(gene, internal_link, linkout_img = FALSE)) %>%  #from fun_helper.R 
            #dplyr::mutate(link = map_chr(gene, internal_link, linkout_img = TRUE)) %>% 
            dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", 
                          # "GLS p-value" = "GLSpvalue", 
                          "R^2" = "r2", "Z-Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index") %>%
            dplyr::select("Query", "Gene", "Name", input$vars_dep_bottom) %>%
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

dissimilarCellsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cells_dep_bottom")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_cells_dep_bottom"),
                                label = "Select columns:",
                                choices = c("Lineage", "Sub-lineage", "Estimate",
                                            "P.Value", "Bonferroni", "Sex", "Age", "Status"),
                                selected = c("Lineage", "Estimate", "Bonferroni"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cells_dep_bottom"))), 
    tags$br()
  )
}

dissimilarCellsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_cells_dep_bottom <- renderText({paste0(nrow(make_cell_sim_table(input = data())$bottom_table), 
                                                         " cells with dissimilar dependency profiles as ", str_c(data()$content, collapse = ", "))})      
      output$cells_dep_bottom <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data())$bottom_table) > 0,
               "No data found for this cell line"))
        DT::datatable(
          make_cell_sim_table(input = data())$bottom_table %>%
            dplyr::mutate(cell2_name = map_chr(cell2_name, cell_linkr, type = "cell") #from fun_helper.R
            ) %>%
            dplyr::rename("Query" = "cell1_name", "Cell" = "cell2_name", "Lineage" = "lineage", 
                          "Sub-lineage" = "lineage_subtype", "Estimate" = "coef", 
                          "P.Value" = "pval", "Bonferroni" = "bonferroni",
                          "Sex" = "sex", "Age" = "age", "Status" = "status") %>% 
            dplyr::arrange(Bonferroni) %>% 
            dplyr::select(Query, Cell, input$vars_cells_dep_bottom) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

dissimilarExpCellsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cells_exp_bottom")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_cells_exp_bottom"),
                                label = "Select columns:",
                                choices = c("Lineage", "Sub-lineage", "Estimate",
                                            "P.Value", "Bonferroni", "Sex", "Age", "Status"),
                                selected = c("Lineage", "Estimate", "Bonferroni"),
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("cells_exp_bottom"))), 
    tags$br()
  )
}

dissimilarExpCellsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_cells_exp_bottom <- renderText({paste0(nrow(make_cell_sim_table(input = data(),
                                                                                  similarity = "expression")$bottom_table),
                                                         " cells with dissimilar expression profiles as ", 
                                                         str_c(data()$content, collapse = ", "))})      
      output$cells_exp_bottom <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data(),
                                        similarity = "expression")$bottom_table) > 0,
               "No data found for this cell line"))
        DT::datatable(
          make_cell_sim_table(input = data(),
                              similarity = "expression")$bottom_table %>%
            dplyr::mutate(cell2_name = map_chr(cell2_name, cell_linkr, type = "cell") #from fun_helper.R
            ) %>%
            dplyr::rename("Query" = "cell1_name", "Cell" = "cell2_name", "Lineage" = "lineage", 
                          "Sub-lineage" = "lineage_subtype", "Estimate" = "coef", 
                          "P.Value" = "pval", "Bonferroni" = "bonferroni",
                          "Sex" = "sex", "Age" = "age", "Status" = "status") %>% 
            dplyr::arrange(Bonferroni) %>% 
            dplyr::select(Query, Cell, input$vars_cells_exp_bottom) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}

dissimilarPathwaysTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_neg_enrich")))),
    fluidRow(DT::dataTableOutput(outputId = ns("neg_enrich")))
  )
}

dissimilarPathwaysTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_neg_enrich <- renderText({paste0("Pathways of genes with inverse dependencies as ", str_c(data()$content, collapse = ", "))})
      output$neg_enrich <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% master_negative$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_bottom(master_negative, input = data()) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)),
          options = list(pageLength = 25))
      })      
    }
  )
}

dissimilarCompoundsTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compound_dep_bottom")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_compound_dep_bottom"), 
                                "Select columns:",
                                c("Mechanism", "R^2", "Z-Score"), #"Co-publication Count", "Co-publication Index"
                                selected = c("Z-Score"), #, "Co-publication Count"
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("compound_dep_bottom"))), 
    tags$br()
  )
}

dissimilarCompoundsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_compound_dep_bottom <- renderText({paste0(nrow(make_compound_table(input = data()$content, top = FALSE)), " compounds with inverse dependencies as ", str_c(data()$content, collapse = ", "))})      
      output$compound_dep_bottom <- DT::renderDataTable({
        shiny::validate(
          need(nrow(make_compound_table(input = data()$content, top = FALSE)) != 0, "No data found for this compound"))
        DT::datatable(
          make_compound_table(input = data()$content, top = FALSE) %>%
            dplyr::mutate(name = map_chr(name, drug_linkr), #from fun_helper.R
                          r2 = round(r2, 2), 
                          z_score = round((r2 - mean_virtual_prism_cor)/sd_virtual_prism_cor, 1) 
            ) %>% 
            dplyr::rename("R^2" = "r2", "Z-Score" = "z_score", "Mechanism" = "moa") %>% 
            dplyr::select("Query" = "fav_drug", "Drug" = "name", input$vars_compound_dep_bottom),
          escape = FALSE,
          options = list(pageLength = 25))
      })
    }
  )
}
##Metabolites-----
#module that displays a table for metabolites

metabolitesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_metabolites_table")))),
    DT::dataTableOutput(outputId = ns("metabolites_table"))
  )
}

metabolitesTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_metabolites_table <- renderText({paste0("Metabolite table for ", str_c(data()$content, collapse = ", "))})
      output$metabolites_table <- DT::renderDataTable({
        if(data()$type == "gene") {
          shiny::validate(
            need(hmdb_proteins %>% 
                   filter(fav_gene %in% data()$content) %>% 
                   nrow() > 0, 
                 "No metabolites found that associate with this gene.")
          )
          DT::datatable(make_metabolite_table(input = data()) %>% 
                          dplyr::mutate(metabolite_name = map_chr(metabolite_name, metabolite_linkr)) %>% 
                          dplyr::select('Gene Name' = gene_name, 
                                        'Metabolite' = metabolite_name,
                                        'HMDB ID' = metabolite_accession),
                        escape = FALSE,
                        options = list(pageLength = 10))
        } else if(data()$type == "cell") {
          shiny::validate(
            need(cell_metabolites %>% 
                   left_join(expression_names, by = c("DepMap_ID" = "X1")) %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No metabolites found that associate with this cell line.")
          )
          DT::datatable(make_metabolite_table(input = data()) %>%
                          dplyr::mutate(metabolite = map_chr(metabolite, metabolite_linkr)) %>% 
                          dplyr::select('Cell Line' = cell_line, 
                                        'Metabolite' = metabolite, 
                                        'Value' = value) %>% 
                          mutate(Value = round(Value, 3)),
                        escape = FALSE,
                        options = list(pageLength = 10))
        }
      })
    }
  )
}

##Drug Tables -----
#search gene, find drugs
geneDrugsTable <- function(id) { #GENE QUERY
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_gene_drugs_table")))),
    fluidRow(DT::dataTableOutput(outputId = ns("gene_drugs_table"))), 
    br(),
    fluidRow(textOutput(ns("text_gene_drugs_table")))
  )
}

geneDrugsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_gene_drugs_table <- renderText({paste0("Drug sensitivity table for ", str_c(data()$content, collapse = ", "))})
      output$text_gene_drugs_table <- renderText({paste0("Drugs annotated for ", str_c(data()$content, collapse = ", "))})
      output$gene_drugs_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% gene_drugs_table$fav_gene, "No drug data found for this gene."))
        DT::datatable(make_gene_drugs_table(input = data()) %>% 
                        dplyr::mutate(fav_drug = map_chr(fav_drug, drug_linkr), 
                                      moa = map_chr(moa, moa_linkr)) %>% #from fun_helper.R
                        dplyr::rename(Gene = fav_gene, Drug = fav_drug, Mechanism = moa), 
                      escape = FALSE)
      })
    }
  )
}

cellDrugsTable <- function(id) { 
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_cell_drugs_table")))),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_drugs_table"))), 
    br(),
    fluidRow(textOutput(ns("text_cell_drugs_table")))
  )
}

cellDrugsTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_cell_drugs_table <- renderText({paste0("Drug sensitivity table for ", str_c(data()$content, collapse = ", "))})
      output$text_cell_drugs_table <- renderText({paste0("Drugs annotated for ", str_c(data()$content, collapse = ", "))})
      output$cell_drugs_table <- DT::renderDataTable({
        shiny::validate(
          need(prism_long %>% 
                 dplyr::rename(X1 = 1) %>% 
                 left_join(expression_meta, by = "X1") %>% 
                 pull(cell_line) %in% data()$content,
               "No drug data found for this cell line."))
        DT::datatable(make_cell_drugs_table(input = data()) %>% 
                        dplyr::mutate(name = map_chr(name, drug_linkr),
                                      log2fc = round(log2fc, 3)) %>% #from fun_helper.R
                        dplyr::rename(`Cell Line` = cell_line, Drug = name, LogFC = log2fc), 
                      escape = FALSE)
      })
    }
  )
}

#cor tables
geneDrugsCorTable <- function(id) { #GENE QUERY
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_gene_drugs_cor_table")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_gene_drugs_cor_table"), 
                                "Select columns:",
                                c("Known", "R^2", "Z-Score"), 
                                selected = c("Z-Score"), 
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("gene_drugs_cor_table"))), 
    br(),
    fluidRow(textOutput(ns("text_gene_drugs_cor_table")))
  )
}

geneDrugsCorTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_gene_drugs_cor_table <- renderText({glue::glue('Drug correlation table for {str_c(data()$content, collapse = ", ")}')})
      output$text_gene_drugs_cor_table <- renderText({glue::glue('When {str_c(data()$content, collapse = ", ")} is knocked out, a subset of cells die. These are the drugs that show the same cell killing profile.')})
      output$gene_drugs_cor_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% gene_drugs_cor_table$fav_gene, "No drug targets found for this gene in this dataset."))
        #gene_drug_list will grab a char vec of drugs known to target genes, so we can check if corr's are known; added ifelse logic so index grab doesn't break
        gene_drug_tibble <- gene_drugs_table %>% filter(fav_gene %in% data()$content)
        gene_drug_list <- ifelse(nrow(gene_drug_tibble) != 0, gene_drug_tibble %>% unnest(data) %>% pull(fav_drug), "")
        DT::datatable(make_gene_drugs_cor_table(input = data()) %>% 
                        dplyr::mutate(
                          known = case_when(
                            drug %in% gene_drug_list ~ "TRUE", 
                            TRUE ~ "FALSE"),
                          drug = map_chr(drug, drug_linkr) #from fun_helper.R
                        ) %>% 
                        dplyr::rename(Gene = fav_gene, Drug = drug, Mechanism = moa, 'Known' = known, 'Z-Score' = z_score, 'R^2'=r2) %>% 
                        dplyr::select("Gene", "Drug", "Mechanism", input$vars_gene_drugs_cor_table), 
                      escape = FALSE)
      })
    }
  )
}

#CELL-----
cellSummaryTable <- function (id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("cell_sumary_title")))),
    fluidRow(DT::dataTableOutput(outputId = ns("cell_table")))
  )
}

cellSummaryTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_sumary_title <- renderText({paste0("Lineage table for ", str_c(data()$content, collapse = ", "))})
      output$cell_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% expression_meta$cell_line, 
          "No data found for this cell line.")
        )
        DT::datatable(
          make_cell_line_table(input = data()) %>% 
            dplyr::mutate(`Cell Line` = map_chr(`Cell Line`, cell_linkr, type = "cell")),
                      escape = FALSE,
                      options = list(paging = FALSE, 
                                     searching = FALSE,
                                     pageLength = 10))
      })      
    }
  )
}

#COMPOUND-----
##Pubmed-----
pubmedCompoundTable <- function(id) {
  ns <- NS(id)
  DT::dataTableOutput(outputId = ns("pubmed_compound_table"))
}

pubmedCompoundTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pubmed_compound_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% pubmed$name, ""))
        withProgress(message = 'Building a smart table...', {
          DT::datatable(make_pubmed_table(pubmed, input = data()) %>% 
                          dplyr::mutate(pmid = map_chr(pmid, pubmed_linkr, number_only = TRUE) #from fun_helper.R
                          ) %>% 
                          dplyr::mutate(pmcid = map_chr(pmcid, pmc_linkr) #from fun_helper.R
                          ) %>% 
                          dplyr::rename(Name = name, 'Pubmed ID' = pmid, Year = year, PMCID = pmcid),
                        escape = FALSE,
                        options = list(pageLength = 10))
        })
      })
    }
  )
}

pubmedCellLineTable <- function(id) {
  ns <- NS(id)
  DT::dataTableOutput(outputId = ns("pubmed_cell_line_table"))
}

pubmedCellLineTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pubmed_cell_line_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% pubmed$name, ""))
        withProgress(message = 'Building a smart table...', {
          DT::datatable(make_pubmed_table(pubmed, input = data()) %>% 
                          dplyr::mutate(pmid = map_chr(pmid, pubmed_linkr, number_only = TRUE) #from fun_helper.R
                          ) %>% 
                          dplyr::mutate(pmcid = map_chr(pmcid, pmc_linkr) #from fun_helper.R
                          ) %>% 
                          dplyr::rename(Name = name, 'Pubmed ID' = pmid, Year = year, PMCID = pmcid),
                        escape = FALSE,
                        options = list(pageLength = 10))
        })
      })
    }
  )
}

##Drugs-----

drugGenesTable <- function(id) { #DRUG QUERY
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_drug_genes_table")))),
    fluidRow(DT::dataTableOutput(outputId = ns("drug_genes_table"))), 
    br(),
    fluidRow(textOutput(ns("text_drug_genes_table")))
  )
}

drugGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_drug_genes_table <- renderText({glue::glue('Gene target table for {str_c(data()$content, collapse = ", ")}')})
      output$text_drug_genes_table <- renderText({glue::glue('Genes annotated to be targeted by {str_c(data()$content, collapse = ", ")}')})
      output$drug_genes_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% drug_genes_table$fav_drug, "No gene data found for this compound."))
        DT::datatable(make_drug_genes_table(drug = data()$content) %>% 
                        dplyr::mutate(fav_gene = map_chr(fav_gene, internal_link)) %>% #from fun_helper.R
                        dplyr::rename(Drug = fav_drug, Gene = fav_gene, Name = approved_name), 
                      escape = FALSE)
      })
    }
  )
}

drugGenesCorTable <- function(id) { #DRUG QUERY
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_drug_genes_cor_table")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_drug_genes_cor_table"), 
                                "Select columns:",
                                c("Known", "R^2", "Z-Score"), 
                                selected = c("Z-Score"), 
                                inline = TRUE)),
    fluidRow(DT::dataTableOutput(outputId = ns("drug_genes_cor_table"))), 
    br(),
    fluidRow(textOutput(ns("text_drug_genes_cor_table")))
  )
}

drugGenesCorTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_drug_genes_cor_table <- renderText({glue::glue('Gene correlation table for {str_c(data()$content, collapse = ", ")}')})
      output$text_drug_genes_cor_table <- renderText({glue::glue('When {str_c(data()$content, collapse = ", ")} is placed on cells, a subset dies. These are the genes that show the same cell dependency profile.')})
      output$drug_genes_cor_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% drug_genes_cor_table$fav_drug, "No gene data found for this compound in this dataset."))
        #drug_gene_list will grab a char vec of genes known to be targeted compound, so we can check if corr's are known; added ifelse logic so index grab doesn't break
        drug_gene_tibble <- drug_genes_table %>% filter(fav_drug %in% data()$content)
        drug_gene_list <- ifelse(nrow(drug_gene_tibble) != 0, drug_gene_tibble %>% unnest(data) %>% pull(fav_gene), "")
        DT::datatable(make_drug_genes_cor_table(drug = data()$content) %>% 
                        dplyr::mutate(
                          known = case_when(
                            gene %in% drug_gene_list ~ "TRUE", 
                            TRUE ~ "FALSE"),
                          gene = map_chr(gene, internal_link) #from fun_helper.R
                        ) %>% 
                        dplyr::rename(Drug = fav_drug, Gene = gene, Name = approved_name, 'Known' = known, 'Z-Score' = z_score, 'R^2'=r2) %>% 
                        dplyr::select("Drug", "Gene", "Name", input$vars_drug_genes_cor_table), 
                      escape = FALSE)
      })
    }
  )
}


##Metabolites----

metaboliteGenesTable <- function(id) {
  ns <- NS(id)
  DT::dataTableOutput(outputId = ns("metabolite_genes_table"))
}

metaboliteGenesTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$metabolite_genes_table <- DT::renderDataTable({
        shiny::validate(
          need(data()$content %in% hmdb_metabolites$fav_metabolite, 
               "No genes found that associate with this metabolite"))
        DT::datatable(make_metabolite_table(input = data()) %>% 
                        dplyr::mutate(gene_name = map_chr(gene_name, internal_link)) %>% 
                        dplyr::select('Metabolite' = metabolite_name, 'Gene Name' = gene_name),
                      escape = FALSE,
                      options = list(pageLength = 10))
      })
    }
  )
}

