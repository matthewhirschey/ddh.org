
# conditionalPanel that only shows field is not 0
# using array bracket notation to be compatible with namespaces
notZeroConditionPanel <- function(fieldname, ...) {
  condition_str <- paste0("input['", fieldname, "'] != 0")
  conditionalPanel(condition = condition_str, ...)  
}

#Cell Dep Table -----
cellDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_dep_table")))),
    fluidRow(checkboxInput(inputId = ns("dep_filter_click"), label = "Filter dependency table", value = FALSE)),
    fluidRow(dataTableOutput(outputId = ns("target_achilles"))))
}

cellDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$target_achilles <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        DT::datatable(make_achilles_table(achilles, expression_names, data()$gene_symbols), 
                      filter = if(input$dep_filter_click == FALSE) {'none'} else {'top'}, 
                      options = list(pageLength = 10))
      })
    }
  )
}

#Similar----
similarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_top")))),
    fluidRow(
      column(6, checkboxGroupInput(inputId = ns("vars_dep_top"), 
                                   "Select columns:",
                                   c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                   selected = c("Z Score", "Co-publication Count"), 
                                   inline = TRUE)),
      column(6, fluidRow(sliderInput(inputId = ns("num_sim_genes"),
                                     "Remove genes with <n associations:",
                                     min = 100,
                                     max = 1000,
                                     value = 1000,
                                     step = 100)), 
             fluidRow(column(3, actionButton(inputId = ns("censor"), "Submit")),
                      column(3, actionButton(inputId = ns("reset"), "Reset"))))),
    hr(),
    fluidRow(dataTableOutput(outputId = ns("dep_top"))), 
    tags$br(),
    fluidRow(actionLink(inputId = ns("sim_pathway_click"), "View enriched pathways for positive correlations")), #add conditional panel for sim_pathways
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("sim_pathway_click"), "'] != 0"), 
                     tags$br(),
                     fluidRow(h4(textOutput(ns("text_pos_enrich")))),
                     fluidRow(dataTableOutput(outputId = ns("pos_enrich"))))
  )
}

similarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      #censor reactive values
      censor_status <- reactiveValues(choice = FALSE, 
                                      num_sim_genes = 1000, 
                                      num = "All" #nrow(make_top_table(gene_symbol = data()))
      )
      observeEvent(input$censor, {
        censor_status$choice <- TRUE
        censor_status$num_sim_genes <- input$num_sim_genes
        censor_status$num <- nrow(make_top_table(gene_symbol = data()$gene_symbols) %>% censor(censor_genes, censor_status$choice, censor_status$num_sim_genes))
      })
      
      observeEvent(input$reset, {
        censor_status$choice <- FALSE
        censor_status$num_sim_genes <- 1000
        updateSliderInput(session, inputId = "num_sim_genes", value = 1000)
        censor_status$num <- nrow(make_top_table(gene_symbol = data()$gene_symbols))
      })
      observeEvent(input$sim_pathway_click, { #event to store the 'click'
      })
      output$text_dep_top <- renderText({paste0(censor_status$num, " genes with similar dependencies as ", str_c(data()$gene_symbols, collapse = ", "))})      
      output$dep_top <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% master_top_table$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_top_table(master_top_table, data()$gene_symbols) %>% 
            dplyr::mutate(link = paste0("<center><a href='?show=gene&query_type=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
            dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_top) %>%
            censor(censor_genes, censor_status$choice, censor_status$num_sim_genes),
          escape = FALSE,
          options = list(pageLength = 25))})
      output$text_pos_enrich <- renderText({paste0("Pathways of genes with similar dependencies as ", str_c(data()$gene_symbols, collapse = ", "))})
      output$pos_enrich <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% master_positive$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_top(master_positive, data()$gene_symbols),
          options = list(pageLength = 25))
      })  
    }
  )
}

#Dissimilar----
dissimilarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_bottom")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_dep_bottom"), 
                                "Select columns:",
                                c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                selected = c("Z Score", "Co-publication Count"), 
                                inline = TRUE)),
    fluidRow(dataTableOutput(outputId = ns("dep_bottom"))), 
    tags$br(),
    fluidRow(actionLink(inputId = ns("dsim_pathway_click"), "View enriched pathways below")), #add conditional panel for sim_pathways
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("dsim_pathway_click"), "'] != 0"), 
                     tags$br(),
                     fluidRow(h4(textOutput(ns("text_neg_enrich")))),
                     fluidRow(dataTableOutput(outputId = ns("neg_enrich")))
                     )
  )
}

dissimilarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      observeEvent(input$dsim_pathway_click, { #event to store the 'click'
      })
      output$text_dep_bottom <- renderText({paste0(nrow(make_bottom_table(gene_symbol = data()$gene_symbols)), " genes with inverse dependencies as ", str_c(data()$gene_symbols, collapse = ", "))})      
      output$dep_bottom <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% master_bottom_table$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_bottom_table(master_bottom_table, data()$gene_symbols) %>%
            dplyr::mutate(link = paste0("<center><a href='?show=gene&query_type=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
            dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_bottom),
          escape = FALSE,
          options = list(pageLength = 25))
      })
      output$text_neg_enrich <- renderText({paste0("Pathways of genes with inverse dependencies as ", str_c(data()$gene_symbols, collapse = ", "))})
      output$neg_enrich <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% master_negative$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_bottom(master_negative, data()$gene_symbols),
          options = list(pageLength = 25))
      })      
    }
  )
}

#Browse Pathways----
# module that displays a table of pathways when an link is clicked

browsePathwaysLink <- function (id) {
  ns <- NS(id)
  actionLink(inputId = ns("pathway_click"), "browse the pathways")
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
  notZeroConditionPanel(ns("pathway_click"),
                   tags$br(),
                   h4("GO Biological Processes"),
                   dataTableOutput(outputId = ns("pathway_table")))
}

browsePathwaysPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_table <- DT::renderDataTable({
        DT::datatable(make_pathway_table(pathways) %>% dplyr::rename(Pathway = pathway, GO = go, Genes = genes), 
                      options = list(pageLength = 10))
      })      
    }
  )
}

#Cell anatogram----
# module that displays a table for cell anatogram

cellAnatogramTable <- function(id) {
  ns <- NS(id)
    dataTableOutput(outputId = ns("cellanatogram_table"))
}

cellAnatogramTableServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogram_table <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% subcell$gene_name, ""))
        DT::datatable(make_cellanatogram_table(subcell, data()$gene_symbols),
                      options = list(pageLength = 10))
      })
    }
  )
}

#Cell Expression Table -----
cellExpressionTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_exp_table")))),
    fluidRow(dataTableOutput(outputId = ns("cell_expression"))))
}

cellExpressionTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_exp_table <- renderText({paste0("Cell expression table for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cell_expression <- DT::renderDataTable({
        validate(
          need(data()$gene_symbols %in% colnames(expression), "No expression data found for this gene."))
        make_expression_table(gene_symbol = data()$gene_symbols)
      })
    }
  )
}