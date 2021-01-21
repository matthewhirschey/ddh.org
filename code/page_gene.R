#this page is for the master gene query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page
genePage <- function (id, type) {
  ns <- NS(id)
  
  #this block is the logic to define the summary_var variable, to display the proper summary module
  if(type == "gene"){
    summary_var <- geneSummaryText(ns("summary"))
    gene_var <- geneText(ns("gene_var"))
    protein_summary <- proteinText(ns("protein_summary"))
  } else if (type == "pathway"){
    summary_var <- pathwaySummaryText(ns("summary"))
    gene_var <- nameText(ns("gene_var")) #nameText is a dummy placeholder
    protein_summary <- nameText(ns("protein_summary")) #nameText is a dummy placeholder
  } else if (type == "gene_list") {
    summary_var <- geneListSummaryText(ns("summary")) 
    gene_var <- nameText(ns("gene_var")) #nameText is a dummy placeholder
    protein_summary <- nameText(ns("protein_summary")) #nameText is a dummy placeholder
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      id=ns("geneNavBar"),
      tabPanel("DASHBOARD",
               fluidRow(
                 column(width = 4, 
                        actionLink(inputId = ns("link_to_cellDependenciesPlotDash"), cellDependenciesPlotDash(ns("depdash")))
                 )
               )
      ),
      navbarMenu("SUMMARY",
                 tabPanel(title = "Gene", value = "summary_gene",
                          gene_var #summary variable for alt descriptions
                          ), #change to navbarMenu when you have a submenu
                 tabPanel(title = "Protein", 
                          protein_summary, value = "summary_protein"
                          )
      ),
      navbarMenu(title = "EXPRESSION",
                 tabPanel("Gene", value = "expression_gene",
                          cellAnatogramPlot(ns("exp")),
                          cellAnatogramTable(ns("exp")), 
                          tags$hr(), 
                          cellExpressionPlot(ns("cell_exp")),
                          tags$br(),
                          fluidRow(tags$strong(plot_cellexp_title), plot_cellexp_legend,
                                   actionLink(inputId = ns("cell_exp_click"), "View raw expression data table")), #add conditional panel for raw data
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("cell_exp_click"), "'] != 0"), 
                                           cellExpressionTable(ns("cell_exp")))
      )),
      navbarMenu(title = "DEPENDENCIES",
                 tabPanel("Plots", value = "dependencies_plots", 
                          cellDependenciesPlot(ns("dep")),
                          tags$hr(),
                          cellBinsPlot(ns("dep")),
                          tags$hr(),
                          cellDepsLinPlot(ns("dep")), 
                          cellDepsSubLinPlot(ns("dep")), 
                          tags$hr(),
                          fluidRow(actionLink(inputId = ns("cell_dep_click"), "View raw dependency data table below")), # conditional panel for raw data
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("cell_dep_click"), "'] != 0"), 
                                           cellDependenciesTable(ns("dep")))
                 ),
                 tabPanel("Co-essentiality", value = "dependencies_co-essentiality",
                          fluidRow(actionLink(inputId = ns("corrplot_click"), "Build gene correlation plot")), # conditional panel for correlation plot
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("corrplot_click"), "'] != 0"), 
                                           geneCorrelationPlot(ns("corrplot")), 
                                           tags$hr()),
                          similarGenesTable(ns("sim")),
                          tags$br(),
                          fluidRow(actionLink(inputId = ns("cell_ess_click"), "View inverse dependency data table below")), # conditional panel for raw data
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("cell_ess_click"), "'] != 0"), 
                                           dissimilarGenesTable(ns("dsim"))), 
                          tags$br()),
                 tabPanel("Graph", value = "dependencies_graph",
                          geneNetworkGraph(ns("graph")))
      ),
      tabPanel("METHODS",
               includeHTML(here::here("code","methods.html"))),
      tabPanel("DOWNLOADS",
               downloadReportPanel(ns("download"))),
      formContent=querySearchInput(ns("search")))
  )
}

genePageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {
      
      #this block is the logic to define the summary_var variable, to display the proper summary module
      if(type == "gene"){
        data <- reactive({
          gene_symbol <- getQueryString()$symbol
          list(
            type=type,
            id=gene_symbol,
            gene_symbols=gene_symbol
          )
        })
        summary_var <- geneSummaryTextServer("summary", data)
        gene_var <- geneTextServer("gene_var", data)
        protein_summary <- proteinTextServer("protein_summary", data)
      } else if (type == "pathway"){
        data <- reactive({
          pathway_go <- getQueryString()$go
          pathway_row <- pathways %>%
            filter(go == pathway_go)
          list(
            type=type,
            id=pathway_go,
            gene_symbols=pathway_row$data[[1]]$gene)
        })
        summary_var <- pathwaySummaryTextServer("summary", data)
        gene_var <- nameTextServer("gene_var", data)
        protein_summary <- nameTextServer("protein_summary", data)
      } else if (type == "gene_list") {
        data <- reactive({
          custom_gene_list <- getQueryString()$custom_gene_list
          gene_symbols <- c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_gene_list,
            gene_symbols=gene_symbols
          )
        })
        summary_var <- geneListSummaryTextServer("summary", data)
        gene_var <- nameTextServer("gene_var", data)
        protein_summary <- nameTextServer("protein_summary", data)
      } else {
        stop("fix your summary argument")
      }
      
      querySearchServer("search")
      # DASHBOARD
      observeEvent(input$link_to_cellDependenciesPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "dependencies_plots")
      })
      cellDependenciesPlotDashServer("depdash", data)
      # SUMMARY
      summary_var
      # Gene
      gene_var
      # Protein
      protein_summary
      # EXPRESSION
      cellAnatogramPlotServer("exp", data)
      cellAnatogramTableServer("exp", data)
      cellExpressionPlotServer("cell_exp", data)
      observeEvent(input$cell_exp_click, { #event to store the expression table 'click'
      })
      cellExpressionTableServer("cell_exp", data)
      # DEPENDENCIES
      observeEvent(input$corrplot_click, { #event to store the correlation plot 'click'
      })
      geneCorrelationPlotServer("corrplot", data)
      cellDependenciesPlotServer("dep", data)
      cellBinsPlotServer("dep", data)
      cellDepsLinPlotServer("dep", data)
      observeEvent(input$cell_dep_click, { #event to store the dependency table 'click'
      })
      cellDepsSubLinPlotServer("dep", data)
      cellDependenciesTableServer("dep", data)
      
      # Similar - Genes
      similarGenesTableServer("sim", data)
      # Dissimilar - Genes
      observeEvent(input$cell_ess_click, { #event to store the essentiality 'click'
      })
      dissimilarGenesTableServer("dsim", data)
      # Graph
      geneNetworkGraphServer("graph", data)
      # DOWNLOAD
      downloadReportPanelServer("download", type, data)
    }
  )
}
