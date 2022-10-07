compoundPage <- function (id, subtype) {
  ns <- NS(id)
  
  summary_var <- ""
  
  #this block is the logic to define the summary_var variable, to display the proper summary module
  if(subtype == "compound"){
    title_var <- compoundTitle(ns("title_var"))
    summary_var <- compoundSummaryText(ns("summary"))
    summary_table <- "hi"#pathwayList(ns("gene_pathways"))
  } else if (subtype == "moa"){
    title_var <- "hi"#pathwayTitle(ns("title_var"))
    summary_var <- moaSummaryText(ns("summary"))
    summary_table <- "hi"#pathwayGeneList(ns("gene_pathways"))
  } else if (subtype == "compound_list") {
    title_var <- "hi"#geneListTitle(ns("title_var"))
    summary_var <- compoundListSummaryText(ns("summary"))
    summary_table <- "hi"#pathwayList(ns("gene_pathways"))
  } else if (subtype == "metabolite") {
    title_var <- metaboliteTitle(ns("title_var"))
    summary_var <- "TODO" # metaboliteSummaryText(ns("summary"))
    summary_table <- "TODO"
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      id=ns("compoundNavBar"),
      ## DASHBOARD-----
      tabPanel("SUMMARY",
               title_var,
               flowLayout(
                 actionLink(inputId = ns("link_to_compoundStructureDash"), compoundStructureDash(ns("structuredash"))), 
                 actionLink(inputId = ns("link_to_pubmedCompoundPlotDash"), pubmedCompoundPlotDash(ns("pubmedcompounddash"))),
                 private(actionLink(inputId = ns("link_to_compoundDependenciesPlotDash"), compoundDependenciesPlotDash(ns("compounddepdash")))), 
                 private(actionLink(inputId = ns("link_to_compoundDependenciesTableDash"), compoundDependenciesTableDash(ns("compoundtabledash")))), 
                 private(actionLink(inputId = ns("link_to_compoundDependenciesGraphDash"), compoundDependenciesGraphDash(ns("compoundgraphdash")))), 
                 private(actionLink(inputId = ns("link_to_drugGenesCorTableDash"), drugGenesCorTableDash(ns("druggenescortabledash"))))
               )
      ),
      ## INFO-----
      navbarMenu("INFO",
                 tabPanel(title = "Compound", value = "about_compound",
                          fluidRow(
                            column(8, summary_var), #summary variable for alt descriptions
                            column(4, compoundStructure(ns("structure")))
                          ),
                          tags$br(),
                          tags$hr(),
                          fluidRow(
                            #summary_table
                          )
                 ),
                 tabPanel(title = "Literature", value = "about_literature",
                          #co-citation graph? 
                          pubmedCompoundPlot(ns("pubmed_compound")),  #Summary stats
                          pubmedCompoundTable(ns("pubmed_compound"))
                 )
      ),
      ## EXPRESSION-----
      navbarMenu(title = "EXPRESSION",
                 tabPanel("Sub-cellular", value = "expression_sub",
                          "still thinking", 
                          private_msg()
                 ),
                 tabPanel("Cell Line", value = "expression_cell",
                          "still thinking", 
                          private_msg()
                 ),
                 tabPanel("Tissue", value = "expression_tissue",
                          "still thinking", 
                          private_msg()
                 ),
                 tabPanel("Coexpression", value = "expression_co",
                          "still thinking", 
                          private_msg()
                 )
      ),
      # COMPOUNDS-----
      navbarMenu(title = "COMPOUNDS",
                 tabPanel("Genes", value = "genes",
                          private(drugGenesTable(ns("drug_genes"))), 
                          private_msg()
                 ),
                 tabPanel("Metabolites", value = "metabolites",
                          private(metaboliteGenesTable(ns("metabolite_genes"))), 
                          private_msg()),
                 tabPanel("Graph", value = "compound_bipartite_graph",
                          private(compoundBipartiteGraph(ns("compound_bipartite_graph"))), 
                          private_msg())
      ),
      ## DEPENDENCIES (action)-----
      navbarMenu(title = "DEPENDENCIES",
                 tabPanel("Plots", value = "compound_dependencies_plots", 
                          private(compoundDependenciesPlot(ns("dep"))),
                          private_msg(),
                          tags$hr(),
                          private(compoundBinsPlot(ns("dep"))),
                          private_msg(),
                          tags$hr(),
                          private(compoundLinPlot(ns("dep"))), 
                          private(compoundSubLinPlot(ns("dep"))), 
                          private_msg(),
                          tags$hr(),
                          private(fluidRow(actionLink(inputId = ns("compound_dep_click"), "View raw dependency data table below"))), # conditional panel for raw data
                          tags$br(),
                          private(conditionalPanel(condition = paste0("input['", ns("compound_dep_click"), "'] != 0"),
                                                   compoundDependenciesTable(ns("dep"))))
                 ),
                 tabPanel("Co-essentiality", value = "compound_dependencies_correlation",
                          private(fluidRow(actionLink(inputId = ns("compound_corrplot_click"), "Build compound correlation plot"))), # conditional panel for correlation plot
                          private_msg(),
                          tags$br(), 
                          private(conditionalPanel(condition = paste0("input['", ns("compound_corrplot_click"), "'] != 0"), 
                                                   compoundCorrelationPlot(ns("compound_corrplot")), 
                                                   tags$hr())),
                          private(similarCompoundsTable(ns("compound_sim"))),
                          tags$br(),
                          private(fluidRow(actionLink(inputId = ns("cell_ess_compound_click"), "View inverse dependency data table below"))), # conditional panel for raw data
                          tags$br(),
                          private(conditionalPanel(condition = paste0("input['", ns("cell_ess_compound_click"), "'] != 0"),
                                                   dissimilarCompoundsTable(ns("compound_dsim")))),
                          tags$br()
                 ),
                 tabPanel("Graph", value = "compound_dependencies_graph",
                          private(compoundNetworkGraph(ns("compound_graph"))),
                          private_msg()
                 ),
                 tabPanel("Genes", value = "dependencies_genes", 
                          private(drugGenesCorTable(ns("drug_genes_cor"))),
                          private_msg()
                 )
      ),
      ## DOWNLOADS-----
      downloadTab(ns("download")), #pull download tab from module
      
      ## SEARCH-----
      formContent=querySearchInput(ns("search"))
    )
  )
}

compoundPageServer <- function(id, subtype) {
  moduleServer(
    id,
    function(input, output, session) {
      type <- "compound"
      if(subtype == "compound") {
        data <- reactive({
          compound_name <- getQueryString()$query
          list(
            type=type,
            subtype=subtype,
            query=compound_name,
            content=compound_name
          )
        })
        title_var <- compoundTitleServer("title_var", data)
        summary_var <- compoundSummaryTextServer("summary", data)
      }
      if(subtype == "metabolite") {
        data <- reactive({
          metabolite_name <- getQueryString()$query
          list(
            type=type,
            subtype=subtype,
            query=metabolite_name,
            id=metabolite_name
          )
        })
        title_var <- metaboliteTitleServer("title_var", data)
        # TODO summary_var <- metaboliteSummaryTextServer("summary", data)
      }
      if(subtype == "moa") {
        data <- reactive({
          moa_str <- getQueryString()$query
          if (!is.null(moa_str)) {
            prism_rows <- 
              prism_names %>%
              filter(moa == moa_str)
            list(
              type=type,
              subtype=subtype,
              query=moa_str,
              content=prism_rows$name
            )
          }
        })
        #title_var <- moaTitleServer("title_var", data)
        #summary_var <- moaSummaryTextServer("summary", data)
        
      }
      if(subtype == "compound_list") {
        data <- reactive({
          custom_compound_list <- getQueryString()$query
          compounds_list <- c(str_split(custom_compound_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            subtype=subtype,
            query=custom_compound_list,
            content=compounds_list
          )
        })
        #title_var <- compoundListTitleServer("title_var", data)
        #summary_var <- compoundListSummaryTextServer("summary", data)
        
      }
      ## SEARCH SERVER-----
      querySearchServer("search")
      ## DASHBOARD SERVER-----
      #name pulls from if/else variables above
      #compound plot
      observeEvent(input$link_to_compoundStructureDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "about_compound")
      })
      compoundStructureDashServer("structuredash", data)
      
      #pubmed compound plot
      observeEvent(input$link_to_pubmedCompoundPlotDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "about_literature")
      })
      pubmedCompoundPlotDashServer("pubmedcompounddash", data)
      
      #dep plot
      observeEvent(input$link_to_compoundDependenciesPlotDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "compound_dependencies_plots")
      })
      compoundDependenciesPlotDashServer("compounddepdash", data)
      
      #compound dep table
      observeEvent(input$link_to_compoundDependenciesTableDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "compound_dependencies_correlation")
      })
      compoundDependenciesTableDashServer("compoundtabledash", data)
      
      #compound graph
      observeEvent(input$link_to_compoundDependenciesGraphDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "compound_dependencies_graph")
      })
      compoundDependenciesGraphDashServer("compoundgraphdash", data)
      
      #dep drug genes graph
      observeEvent(input$link_to_drugGenesCorTableDash, {
        updateNavbarPage(session, inputId = "compoundNavBar", selected = "dependencies_genes")
      })
      drugGenesCorTableDashServer("druggenescortabledash", data)
      
      ## INFO SERVER-----
      # Compound Structure
      compoundStructureServer("structure", data)
      # Literature
      pubmedCompoundPlotServer("pubmed_compound", data)
      pubmedCompoundTableServer("pubmed_compound", data)
      
      ## EXPRESSION SERVER-----
      
      # MOLECULES SERVER-----
      private({drugGenesTableServer("drug_genes", data)})
      private({metaboliteGenesTableServer("metabolite_genes", data)})
      private({compoundBipartiteGraphServer("compound_bipartite_graph", data)})
      
      # DEPENDENCIES SERVER----
      private({compoundDependenciesPlotServer("dep", data)})
      private({compoundBinsPlotServer("dep", data)})
      private({compoundLinPlotServer("dep", data)})
      private({compoundSubLinPlotServer("dep", data)})
      observeEvent(input$compound_dep_click, { #event to store the dependency table 'click'
      })
      private({compoundDependenciesTableServer("dep", data)})
      # Similar - Genes
      observeEvent(input$compound_corrplot_click, { #event to store the correlation plot 'click'
      })
      private({compoundCorrelationPlotServer("compound_corrplot", data)})
      private({similarCompoundsTableServer("compound_sim", data)})
      # Dissimilar - Genes
      observeEvent(input$cell_ess_compound_click, { #event to store the essentiality 'click'
      })
      private({dissimilarCompoundsTableServer("compound_dsim", data)})
      # Graph
      private({compoundNetworkGraphServer("compound_graph", data)})
      # Drug Cor
      private({drugGenesCorTableServer("drug_genes_cor", data)})
      
      # DOWNLOAD SERVER-----
      downloadTabServer("download", data, privateMode) #pull download tab from module
      
    }
  )
}
