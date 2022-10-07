cellPage <- function (id, subtype) {
  ns <- NS(id)
  
  if (subtype == "cell") {
    title_var <- cellTitle(ns("title"))
    summary_var <- summaryListText(ns("summary"))
  } else if (subtype == "lineage") {
    title_var <- lineageTitle(ns("title"))
    summary_var <- lineageSummaryText(ns("summary"))
  } else if (subtype == "lineage_subtype") {
    title_var <- lineageSubtypeTitle(ns("title"))
    summary_var <- lineageSubtypeSummaryText(ns("summary"))
  } else if (subtype == "cell_list") {
    title_var <- summaryListTitle(ns("title"))
    summary_var <- summaryListText(ns("summary"))
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      id=ns("cellNavBar"),
      ## DASHBOARD-----
      tabPanel("SUMMARY",
               title_var,
               cardLayout(
                 actionLink(inputId = ns("link_tocellImagePlotDash"), cellImageDash(ns("cellimagedash"))),
                 actionLink(inputId = ns("link_to_pubmedCellPlotDash"), pubmedPlotDash(ns("pubmedcelldash"))),
                 private(actionLink(inputId = ns("link_to_cellExpressionPlotDash"), cellExpressionPlotDash(ns("cellexpressiondash")))), 
                 private(actionLink(inputId = ns("link_to_cellLineCoexpression"), cellLineExpressionPosTableDash(ns("cellexptabledash")))),
                 private(actionLink(inputId = ns("link_to_cellLineExpressionGraphDash"), cellExpressionGraphDash(ns("cellexpgraphdash")))),
                 private(actionLink(inputId = ns("link_to_cellLineFunctionalPlotDash"), cellFunctionalPlotDash(ns("cellfunplotdash")))),
                 private(actionLink(inputId = ns("link_to_cellLineMetaboliteTable"), cellMetabolitesTableDash(ns("cellmetabolitetabledash")))),
                 private(actionLink(inputId = ns("link_to_cellLineDrugTable"), cellDrugsTableDash(ns("celldrugtabledash")))),
                 private(actionLink(inputId = ns("link_to_cellLineDependenciesPlotDash"), cellDependenciesPlotDash(ns("celldepdash")))),
                 private(actionLink(inputId = ns("link_to_cellLineCoessentiality"), cellLineDependenciesPosTableDash(ns("celldeptabledash")))), 
                 private(actionLink(inputId = ns("link_to_cellLineDependenciesGraphDash"), cellDependenciesGraphDash(ns("celldepgraphdash"))))
               )
      ),
      ## INFO-----
      navbarMenu("INFO",
                 tabPanel(title = "Cell Line", value = "about_cell",
                          shinyjs::useShinyjs(),
                          fluidRow(
                            column(8, summary_var), #summary variable for alt descriptions
                            column(4, cellImage(ns("cell_image")))
                          ),
                          tags$br(),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("cell_summary_table_click"), cellSummaryTableTab(ns("celltab")))
                            )
                          ),
                          tags$br(),
                          #conditional for cell table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_tabcard"),
                                style = "padding-left:1%",
                                cellSummaryTable(ns("cell_table_tab")),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel(title = "Literature", value = "cell_literature",
                          #"co-citation graph?"
                          pubmedPlot(ns("pubmedcell")),  #Summary stats
                          pubmedTable(ns("pubmedcelltable"))
                 )
      ),
      ## EXPRESSION-----
      navbarMenu(title = "EXPRESSION",
                 tabPanel("Gene", value = "expression_cell", 
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            private(cellGeneExpressionPlot(ns("cellLine_gene"))),
                            private_msg()
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              private(actionLink(inputId = ns("gene_exp_table_click"), cellGeneExpressionTableTab(ns("gene_exp_table_tab")))),
                              private(actionLink(inputId = ns("protein_exp_plot_click"), cellProteinExpressionPlotTab(ns("protein_exp_plot_tab")))),
                              private(actionLink(inputId = ns("protein_exp_table_click"), cellProteinExpressionTableTab(ns("protein_exp_table_tab")))),
                              private(actionLink(inputId = ns("gene_protein_plot_click"), cellGeneProteinPlotTab(ns("gene_protein_plot_tab"))))
                            )
                          ),
                          tags$br(),
                          #conditional for cellGeneExpressionTableTab
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("gene_exp_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellGeneExpressionTable(ns("cellLine_gene_table"))),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for cellProteinExpressionPlot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("protein_exp_plot_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellProteinExpressionPlot(ns("cellLine_protein"))),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for cellProteinExpressionTable
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("protein_exp_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellProteinExpressionTable(ns("cellLine_protein_table"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for cellGeneProteinPlot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("gene_protein_plot_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellLineGeneProteinPlot(ns("cellLine_geneprotein"))),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel("Co-expression", value = "cell_coexpression",
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            private_msg(),
                            private(cellCoexpressionPlot(ns("cell_coexpression_plot"))) # MAIN PLOT
                                   ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              private(actionLink(inputId = ns("cell_exp_pos_table_click"), cellLineExpressionPosTableTab(ns("cell_exp_pos_table_tab")))), #pos table
                              private(actionLink(inputId = ns("cell_exp_neg_table_click"), cellLineExpressionNegTableTab(ns("cell_exp_neg_table_tab")))), #neg table
                              private(actionLink(inputId = ns("cell_metadata_exp_click"), cellLineMetadataPlotTab(ns("cell_metadata_exp_tab"))))
                            )
                          ),
                          tags$br(),
                          #conditional for pos table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_exp_pos_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(similarExpCellsTable(ns("cell_sim_exp"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for neg table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_exp_neg_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(dissimilarExpCellsTable(ns("cell_diss_exp"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for metadata plot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_metadata_exp_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellMetadataPlot(ns("cell_metadata_exp_plot"))),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel("Graph", value = "cell_expression_graph", 
                          private_msg(),
                          private(geneNetworkGraphExp(ns("cell_exprs_graph")))
                 ),
                 tabPanel("Differential Pathway Expression", value = "cell_functional_plot", 
                          private_msg(),
                          private(cellFunctionalPlot(ns("cell_fun_plot")))
                 )
      ),
      # COMPOUNDS-----
      navbarMenu(title = "COMPOUNDS",
                 tabPanel("Metabolites", value = "metabolites",
                          private_msg(),
                          private(metabolitesTable(ns("cell_metabolites")))
                 ),
                 tabPanel("Drugs", value = "drugs", 
                          private_msg(),
                          private(cellDrugsTable(ns("cell_drugs")))
                 )
      ),
      ## DEPENDENCIES (action)-----
      navbarMenu(title = "DEPENDENCIES",
                 tabPanel("Plots", value = "cell_dependencies_plots",
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            private(cellDependenciesPlot(ns("dep"))),
                            private_msg()
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              private(actionLink(inputId = ns("cell_dep_bar_click"), cellDependenciesBarPlotTab(ns("cell_bar_dep")))), #barplot
                              private(actionLink(inputId = ns("cell_dep_density_click"), cellDependenciesDensityPlotTab(ns("cell_dep_density_tab")))), #density
                              private(actionLink(inputId = ns("cell_dep_table_click"), cellLineDependenciesTableTab(ns("cell_dep_table_tab")))), #table
                              private(actionLink(inputId = ns("cell_expdep_click"), expdepPlotTab(ns("cell_expdep_plot_tab")))) #expdep
                            )
                          ),
                          tags$br(),
                          #conditional for barplot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_barplot_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellDependenciesBarPlot(ns("dep_barplot"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for density
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_density_tabcard"),
                                style = "padding-left:1%",
                                private(cellDependenciesDensityPlot(ns("dep"))),
                                private_msg(),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_dep_table_tabcard"),
                                style = "padding-left:1%",
                                private(cellLineDependenciesTable(ns("dep"))),
                                private_msg(),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for expdep
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_expdep_plot_tabcard"),
                                style = "padding-left:1%",
                                private(expdepPlot(ns("expdep_plot"))),
                                private_msg(),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel("Co-essentiality", value = "cell_dependencies_correlation",
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            private_msg(),
                            private(cellCoessentialityPlot(ns("cell_coessentiality_plot"))) # MAIN PLOT
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              private(actionLink(inputId = ns("cell_dep_pos_table_click"), cellLineDependenciesPosTableTab(ns("cell_dep_pos_table_tab")))), #pos table
                              private(actionLink(inputId = ns("cell_dep_neg_table_click"), cellLineDependenciesNegTableTab(ns("cell_dep_neg_table_tab")))), #neg table
                              private(actionLink(inputId = ns("cell_metadata_dep_click"), cellLineMetadataPlotTab(ns("cell_metadata_dep_tab"))))
                            )
                          ),
                          tags$br(),
                          #conditional for pos table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_dep_pos_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(similarCellsTable(ns("cell_sim"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for neg table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_dep_neg_table_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(dissimilarCellsTable(ns("cell_diss"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for metadata plot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("cell_metadata_dep_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(cellMetadataPlot(ns("cell_metadata_dep_plot"))),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel("Graph", value = "cell_dependencies_graph", 
                          private_msg(),
                          private(geneNetworkGraph(ns("cell_deps_graph")))
                 )
      ),
      ## DOWNLOADS-----
      downloadTab(ns("download")), #pull download tab from module
      
      ## SEARCH-----
      formContent=querySearchInput(ns("search"))
    )
  )
}

cellPageServer <- function(id, subtype) {
  moduleServer(
    id,
    function(input, output, session) {
      type <- "cell"
      if(subtype == "cell") {
        data <- reactive({
          cell_lines <- getQueryString()$query
          list(
            type=type,
            subtype=subtype,
            query=cell_lines,
            content=cell_lines
          )
        })
        title_var <- cellTitleServer("title", data)
        summary_var <- summaryListTextServer("summary", data)
      } else if(subtype == "lineage") {
        data <- reactive({
          lineage_str <- getQueryString()$query
          if (!is.null(lineage_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage == lineage_str)
            list(
              type=type,
              subtype=subtype,
              query=lineage_str,
              content=expression_name_row$cell_line
            )
          }
        })
        title_var <- lineageTitleServer("title", data) 
        summary_var <- lineageSummaryTextServer("summary", data)
      } else if(subtype == "lineage_subtype") {
        data <- reactive({
          lineage_subtype_str <- getQueryString()$query
          if (!is.null(lineage_subtype_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage_subtype == lineage_subtype_str)
            list(
              type=type,
              subtype=subtype,
              query=lineage_subtype_str,
              content=expression_name_row$cell_line
            )
          }
        })
        title_var <- lineageSubtypeTitleServer("title", data) 
        summary_var <- lineageSubtypeSummaryTextServer("summary", data)
        
      } else if(subtype == "cell_list") {
        data <- reactive({
          custom_cell_list <- getQueryString()$query
          cell_lines <- c(str_split(custom_cell_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            subtype=subtype,
            query=custom_cell_list,
            content=cell_lines
          )
        })
        title_var <- summaryListTitleServer("title", data) 
        summary_var <- summaryListTextServer("summary", data)
      } else {
        stop("fix your summary argument")
      }
      
      ## SEARCH SERVER-----
      querySearchServer("search")
      
      ## DASHBOARD SERVER-----
      #name pulls from if/else variables above
      #cell image plot
      observeEvent(input$link_tocellImagePlotDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "about_cell")
      })
      cellImageDashServer("cellimagedash", data)
      
      #pubmed gene plot
      observeEvent(input$link_to_pubmedCellPlotDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_literature")
      })
      pubmedPlotDashServer("pubmedcelldash", data)
      
      #dep plot
      observeEvent(input$link_to_cellLineDependenciesPlotDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_dependencies_plots")
      })
      private({cellDependenciesPlotDashServer("celldepdash", data)})
      
      # cell line dep table
      observeEvent(input$link_to_cellLineCoessentiality, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_dependencies_correlation")
      })
      private({cellLineDependenciesPosTableDashServer("celldeptabledash", data)})
      
      # cell line exp table
      observeEvent(input$link_to_cellLineCoexpression, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_coexpression")
      })
      private({cellLineExpressionPosTableDashServer("cellexptabledash", data)})
      
      # gene expression rug plot graph
      observeEvent(input$link_to_cellExpressionPlotDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "expression_cell")
      })

      private({cellExpressionPlotDashServer("cellexpressiondash", data)})

      ## INFO SERVER-----
      # Cell Image
      cellImageServer("cell_image", data)
      
      # Cell summary table
      observeEvent(input$cell_summary_table_click, { #store click
        shinyjs::show("cell_tabcard")
      })
      #serves the card for the table
      cellSummaryTableTabServer("celltab", data)
      #serves the table
      cellSummaryTableServer("cell_table_tab", data)
      
      # Literature
      pubmedPlotServer("pubmedcell", data)
      pubmedTableServer("pubmedcelltable", data)
      
      ## EXPRESSION SERVER-----
      # GENE
      private({cellGeneExpressionPlotServer("cellLine_gene", data)})
      # CONDITIONAL cellGeneExpressionTable
      observeEvent(input$gene_exp_table_click, { #store click
        shinyjs::show("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      private({cellGeneExpressionTableTabServer("gene_exp_table_tab", data)})
      #serves the data plots
      private({cellGeneExpressionTableServer("cellLine_gene_table", data)})
      
      # CONDITIONAL cellProteinExpressionPlot
      observeEvent(input$protein_exp_plot_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::show("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      private({cellProteinExpressionPlotTabServer("protein_exp_plot_tab", data)})
      #serves the data plots
      private({cellProteinExpressionPlotServer("cellLine_protein", data)})
      
      # CONDITIONAL cellProteinExpressionTable
      observeEvent(input$protein_exp_table_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::show("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      private({cellProteinExpressionTableTabServer("protein_exp_table_tab", data)})
      #serves the data plots
      private({cellProteinExpressionTableServer("cellLine_protein_table", data)})
      
      # CONDITIONAL cellGeneProteinPlotTab
      observeEvent(input$gene_protein_plot_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::show("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      private({cellGeneProteinPlotTabServer("gene_protein_plot_tab", data)})
      #serves the data plots
      private({cellLineGeneProteinPlotServer("cellLine_geneprotein", data)})
      
      ## COEXPRESSION
      private({cellCoexpressionPlotServer("cell_coexpression_plot", data)})
      
      # CONDITIONAL pos table
      observeEvent(input$cell_exp_pos_table_click, { #store click
        shinyjs::show("cell_exp_pos_table_tabcard")
        shinyjs::hide("cell_exp_neg_table_tabcard")
        shinyjs::hide("cell_metadata_exp_tabcard")
      })
      #serves the card for the image
      private({cellLineExpressionPosTableTabServer("cell_exp_pos_table_tab", data)})
      #serves the table
      private({similarExpCellsTableServer("cell_sim_exp", data)})
      
      # CONDITIONAL neg table
      observeEvent(input$cell_exp_neg_table_click, { #store click
        shinyjs::hide("cell_exp_pos_table_tabcard")
        shinyjs::show("cell_exp_neg_table_tabcard")
        shinyjs::hide("cell_metadata_exp_tabcard")
      })
      #serves the card for the image
      private({cellLineExpressionNegTableTabServer("cell_exp_neg_table_tab", data)})
      #serves the table
      private({dissimilarExpCellsTableServer("cell_diss_exp", data)})
      
      # CONDITIONAL metadata plot
      observeEvent(input$cell_metadata_exp_click, { #store click
        shinyjs::hide("cell_exp_pos_table_tabcard")
        shinyjs::hide("cell_exp_neg_table_tabcard")
        shinyjs::show("cell_metadata_exp_tabcard")
      })
      # serves the card
      private({cellLineMetadataPlotTabServer("cell_metadata_exp_tab", data, "expression")})
      # serves the plot
      private({cellMetadataPlotServer("cell_metadata_exp_plot", data, "expression")})
      
      # Graph
      observeEvent(input$link_to_cellLineExpressionGraphDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_expression_graph")
      })
      
      private({cellExpressionGraphDashServer("cellexpgraphdash", data)})
      private({geneNetworkGraphExpServer("cell_exprs_graph", data)})
      
      # FUNCTIONAL PLOT
      observeEvent(input$link_to_cellLineFunctionalPlotDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_functional_plot")
      })
      
      private({cellFunctionalPlotDashServer("cellfunplotdash", data)})
      private({cellFunctionalPlotServer("cell_fun_plot", data)})
      
      # COMPOUNDS SERVER-----
      ## cards
      observeEvent(input$link_to_cellLineDrugTable, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "drugs")
      })
      observeEvent(input$link_to_cellLineMetaboliteTable, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "metabolites")
      })
      
      private({cellMetabolitesTableDashServer("cellmetabolitetabledash", data)})
      private({cellDrugsTableDashServer("celldrugtabledash", data)})
      
      ## tables
      private({metabolitesTableServer("cell_metabolites", data)})
      
      # DEPENDENCIES SERVER----
      observeEvent(input$cell_dep_bar_click, { #store click
        shinyjs::show("cell_barplot_tabcard")
        shinyjs::hide("cell_dep_table_tabcard")
        shinyjs::hide("cell_density_tabcard")
        shinyjs::hide("cell_expdep_plot_tabcard")
      })
      
      observeEvent(input$cell_dep_density_click, { #store click
        shinyjs::hide("cell_barplot_tabcard")
        shinyjs::show("cell_density_tabcard")
        shinyjs::hide("cell_dep_table_tabcard")
        shinyjs::hide("cell_expdep_plot_tabcard")
      })
  
      observeEvent(input$cell_dep_table_click, { #store click
        shinyjs::hide("cell_barplot_tabcard")
        shinyjs::show("cell_dep_table_tabcard")
        shinyjs::hide("cell_density_tabcard")
        shinyjs::hide("cell_expdep_plot_tabcard")
      })
      
      observeEvent(input$cell_expdep_click, { #store click
        shinyjs::hide("cell_barplot_tabcard")
        shinyjs::hide("cell_dep_table_tabcard")
        shinyjs::hide("cell_density_tabcard")
        shinyjs::show("cell_expdep_plot_tabcard")
      })
      
      #serves the card for the image
      private({cellDependenciesDensityPlotTabServer("cell_dep_density_tab", data)})
      private({cellLineDependenciesTableTabServer("cell_dep_table_tab", data)})
      private({cellDependenciesBarPlotTabServer("cell_bar_dep", data)})
      private({expdepPlotTabServer("cell_expdep_plot_tab", data)})
      
      #serves plots and tables
      private({cellDependenciesPlotServer("dep", data)})
      private({cellDependenciesBarPlotServer("dep_barplot", data)})
      private({cellDependenciesDensityPlotServer("dep", data)})
      private({cellLineDependenciesTableServer("dep", data)})
      private({expdepPlotServer("expdep_plot", data)})
      # no equivalent for lineage/sublineage
      
      #COESSENTIALITY PAGE
      private({cellCoessentialityPlotServer("cell_coessentiality_plot", data)})
      
      # CONDITIONAL pos table
      observeEvent(input$cell_dep_pos_table_click, { #store click
        shinyjs::show("cell_dep_pos_table_tabcard")
        shinyjs::hide("cell_dep_neg_table_tabcard")
        shinyjs::hide("cell_metadata_dep_tabcard")
      })
      #serves the card for the image
      private({cellLineDependenciesPosTableTabServer("cell_dep_pos_table_tab", data)})
      #serves the table
      private({similarCellsTableServer("cell_sim", data)})
      
      # CONDITIONAL neg table
      observeEvent(input$cell_dep_neg_table_click, { #store click
        shinyjs::hide("cell_dep_pos_table_tabcard")
        shinyjs::show("cell_dep_neg_table_tabcard")
        shinyjs::hide("cell_metadata_dep_tabcard")
      })
      #serves the card for the image
      private({cellLineDependenciesNegTableTabServer("cell_dep_neg_table_tab", data)})
      #serves the table
      private({dissimilarCellsTableServer("cell_diss", data)})
      
      # CONDITIONAL metadata plot
      observeEvent(input$cell_metadata_dep_click, { #store click
        shinyjs::hide("cell_dep_pos_table_tabcard")
        shinyjs::hide("cell_dep_neg_table_tabcard")
        shinyjs::show("cell_metadata_dep_tabcard")
      })
      # serves the card
      private({cellLineMetadataPlotTabServer("cell_metadata_dep_tab", data, "dependency")})
      # serves the plot
      private({cellMetadataPlotServer("cell_metadata_dep_plot", data, "dependency")})
      
      # Graph
      observeEvent(input$link_to_cellLineDependenciesGraphDash, {
        updateNavbarPage(session, inputId = "cellNavBar", selected = "cell_dependencies_graph")
      })
      
      private({cellDependenciesGraphDashServer("celldepgraphdash", data)})
      private({geneNetworkGraphServer("cell_deps_graph", data)})
      
      # # Drug Cor
      private({cellDrugsTableServer("cell_drugs", data)})
      
      # DOWNLOAD SERVER-----
      downloadTabServer("download", data, privateMode) #pull download tab from module
      
    }
  )
}