#this page is for the master gene query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page
genePage <- function (id, subtype) {
  ns <- NS(id)
  
  #this block is the logic to define the summary_var variable, to display the proper summary module
  if(subtype == "gene"){
    title_var <- geneTitle(ns("title_var"))
    gene_var <- summaryListText(ns("gene_var"))
    protein_summary <- proteinText(ns("protein_summary"))
    summary_table <- pathwayList(ns("gene_pathways"))
  } else if (subtype == "pathway"){
    title_var <- pathwayTitle(ns("title_var"))
    gene_var <- pathwayText(ns("gene_var"))
    protein_summary <- pathwayText(ns("protein_summary"))
    summary_table <- pathwayGeneList(ns("gene_pathways"))
  } else if (subtype == "gene_list") {
    title_var <- summaryListTitle(ns("title_var"))
    gene_var <- summaryListText(ns("gene_var"))
    protein_summary <- summaryListText(ns("protein_summary"))
    summary_table <- pathwayList(ns("gene_pathways"))
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      id=ns("geneNavBar"),
      ## DASHBOARD-----
      tabPanel("SUMMARY",
               title_var,
               cardLayout(
                 barcodeDash(ns("barcodedash")),
                 actionLink(inputId = ns("link_to_ideogramPlotDash"), ideogramPlotDash(ns("ideogramdash"))),
                 actionLink(inputId = ns("link_to_structurePlotDash"), structureDash(ns("structuredash"))),
                 actionLink(inputId = ns("link_to_pubmedGenePlotDash"), pubmedPlotDash(ns("pubmedgenedash"))),
                 actionLink(inputId = ns("link_to_cellAnatogramPlotDash"), cellAnatogramPlotDash(ns("cellanatogramdash"))),
                 actionLink(inputId = ns("link_to_cellExpressionPlotDash"), cellExpressionPlotDash(ns("cellexpressiondash"))),
                 actionLink(inputId = ns("link_to_tissueAnatogramPlotDash"), tissueAnatogramPlotDash(ns("tissueanatogramdash"))),
                 actionLink(inputId = ns("link_to_cellDependenciesPlotDash"), cellDependenciesPlotDash(ns("depdash"))),
                 actionLink(inputId = ns("link_to_cellDependenciesTableDash"), cellDependenciesTableDash(ns("deptabledash"))),
                 actionLink(inputId = ns("link_to_cellDependenciesGraphDash"), cellDependenciesGraphDash(ns("depgraphdash"))),
                 private(actionLink(inputId = ns("link_to_geneDrugsCorTableDash"), geneDrugsCorTableDash(ns("genedrugscortabledash"))))
               )
      ),
      ## INFO (person)-----
      navbarMenu("INFO",
                 tabPanel(title = "Gene", value = "about_gene",
                          shinyjs::useShinyjs(),
                          #summary
                          fluidRow(
                            column(8, gene_var), #summary variable for alt descriptions
                            column(4, ideogramPlot(ns("chromo")))
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("go_click"), geneGoTableTab(ns("gotab"))) 
                            )
                          ),
                          tags$br(),
                          #conditional for go table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("go_tabcard"),
                                style = "padding-left:1%",
                                #instead of modules here, module logic is above
                                summary_table,
                                tags$br()
                              )
                            )
                          )
                 ), #end tab panel
                 tabPanel(title = "Protein", value = "about_protein",
                          shinyjs::useShinyjs(),
                          #summary
                          fluidRow(
                            column(8, protein_summary),
                            column(4, proteinStructurePlot(ns("structure")))
                          ),
                          #"ADD: seq, blastP link, pfam, "
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("size_click"), sizePlotTab(ns("sizetab"))), 
                              actionLink(inputId = ns("sequence_click"), sequencePlotTab(ns("sequencetab"))), 
                              actionLink(inputId = ns("signature_click"), signaturePlotTab(ns("signaturetab"))), 
                              actionLink(inputId = ns("structure_click"), structurePlotTab(ns("structuretab")))
                            )
                          ),
                          tags$br(),
                          #conditional for size
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("size_tabcard"),
                                style = "padding-left:1%",
                                proteinSizeText(ns("protein_size_text")),
                                proteinSizePlot(ns("protein_size")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for sequence
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("sequence_tabcard"),
                                style = "padding-left:1%",
                                proteinSeqText(ns("protein_seq_text")),
                                proteinSeq(ns("protein_seq")),
                                proteinDomainPlot(ns("protein_domain_plot")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for signature
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("signature_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(proteinSignatureText(ns("protein_signature_text"))),
                                private(radialPlot(ns("radial_plot"))),
                                private(AABarPlot(ns("aa_bar_plot"))),
                                private(tags$hr()),
                                private(UMAPPlot(ns("umap_plot"))),
                                private(ClusterradialPlot(ns("cluster_radial_plot"))),
                                private(ClusterAABarPlot(ns("cluster_aa_bar_plot"))),
                                private(proteinClusterTable(ns("prot_clust_table"))),
                                private(ClusterEnrichmentPlot(ns("cluster_enrichment_plot"))),
                                private(proteinClusterEnrichmentTable(ns("prot_clust_enrich_table"))),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for structure
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("structure_tabcard"),
                                style = "padding-left:1%",
                                private_msg(),
                                private(proteinStructureText(ns("protein_structure_text"))),
                                private(proteinStructurePlot3d(ns("protein_structure_plot3d"))),
                                tags$br()
                              )
                            )
                          )
                 ), #end tab panel
                 tabPanel(title = "Literature", value = "about_literature",
                          #co-citation graph? 
                          pubmedPlot(ns("pubmed")),  #Summary stats
                          pubmedTable(ns("pubmed"))
                 )
      ),
      ## EXPRESSION (place)-----
      navbarMenu(title = "EXPRESSION",
                 tabPanel("Sub-cellular", value = "expression_sub",
                          cellAnatogramPlot(ns("exp")),
                          fluidRow(actionLink(inputId = ns("anato_facet_click"), "View detailed cell anatograms below")), # conditional panel for raw data
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("anato_facet_click"), "'] != 0"), 
                                           cellAnatogramFacetPlot(ns("exp"))),
                          tags$br(),
                          cellAnatogramTable(ns("exp"))
                 ), 
                 tabPanel("Cell Line", value = "expression_cell", 
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            cellGeneExpressionPlot(ns("cell_gene"))
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("gene_exp_table_click"), cellGeneExpressionTableTab(ns("gene_exp_table_tab"))), 
                              actionLink(inputId = ns("protein_exp_plot_click"), cellProteinExpressionPlotTab(ns("protein_exp_plot_tab"))), 
                              actionLink(inputId = ns("protein_exp_table_click"), cellProteinExpressionTableTab(ns("protein_exp_table_tab"))), 
                              actionLink(inputId = ns("gene_protein_plot_click"), cellGeneProteinPlotTab(ns("gene_protein_plot_tab")))
                            )
                          ),
                          tags$br(),
                          #conditional for cellGeneExpressionTableTab
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("gene_exp_table_tabcard"),
                                style = "padding-left:1%",
                                # cellGeneExpressionTableText(ns("cell_gene_table_text")),
                                cellGeneExpressionTable(ns("cell_gene_table")),
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
                                # cellProteinExpressionPlotText(ns("cell_protein_text")), 
                                cellProteinExpressionPlot(ns("cell_protein")),
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
                                cellProteinExpressionTable(ns("cell_protein_table")),
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
                                #plots go here
                                cellGeneProteinPlot(ns("cell_geneprotein")),
                                tags$br()
                              )
                            )
                          )
                 ), #end tab panel
                 tabPanel("Tissue", value = "expression_tissue", 
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            tissueTitle(ns("tissue_title")), #from shiny_text
                          ),
                          fluidRow(
                            column(6, maleAnatogramPlot(ns("male_anatogram"))),
                            column(6, femaleAnatogramPlot(ns("female_anatogram")))
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("tissue_plot_click"), tissuePlotTab(ns("tissue_plot_tab"))), 
                              actionLink(inputId = ns("tissue_table_click"), tissueTableTab(ns("tissue_table_tab")))
                            )
                          ),
                          tags$br(),
                          #conditional for tissuePlotTab
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("tissue_plot_tabcard"),
                                style = "padding-left:1%",
                                tissuePlotText(ns("tissue_plot_text")),
                                tissuePlot(ns("tissue_plot")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for tissueTableTab
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("tissue_table_tabcard"),
                                style = "padding-left:1%",
                                tissueTableText(ns("tissue_table_text")), 
                                tissueTable(ns("tissue_table")),
                                tags$br()
                              )
                            )
                          )
                 )#, #end tab panel
                 # tabPanel("Coexpression", value = "expression_co",
                 #          "Integrated co-expression analyses")
      ),
      # COMPOUNDS (thing)-----
      navbarMenu(title = "COMPOUNDS",
                 tabPanel("Drugs", value = "drugs",
                          private(geneDrugsTable(ns("gene_drugs"))), 
                          private_msg()),
                 tabPanel("Metabolites", value = "metabolites",
                          private(metabolitesTable(ns("gene_metabolites"))), 
                          private_msg()),
                 tabPanel("Graph", value = "gene_bipartite_graph",
                          private(geneBipartiteGraph(ns("gene_bipartite_graph"))), 
                          private_msg())
      ),
      ## DEPENDENCIES (action)-----
      navbarMenu(title = "DEPENDENCIES",
                 tabPanel("Plots", value = "dependencies_plots", 
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(
                            cellDependenciesPlot(ns("dep"))
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("dep_bar_click"), cellDependenciesBarPlotTab(ns("bar_dep"))), #barplot
                              actionLink(inputId = ns("dep_density_click"), cellDependenciesDensityPlotTab(ns("dep_density_tab"))), #density
                              actionLink(inputId = ns("dep_lineage_click"), cellDepsLinPlotTab(ns("dep_lineage_tab"))), #lineage
                              actionLink(inputId = ns("dep_sublineage_click"), cellDepsSubLinPlotTab(ns("dep_sublineage_tab"))), #sublineage
                              actionLink(inputId = ns("dep_table_click"), cellDependenciesTableTab(ns("dep_table_tab"))), #table
                              actionLink(inputId = ns("expdep_click"), expdepPlotTab(ns("expdep_plot_tab"))) #expdep
                            )
                          ),
                          tags$br(),
                          # conditional for barplot
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("barplot_tabcard"),
                                style = "padding-left:1%",
                                cellDependenciesBarPlot(ns("dep_barplot")), 
                                tags$br()
                              )
                            )
                          ),
                          #conditional for density
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("density_tabcard"),
                                style = "padding-left:1%",
                                cellDependenciesDensityPlot(ns("dep_density")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for lineage
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("lineage_tabcard"),
                                style = "padding-left:1%",
                                # cellDepsLinPlotText(ns("dep_lineage_text")), 
                                cellDepsLinPlot(ns("dep_lineage")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for sublineage
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("sublineage_tabcard"),
                                style = "padding-left:1%",
                                # cellDepsSubLinPlotText(ns("dep_sublineage_text")), 
                                cellDepsSubLinPlot(ns("dep_sublineage")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("dep_table_tabcard"),
                                style = "padding-left:1%",
                                # cellDependenciesTableText(ns("dep_table_text")), 
                                cellDependenciesTable(ns("dep_table")),
                                tags$br()
                              )
                            )
                          ),
                          #conditional for expdep
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("expdep_plot_tabcard"),
                                style = "padding-left:1%",
                                expdepPlot(ns("expdep_plot")),
                                tags$br()
                              )
                            )
                          )
                 ),
                 tabPanel("Co-essentiality", value = "dependencies_co-essentiality",
                          shinyjs::useShinyjs(),
                          #summary plot
                          fluidRow(cellDependenciesCorrPlot(ns("corrplot"))
                          ),
                          tags$hr(),
                          # cards in a fluid row
                          fluidRow(
                            cardLayout(
                              actionLink(inputId = ns("dep_pos_table_click"), cellDependenciesPosTableTab(ns("dep_pos_table_tab"))), #pos table
                              actionLink(inputId = ns("dep_pos_pathways_click"), cellDependenciesPosPathwayTableTab(ns("dep_pos_pathways_tab"))), #pos pathways
                              actionLink(inputId = ns("dep_neg_table_click"), cellDependenciesNegTableTab(ns("dep_neg_table_tab"))), #neg table
                              actionLink(inputId = ns("dep_neg_pathways_click"), cellDependenciesNegPathwayTableTab(ns("dep_neg_pathways_tab"))) #neg pathways
                            )
                          ),
                          tags$br(),
                          #conditional for pos table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("dep_pos_table_tabcard"),
                                style = "padding-left:1%",
                                # similarGenesTableText(ns("sim_text")), 
                                similarGenesTable(ns("sim")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for pos pathways
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("dep_pos_pathways_tabcard"),
                                style = "padding-left:1%",
                                # similarPathwaysTableText(ns("sim_pathways_text")), 
                                similarPathwaysTable(ns("sim_pathways")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for neg table
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("dep_neg_table_tabcard"),
                                style = "padding-left:1%",
                                # dissimilarGenesTableText(ns("dsim_text")), 
                                dissimilarGenesTable(ns("dsim")),
                                tags$br()
                              )
                            )
                          ), 
                          #conditional for neg pathways
                          fluidRow(
                            shinyjs::hidden(
                              div(
                                id = ns("dep_neg_pathways_tabcard"),
                                style = "padding-left:1%",
                                # dissimilarPathwaysTableText(ns("dsim_pathways_text")), 
                                dissimilarPathwaysTable(ns("dsim_pathways")),
                                tags$br()
                              )
                            )
                          ) #end of fluid row
                 ),###end of tabpanel
                 tabPanel("Graph", value = "dependencies_graph", 
                          geneNetworkGraph(ns("graph"))
                 ),
                 tabPanel("Drugs", value = "dependencies_drugs", 
                          private(geneDrugsCorTable(ns("gene_drugs_cor"))), 
                          private_msg())
      ),
      ## DOWNLOADS-----
      downloadTab(ns("download")), #pull download tab from module
      
      ## SEARCH-----
      formContent=querySearchInput(ns("search"))
    )
  )
}

genePageServer <- function(id, subtype) {
  moduleServer(
    id,
    function(input, output, session) {
      type <- "gene"
      #this block is the logic to define the summary_var variable, to display the proper summary module
      if(subtype == "gene"){
        data <- reactive({
          gene_symbols <- getQueryString()$query
          list(
            type=type,
            subtype=subtype,
            query=gene_symbols,
            content=gene_symbols
          )
        })
        title_var <- geneTitleServer("title_var", data)
        gene_var <- summaryListTextServer("gene_var", data)
        protein_summary <- proteinTextServer("protein_summary", data)
        summary_table <- pathwayListServer("gene_pathways", data)
      } else if (subtype == "pathway"){
        data <- reactive({
          pathway_go <- getQueryString()$query
          pathway_row <- pathways %>%
            filter(go %in% pathway_go)
          list(
            type=type,
            subtype=subtype,
            query=pathway_go,
            content=pathway_row$data[[1]]$gene
            )
        })
        title_var <- pathwayTitleServer("title_var", data)
        gene_var <- pathwayTextServer("gene_var", data)
        protein_summary <- pathwayTextServer("protein_summary", data)
        summary_table <- pathwayGeneListServer("gene_pathways", data)
      } else if (subtype == "gene_list") {
        data <- reactive({
          custom_gene_list <- getQueryString()$query
          gene_symbols <- c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            subtype=subtype,
            query=custom_gene_list,
            content=gene_symbols
          )
        })
        title_var <- summaryListTitleServer("title_var", data)
        gene_var <- summaryListTextServer("gene_var", data)
        protein_summary <- summaryListTextServer("protein_summary", data)
        summary_table <- pathwayListServer("gene_pathways", data)
      } else {
        stop("fix your summary argument")
      }
      
      ## SEARCH SERVER-----
      querySearchServer("search")
      
      ## DASHBOARD SERVER-----
      
      #barcode plot
      #no observe event, because it's just a linkout
      barcodeDashServer("barcodedash", data)
      
      #ideogram plot
      observeEvent(input$link_to_ideogramPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "about_gene")
      })
      ideogramPlotDashServer("ideogramdash", data)
      
      #structure plot
      observeEvent(input$link_to_structurePlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "about_protein")
      })
      structureDashServer("structuredash", data)
      
      #pubmed gene plot
      observeEvent(input$link_to_pubmedGenePlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "about_literature")
      })
      pubmedPlotDashServer("pubmedgenedash", data)
      
      #cell anatogram plot
      observeEvent(input$link_to_cellAnatogramPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "expression_sub")
      })
      cellAnatogramPlotDashServer("cellanatogramdash", data)
      
      #cell expression plot
      observeEvent(input$link_to_cellExpressionPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "expression_cell")
      })
      cellExpressionPlotDashServer("cellexpressiondash", data)
      
      #tissue anatogram plot
      observeEvent(input$link_to_tissueAnatogramPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "expression_tissue")
      })
      tissueAnatogramPlotDashServer("tissueanatogramdash", data)
      
      #dep plot
      observeEvent(input$link_to_cellDependenciesPlotDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "dependencies_plots")
      })
      cellDependenciesPlotDashServer("depdash", data)
      
      #dep table
      observeEvent(input$link_to_cellDependenciesTableDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "dependencies_co-essentiality")
      })
      cellDependenciesTableDashServer("deptabledash", data)
      
      #dep graph
      observeEvent(input$link_to_cellDependenciesGraphDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "dependencies_graph")
      })
      cellDependenciesGraphDashServer("depgraphdash", data)
      
      #dep drugs
      #do i need a wrapper around this observent event too?
      observeEvent(input$link_to_geneDrugsCorTableDash, {
        updateNavbarPage(session, inputId = "geneNavBar", selected = "dependencies_drugs")
      })
      private(geneDrugsCorTableDashServer("genedrugscortabledash", data))
      
      ## INFO SERVER-----
      #name pulls from if/else variables above
      
      # GENE
      gene_var
      ideogramPlotServer("chromo", data)
      
      # CONDITIONAL GO PATHWAY
      observeEvent(input$go_click, { #store click
        shinyjs::show("go_tabcard")
        # shinyjs::hide("sequence_tabcard") #template for when we add more
      })
      #serves the card for the image
      geneGoTableTabServer("gotab", data)
      #serves the data plots
      geneGoTableTextServer("go_table_text", data)
      #***summary table for go terms is assigned above***
      
      # PROTEIN
      protein_summary
      proteinStructurePlotServer("structure", data)
      
      # CONDITIONAL SIZE
      observeEvent(input$size_click, { #store click
        shinyjs::show("size_tabcard")
        shinyjs::hide("sequence_tabcard")
        shinyjs::hide("signature_tabcard")
        shinyjs::hide("structure_tabcard")
      })
      #serves the card for the image
      sizePlotTabServer("sizetab", data)
      #serves the data plots
      proteinSizeTextServer("protein_size_text", data)
      proteinSizePlotServer("protein_size", data)
      
      # CONDITIONAL SEQUENCE
      observeEvent(input$sequence_click, { #store click
        shinyjs::hide("size_tabcard")
        shinyjs::show("sequence_tabcard")
        shinyjs::hide("signature_tabcard")
        shinyjs::hide("structure_tabcard")
      })
      #serves the card for the image
      sequencePlotTabServer("sequencetab", data)
      #serves the data plots
      proteinSeqTextServer("protein_seq_text", data)
      proteinSeqServer("protein_seq", data)
      proteinDomainPlotServer("protein_domain_plot", data)
      
      # CONDITIONAL SIGNATURE
      observeEvent(input$signature_click, { #store click
        shinyjs::hide("size_tabcard")
        shinyjs::hide("sequence_tabcard")
        shinyjs::show("signature_tabcard")
        shinyjs::hide("structure_tabcard")
      })
      #serves the card for the image
      signaturePlotTabServer("signaturetab", data)
      #serves the data plots
      private({proteinSignatureTextServer("protein_signature_text", data)})
      private({radialPlotServer("radial_plot", data)})
      private({AABarPlotServer("aa_bar_plot", data)})
      private({UMAPPlotServer("umap_plot", data)})
      private({ClusterradialPlotServer("cluster_radial_plot", data)})
      private({ClusterAABarPlotServer("cluster_aa_bar_plot", data)})
      private({proteinClusterTableServer("prot_clust_table", data)})
      private({ClusterEnrichmentPlotServer("cluster_enrichment_plot", data)})
      private({proteinClusterEnrichmentTableServer("prot_clust_enrich_table", data)})
      
      # CONDITIONAL STRUCTURE
      observeEvent(input$structure_click, { #store click
        shinyjs::hide("size_tabcard")
        shinyjs::hide("sequence_tabcard")
        shinyjs::hide("signature_tabcard")
        shinyjs::show("structure_tabcard")
      })
      #serves the card for the image
      structurePlotTabServer("structuretab", data)
      #serves the data plots
      private({proteinStructureTextServer("protein_structure_text", data)})
      private({proteinStructurePlot3dServer("protein_structure_plot3d", data)})
      
      # Literature
      pubmedPlotServer("pubmed", data)
      pubmedTableServer("pubmed", data)
      
      ## EXPRESSION SERVER-----
      #Subcell
      cellAnatogramPlotServer("exp", data)
      observeEvent(input$anato_facet_click, { #event to store the anatogram facet 'click'
      })
      cellAnatogramFacetPlotServer("exp", data)
      cellAnatogramTableServer("exp", data)
      
      # CELL LINE
      cellGeneExpressionPlotServer("cell_gene", data)
      
      # CONDITIONAL cellGeneExpressionTable
      observeEvent(input$gene_exp_table_click, { #store click
        shinyjs::show("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      cellGeneExpressionTableTabServer("gene_exp_table_tab", data)
      #serves the data plots
      cellGeneExpressionTableTextServer("cell_gene_table_text", data)
      cellGeneExpressionTableServer("cell_gene_table", data)
      
      # CONDITIONAL cellProteinExpressionPlot
      observeEvent(input$protein_exp_plot_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::show("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      cellProteinExpressionPlotTabServer("protein_exp_plot_tab", data)
      #serves the data plots
      cellProteinExpressionPlotTextServer("cell_protein_text", data)
      cellProteinExpressionPlotServer("cell_protein", data)
      
      # CONDITIONAL cellProteinExpressionTable
      observeEvent(input$protein_exp_table_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::show("protein_exp_table_tabcard")
        shinyjs::hide("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      cellProteinExpressionTableTabServer("protein_exp_table_tab", data)
      #serves the data plots
      cellProteinExpressionTableServer("cell_protein_table", data)
      
      # CONDITIONAL cellGeneProteinPlotTab
      observeEvent(input$gene_protein_plot_click, { #store click
        shinyjs::hide("gene_exp_table_tabcard")
        shinyjs::hide("protein_exp_plot_tabcard")
        shinyjs::hide("protein_exp_table_tabcard")
        shinyjs::show("gene_protein_plot_tabcard")
      })
      #serves the card for the image
      cellGeneProteinPlotTabServer("gene_protein_plot_tab", data)
      #serves the data plots
      cellGeneProteinPlotServer("cell_geneprotein", data)
      
      # HUMAN TISSUE
      tissueTitleServer("tissue_title", data)
      maleAnatogramPlotServer("male_anatogram", data)
      femaleAnatogramPlotServer("female_anatogram", data)
      
      # CONDITIONAL tissuePlot
      observeEvent(input$tissue_plot_click, { #store click
        shinyjs::show("tissue_plot_tabcard")
        shinyjs::hide("tissue_table_tabcard")
      })
      #serves the card for the image
      tissuePlotTabServer("tissue_plot_tab", data)
      #serves the data plots
      tissuePlotTextServer("tissue_plot_text", data)
      tissuePlotServer("tissue_plot", data)
      
      # CONDITIONAL tissueTable
      observeEvent(input$tissue_table_click, { #store click
        shinyjs::hide("tissue_plot_tabcard")
        shinyjs::show("tissue_table_tabcard")
      })
      #serves the card for the image
      tissueTableTabServer("tissue_table_tab", data)
      #serves the data plots
      tissueTableTextServer("tissue_table_text", data)
      tissueTableServer("tissue_table", data)
      
      #COEXPRESSION
      #still thinking
      
      # COMPOUNDS SERVER-----
      private({geneDrugsTableServer("gene_drugs", data)})
      private({metabolitesTableServer("gene_metabolites", data)})
      private({geneBipartiteGraphServer("gene_bipartite_graph", data)})
      
      # DEPENDENCIES SERVER----
      cellDependenciesPlotServer("dep", data)
      
      # CONDITIONAL Barplot
      # conditionalPanel approach
      # output$multiple4barplot <- reactive({
      #   if(!is.null(data()$query)) {
      #     gene_len <- length(data()$query)
      #     return(gene_len)
      #   }
      # })
      # outputOptions(output, "multiple4barplot", suspendWhenHidden = FALSE)
      
      # shinyjs approach
      # if(!is.null(data()$query)) {
      # gene_len <- length(data()$query)
      # return(gene_len)
      # }
      # if(gene_len > 1) {
      # multyquery <- TRUE
      # } 
      # else {
      # multyquery <- FALSE
      # }
      
      # if(multyquery) {
      # shinyjs::show("conditional_w_barplot")
      # shinyjs::hide("conditional_wo_barplot")
      # } else {
      # shinyjs::hide("conditional_w_barplot")
      # shinyjs::show("conditional_wo_barplot")
      # }
      
      observeEvent(input$dep_bar_click, { #store click
        shinyjs::show("barplot_tabcard")
        shinyjs::hide("density_tabcard")
        shinyjs::hide("lineage_tabcard")
        shinyjs::hide("sublineage_tabcard")
        shinyjs::hide("dep_table_tabcard")
        shinyjs::hide("expdep_plot_tabcard")
      })
      
      #serves the card for the image
      cellDependenciesBarPlotTabServer("bar_dep", data)
      #serves the data plots
      cellDependenciesBarPlotServer("dep_barplot", data)
      
      # CONDITIONAL density
      observeEvent(input$dep_density_click, { #store click
        shinyjs::hide("barplot_tabcard")
        shinyjs::show("density_tabcard")
        shinyjs::hide("lineage_tabcard")
        shinyjs::hide("sublineage_tabcard")
        shinyjs::hide("dep_table_tabcard")
        shinyjs::hide("expdep_plot_tabcard")
      })
      #serves the card for the image
      cellDependenciesDensityPlotTabServer("dep_density_tab", data)
      #serves the data plots
      cellDependenciesDensityPlotServer("dep_density", data)
      
      # CONDITIONAL lineage
      observeEvent(input$dep_lineage_click, { #store click
        shinyjs::hide("barplot_tabcard")
        shinyjs::hide("density_tabcard")
        shinyjs::show("lineage_tabcard")
        shinyjs::hide("sublineage_tabcard")
        shinyjs::hide("dep_table_tabcard")
        shinyjs::hide("expdep_plot_tabcard")
      })
      #serves the card for the image
      cellDepsLinPlotTabServer("dep_lineage_tab", data)
      #serves the data plots
      cellDepsLinPlotTextServer("dep_lineage_text", data)
      cellDepsLinPlotServer("dep_lineage", data)
      
      # CONDITIONAL sublineage
      observeEvent(input$dep_sublineage_click, { #store click
        shinyjs::hide("barplot_tabcard")
        shinyjs::hide("density_tabcard")
        shinyjs::hide("lineage_tabcard")
        shinyjs::show("sublineage_tabcard")
        shinyjs::hide("dep_table_tabcard")
        shinyjs::hide("expdep_plot_tabcard")
      })
      #serves the card for the image
      cellDepsSubLinPlotTabServer("dep_sublineage_tab", data)
      #serves the data plots
      cellDepsSubLinPlotTextServer("dep_sublineage_text", data)
      cellDepsSubLinPlotServer("dep_sublineage", data)
      
      # CONDITIONAL table
      observeEvent(input$dep_table_click, { #store click
        shinyjs::hide("barplot_tabcard")
        shinyjs::hide("density_tabcard")
        shinyjs::hide("lineage_tabcard")
        shinyjs::hide("sublineage_tabcard")
        shinyjs::show("dep_table_tabcard")
        shinyjs::hide("expdep_plot_tabcard")
      })
      #serves the card for the image
      cellDependenciesTableTabServer("dep_table_tab", data)
      #serves the data table
      cellDependenciesTableTextServer("dep_table_text", data)
      cellDependenciesTableServer("dep_table", data)    
      
      # CONDITIONAL expdep
      observeEvent(input$expdep_click, { #store click
        shinyjs::hide("barplot_tabcard")
        shinyjs::hide("density_tabcard")
        shinyjs::hide("lineage_tabcard")
        shinyjs::hide("sublineage_tabcard")
        shinyjs::hide("dep_table_tabcard")
        shinyjs::show("expdep_plot_tabcard")
      })
      #serves the card for the image
      expdepPlotTabServer("expdep_plot_tab", data)
      #serves the data plots
      expdepPlotServer("expdep_plot", data)
      
      #COESSENTIALITY PAGE
      cellDependenciesCorrPlotServer("corrplot", data)
      
      # CONDITIONAL pos table
      observeEvent(input$dep_pos_table_click, { #store click
        shinyjs::show("dep_pos_table_tabcard")
        shinyjs::hide("dep_pos_pathways_tabcard")
        shinyjs::hide("dep_neg_table_tabcard")
        shinyjs::hide("dep_neg_pathways_tabcard")
        shinyjs::hide("dep_corr_tabcard")
      })
      #serves the card for the image
      cellDependenciesPosTableTabServer("dep_pos_table_tab", data)
      #serves the data plots
      similarGenesTableTextServer("sim_text", data)
      similarGenesTableServer("sim", data)
      
      # CONDITIONAL pos pathways
      observeEvent(input$dep_pos_pathways_click, { #store click
        shinyjs::hide("dep_pos_table_tabcard")
        shinyjs::show("dep_pos_pathways_tabcard")
        shinyjs::hide("dep_neg_table_tabcard")
        shinyjs::hide("dep_neg_pathways_tabcard")
        shinyjs::hide("dep_corr_tabcard")
      })
      #serves the card for the image
      cellDependenciesPosPathwayTableTabServer("dep_pos_pathways_tab", data)
      #serves the data plots
      similarPathwaysTableTextServer("sim_pathways_text", data)
      similarPathwaysTableServer("sim_pathways", data)
      
      # CONDITIONAL neg table
      observeEvent(input$dep_neg_table_click, { #store click
        shinyjs::hide("dep_pos_table_tabcard")
        shinyjs::hide("dep_pos_pathways_tabcard")
        shinyjs::show("dep_neg_table_tabcard")
        shinyjs::hide("dep_neg_pathways_tabcard")
        shinyjs::hide("dep_corr_tabcard")
      })
      #serves the card for the image
      cellDependenciesNegTableTabServer("dep_neg_table_tab", data)
      #serves the data plots
      dissimilarGenesTableTextServer("dsim_text", data)
      dissimilarGenesTableServer("dsim", data)
      
      # CONDITIONAL neg pathways
      observeEvent(input$dep_neg_pathways_click, { #store click
        shinyjs::hide("dep_pos_table_tabcard")
        shinyjs::hide("dep_pos_pathways_tabcard")
        shinyjs::hide("dep_neg_table_tabcard")
        shinyjs::show("dep_neg_pathways_tabcard")
        shinyjs::hide("dep_corr_tabcard")
      })
      #serves the card for the image
      cellDependenciesNegPathwayTableTabServer("dep_neg_pathways_tab", data)
      #serves the data plots
      dissimilarPathwaysTableTextServer("dsim_pathways_text", data)
      dissimilarPathwaysTableServer("dsim_pathways", data)
      
      #Graph
      geneNetworkGraphServer("graph", data)
      
      
      # Drug Cor
      private({geneDrugsCorTableServer("gene_drugs_cor", data)})
      
      # DOWNLOAD SERVER-----
      downloadTabServer("download", data, privateMode) #pull download tab from module
    }
  )
}

