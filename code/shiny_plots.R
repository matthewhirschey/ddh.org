# GENE QUERY-----
## Ideogram plot--------------------------------------------------------
ideogramPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_ideogram_plot")))),
    uiOutput(outputId = ns("conditional_ideogramplot")),
    tags$br(),
    fluidRow(ddh::make_legend("make_ideogram"))
  )
}

ideogramPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_ideogramplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "ideogram_plot"
        fieldnameMultiple <- "ideogram_multiple_plot"
        ImgFun <- "make_ideogram"
        query_type <- "gene"
        image_type <- "plot"
        SpinnerColor <- query_type
        PlotSizes <- paste0(plot_size_finder(ImgFun), "px")
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$ideogram_plot <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene, "No data found for this gene."))
        list(src= load_image(input = data(), 
                             fun_name = "make_ideogram", 
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_ideogram")$plot_width,
             height = plot_size_finder("make_ideogram")$plot_height
        )
      }, deleteFile = FALSE)
      output$text_ideogram_plot <- renderText({paste0("Chromosome location for ", str_c(data()$content, collapse = ", "))})
      output$ideogram_multiple_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% gene_location$approved_symbol, "No data found for this gene."))
        make_ideogram(location_data = gene_location, chromosome_data = chromosome, input = data())
      })      
    }
  )
}

# PROTEIN SIZE PLOT --------------------------------------------------------
proteinSizePlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_protein_size_plot")))),
    uiOutput(outputId = ns("conditional_proteinsizeplot")),
    tags$br(),
    fluidRow(ddh::make_legend("make_proteinsize"))
  )
}

proteinSizePlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_proteinsizeplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "proteinsizeplot"
        fieldnameMultiple <- "proteinsizeplot_multiple"
        ImgFun <- "make_proteinsize"
        query_type <- "gene"
        image_type <- "plot"
        SpinnerColor <- "protein"
        PlotSizes <- paste0(plot_size_finder(ImgFun), "px")
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        }
      })
      output$proteinsizeplot <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene, "No data found for this gene."))
        list(src= load_image(input = data(), 
                             fun_name = "make_proteinsize", 
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_proteinsize")$plot_width,
             height = plot_size_finder("make_proteinsize")$plot_height
        )
      }, deleteFile = FALSE)
      output$text_protein_size_plot <- renderText({paste0("Protein size for ", str_c(data()$content, collapse = ", "))})
      output$proteinsizeplot_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "No protein mass found"))
        make_proteinsize(input = data())
      },
      height = function() length(data()$content) * 120 + 60)
    }
  )
}

# PROTEIN STRUCTURE PLOT --------------------------------------------------------
proteinStructurePlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(imageOutput(outputId = ns("protein_structure"), height = "auto") %>% 
               withSpinnerColor(plot_type = "protein") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_structure"))
  )
}

proteinStructurePlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_protein_structure <- renderText({paste0(str_c(data()$content[[1]], " Predicted Structure", collapse = ", "))})
      output$protein_structure <- renderImage({
        shiny::validate(
        need(data()$content %in% achilles_long$gene, "No structure found for this protein")
      )
      list(src= load_image(input = data(), fun_name = "make_structure"),
           width="100%")
    }, deleteFile = FALSE)
    }
  )
}

# PROTEIN 3D STRUCTURE PLOT --------------------------------------------------------
proteinStructurePlot3d <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_protein_structure3d")))),
    tags$br(),
    fluidRow(
      column(8,
             r3dmol::r3dmolOutput(outputId = ns("protein_structure3D"),
                                  height = "600px") %>% 
               withSpinnerColor(plot_type = "protein") #see shiny_helper.R
      ),
      column(4,
             uiOutput(ns("gene_ui")),
             h4("Visualization parameters:"),
             shinyWidgets::prettySwitch(inputId = ns("color3dstructure"), 
                                        "Color structure", value = FALSE),
             
             shinyWidgets::prettySwitch(inputId = ns("ribbon3dstructure"), 
                                        "Ribbon structure", value = FALSE),
             
             shinyWidgets::prettySwitch(inputId = ns("selection3dstructure"), 
                                        "Select a region", value = FALSE),
             
             conditionalPanel(condition = paste0("input['", ns("selection3dstructure"), "']"),
                              textInput(ns("residue"), 
                                        label = "Highlight residues:",
                                        placeholder = "e.g., 1:10"),
                              textInput(ns("chain"), 
                                        label = "Chain:",
                                        placeholder = "e.g., A")
                              ),
             actionButton(inputId = ns("update3d"), 
                          label = "Update plot")
      )
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_structure3d"),
             actionLink(inputId = ns("pdb_table_click"), " View more PDB models for my protein/s")),
    ## TABLE
    conditionalPanel(condition = paste0("input['", ns("pdb_table_click"), "'] != 0"),
                     fluidRow(h4(textOutput(ns("title_structure3d_table")))),
                     DT::dataTableOutput(outputId = ns("structure3d_table"))
                     )
  )
}

proteinStructurePlot3dServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$gene_ui <- renderUI({
        
        gene_options <- make_structure3d_table(input = data()) %>%
          pull(gene_name) %>%
          unique()
        
        selectizeInput(session$ns("gene3dstructure"),
                       "Protein",
                       choices = gene_options[order(match(gene_options, data()$content))]
        )
      })
      
      ## TABLE
      output$title_structure3d_table <- renderText({paste0("PDB table for ", str_c(input$gene3dstructure, collapse = ", "))})
      output$structure3d_table <- DT::renderDataTable({
        shiny::validate(
          need(proteins %>% 
                 filter(gene_name %in% data()$content) %>% 
                 left_join(uniprot_pdb_table, by = c("uniprot_id" = "uniprot")) %>%
                 unnest(data) %>% 
                 nrow() > 0, 
               "No structure found for this protein/s")
        )
        
        DT::datatable(make_structure3d_table(input = data()) %>%
                        mutate(uniprot_id = uniprot_linkr2(uniprot_id),
                               pdb = pdb_linkr(pdb)) %>% 
                        dplyr::select('Query' = gene_name, 
                                      'UniProt ID' = uniprot_id, 
                                      'PDB ID' = pdb,
                                      'Name' = title, 
                                      'Organism' = organism) %>% 
                        dplyr::filter(Query == input$gene3dstructure),
                      escape = FALSE,
                      options = list(pageLength = 10),
                      selection = c("single"),
        )
      })
      
      ## PLOT
      rv <- reactiveValues(gene3dstructure = NULL,
                           pdb3dstructure = NULL,
                           color3dstructure = FALSE,
                           ribbon3dstructure = FALSE,
                           selection3dstructure = FALSE,
                           residue = 1:10,
                           chain = "A")
      
      #update value upon call
      observeEvent(input$update3d, {
        
        if(!is.null(input$structure3d_table_rows_selected)){
          pdb3dstructure <- make_structure3d_table(input = data()) %>% 
            dplyr::slice(input$structure3d_table_rows_selected) %>% 
            pull(pdb)
        } else {
          pdb3dstructure <- NULL
        }
        
        rv$gene3dstructure <- input$gene3dstructure 
        rv$pdb3dstructure <- pdb3dstructure
        rv$color3dstructure <- input$color3dstructure
        rv$ribbon3dstructure <- input$ribbon3dstructure
        rv$selection3dstructure <- input$selection3dstructure
        rv$residue <- input$residue
        rv$chain <- input$chain
      })
      
      plot3dprotein <- reactive({
        make_structure3d(input = data(),
                         gene_id = rv$gene3dstructure,
                         pdb_id = rv$pdb3dstructure,
                         color = rv$color3dstructure,
                         ribbon = rv$ribbon3dstructure,
                         selection = rv$selection3dstructure,
                         resi = rv$residue,
                         chain = rv$chain)
      })
      
      output$text_protein_structure3d <- renderText({paste0("Predicted 3D Structure for ", str_c(input$gene3dstructure, collapse = ", "))})
      output$protein_structure3D <- r3dmol::renderR3dmol({
        shiny::validate(
          need(proteins %>% 
                 filter(gene_name %in% data()$content) %>% 
                 left_join(uniprot_pdb_table, by = c("uniprot_id" = "uniprot")) %>%
                 unnest(data) %>% 
                 nrow() > 0, 
               "No structure found for this protein")
        )
        plot3dprotein()
        })
    }
  )
}

# AA RADIAL PLOT --------------------------------------------------------
radialPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_radial_plot")))),
    shinyWidgets::prettySwitch(inputId = ns("mean_relative"), 
                               "Show relative frequency to the mean", 
                               value = TRUE),
    fluidRow(plotOutput(outputId = ns("radial_plot"), height = "auto") %>% 
               withSpinnerColor(plot_type = "protein") #see shiny_helper.R
    )
  )
}

radialPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_radial_plot <- renderText({paste0("Amino Acid Signature for ", str_c(data()$content, collapse = ", "))})
      output$radial_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        make_radial(input = data(),
                    relative = input$mean_relative,
                    cluster = FALSE,
                    barplot = FALSE)
      },
      height = 550)
    }
  )
}

# AA BAR PLOT --------------------------------------------------------
AABarPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(ddh::make_legend("make_radial"), # this belongs to RADIAL PLOT
             actionLink(inputId = ns("aa_bar_click"), " View bar plot of amino acid signatures")),
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("aa_bar_click"), "'] != 0"), 
                     fluidRow(h4(textOutput(ns("text_aa_bar_plot")))),
                     shinyWidgets::prettySwitch(inputId = ns("bar_mean_relative"), 
                                                "Show relative frequency to the mean", 
                                                value = TRUE),
                     fluidRow(plotOutput(outputId = ns("aa_bar_plot"), height = "auto") %>% 
                                withSpinnerColor(plot_type = "protein") #see shiny_helper.R
                     ),
                     tags$br(),
                     fluidRow(ddh::make_legend("make_radial_bar"))
    )
  )
}

AABarPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_aa_bar_plot <- renderText({paste0("Amino Acid Signature for ", str_c(data()$content, collapse = ", "))})
      output$aa_bar_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        make_radial(input = data(),
                    relative = input$bar_mean_relative,
                    cluster = FALSE,
                    barplot = TRUE)
      },
      height = 550)
    }
  )
}

# CLUSTER AA RADIAL PLOT --------------------------------------------------------
ClusterradialPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("cluster_text_radial_plot")))),
    shinyWidgets::prettySwitch(inputId = ns("cluster_mean_relative"), 
                               "Show relative frequency to the mean", 
                               value = TRUE),
    fluidRow(plotOutput(outputId = ns("cluster_radial_plot"), height = "auto") %>% 
               withSpinnerColor(plot_type = "protein") #see shiny_helper.R
    ),
    tags$br()
  )
}

ClusterradialPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cluster_text_radial_plot <- renderText({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        
        title_text <- paste0("Amino Acid Signature for Cluster ", 
                             str_c(clust_num, collapse = ", "))
        
        return(title_text)
      })
      output$cluster_radial_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        make_radial(input = data(),
                    relative = input$cluster_mean_relative,
                    cluster = TRUE,
                    barplot = FALSE)
      },
      height = 550)
    }
  )
}

# CLUSTER AA BAR PLOT --------------------------------------------------------
ClusterAABarPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(ddh::make_legend("make_radial"), # this belongs to RADIAL PLOT
             actionLink(inputId = ns("aa_bar_cluster_click"), " View bar plot of cluster amino acid signatures")),
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("aa_bar_cluster_click"), "'] != 0"), 
                     fluidRow(h4(textOutput(ns("cluster_text_aa_bar_plot")))),
                     shinyWidgets::prettySwitch(inputId = ns("cluster_bar_mean_relative"), 
                                                "Show relative frequency to the mean", 
                                                value = TRUE),
                     fluidRow(plotOutput(outputId = ns("cluster_aa_bar_plot"), height = "auto") %>% 
                                withSpinnerColor(plot_type = "protein") #see shiny_helper.R
                     ),
                     tags$br(),
                     fluidRow(ddh::make_legend("make_radial_bar"))
    )
  )
}

ClusterAABarPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cluster_text_aa_bar_plot <- renderText({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        
        title_text <- paste0("Amino Acid Signature for Cluster ", 
                             str_c(clust_num, collapse = ", "))
        
        return(title_text)
      })
      output$cluster_aa_bar_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        make_radial(input = data(),
                    relative = input$cluster_bar_mean_relative,
                    cluster = TRUE,
                    barplot = TRUE)
      },
      height = 550)
    }
  )
}

# UMAP PLOT --------------------------------------------------------
UMAPPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(actionLink(inputId = ns("umap_click"), " View UMAP embeddings plot")),
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("umap_click"), "'] != 0"),
                     fluidRow(h4(textOutput(ns("text_umap_plot")))),
                     checkboxInput(inputId = ns("show_all_umap"), 
                                   label = "Show only selected clusters", value = FALSE),
                     checkboxInput(inputId = ns("labels_umap"), 
                                   label = "Show labels", value = FALSE),
                     fluidRow(plotOutput(outputId = ns("umap_plot"), height = "auto") %>% 
                                withSpinnerColor(plot_type = "protein") #see shiny_helper.R
                     ),
                     tags$br(),
                     fluidRow(ddh::make_legend("make_umap_plot"))
    )
  )
}

UMAPPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_umap_plot <- renderText({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        
        title_text <- paste0("UMAP Embeddings for Cluster ", 
                             str_c(clust_num, collapse = ", "))
        
        return(title_text)
      })
      output$umap_plot <- renderPlot({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        make_umap_plot(input = data(),
                       show_subset = input$show_all_umap,
                       labels = input$labels_umap)
      },
      height = 550)
    }
  )
}

# CLUSTER ENRICHMENT PLOT --------------------------------------------------------
ClusterEnrichmentPlot <- function(id) {
  ns <- NS(id)
  uiOutput(outputId = ns("conditional_clusterenrichmentplot"))
}

ClusterEnrichmentPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_clusterenrichmentplot <- renderUI({
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
            fluidRow(h4(textOutput(session$ns("text_cluster_enrich_plot")))),
            selectizeInput(session$ns("ontology_info_plot"),
                           "Ontology",
                           choices = c("Biological Process (GO)" = "BP",
                                       "Molecular Function (GO)" = "MF",
                                       "Cellular Component (GO)" = "CC"), 
                           selected = "BP"),
            sliderInput(session$ns("num_terms"), "Number of terms to show",
                        min = 10, max = 40, value = 20),
            fluidRow(plotOutput(outputId = session$ns("cluster_enrichment_plot"), height = "auto") %>% 
                       withSpinnerColor(plot_type = "protein") #see shiny_helper.R
            ),
            tags$br(),
            fluidRow(ddh::make_legend("make_cluster_enrich"))
          )
        } else {
          return(NULL)
        }
        
      })
      output$text_cluster_enrich_plot <- renderText({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"))
        
        if(length(clust_num) == 1) {
          title_text <- paste0("Enrichment Analysis Plot for Cluster ", 
                               clust_num)
        } else{
          title_text <- NULL
        }
        
        return(title_text)
      })
      output$cluster_enrichment_plot <- renderPlot({
        
        clust_num <- make_clustering_table(input = data()) %>%
          dplyr::pull(clust) %>%
          unique()
        
        shiny::validate(
          need(data()$content %in% proteins$gene_name, "Protein not found"),
          need(length(clust_num) == 1, 
               "More than one cluster identified in the table above. Cluster enrichment analysis is only available for individual gene queries or multiple gene queries belonging to the same cluster."))
        make_cluster_enrich(input = data(),
                            ontology = input$ontology_info_plot,
                            num_terms = input$num_terms)
      },
      height = 550)
    }
  )
}

# PROTEIN DOMAIN PLOT --------------------------------------------------------
proteinDomainPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_protein_domain_plot")))),
    fluidRow(column(uiOutput(outputId = ns("dom_choice")), width = 6), 
             column(uiOutput(outputId = ns("ptm_choice")), width = 6)
    ),
    fluidRow(plotOutput(outputId = ns("protein_domain_plot"), height = "auto") %>% 
               withSpinnerColor(plot_type = "protein") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_protein_domain")),
    tags$br()
  )
}

proteinDomainPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_protein_domain_plot <- renderText({paste0("Protein Domain Plot for ", str_c(data()$content, collapse = ", "))})
      output$dom_choice <- renderUI({
        
        protein_domain_data <- protein_domains %>%
          dplyr::filter(gene_name %in% data()$content)
        
        prots_dr <- protein_domain_data %>%
          filter(type %in% c("DOMAIN", "REGION")) %>%
          drop_na(description)
        
        selectizeInput(session$ns("dom_var"), "Protein Features (select):", 
                       choices = prots_dr %>% 
                         distinct(description) %>% 
                         pull(description),
                       multiple = TRUE,
                       selected = prots_dr %>%  
                         distinct(description) %>% 
                         dplyr::slice(1) %>% 
                         pull(description)
        ) 
      })
      output$ptm_choice <- renderUI({
        
        protein_domain_data <- protein_domains %>%
          dplyr::filter(gene_name %in% data()$content)
        
        prots_ptm <- protein_domain_data %>%
          filter(category == "PTM") %>% 
          drop_na(description)
        
        selectizeInput(session$ns("ptm_var"), "PTMs (select):", 
                       choices = prots_ptm %>%  
                         distinct(description) %>% 
                         pull(description),
                       multiple = TRUE,
                       selected = prots_ptm %>%  
                         distinct(description) %>% 
                         dplyr::slice(1) %>% 
                         pull(description)
        ) 
      })
      output$protein_domain_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% protein_domains$gene_name, "No data found for this protein."))
        make_protein_domain(input = data(),
                            dom_var = input$dom_var,
                            ptm_var = input$ptm_var)
      },
      height = function() length(data()$content) * 120 + 120
      )
    }
  )
}

# ANATOGRAM PLOTS -----
cellAnatogramPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_subcell_exp_plot")))),
    fluidRow(plotOutput(outputId = ns("cellanatogram")) %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
  )
  
}

cellAnatogramPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_subcell_exp_plot <- renderText({paste0("Sub-cellular expression of ", str_c(data()$content, collapse = ", "))})
      output$cellanatogram <- renderPlot({
        shiny::validate(
          need(data()$content %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, input = data())
      })
    }
  )
}

cellAnatogramFacetPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("cellanatogramfacet")) %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
  )
  
}

cellAnatogramFacetPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogramfacet <- renderPlot({
        shiny::validate(
          need(data()$content %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogramfacet(subcell, data()$content)
      })
    }
  )
}

maleAnatogramPlot <- function(id) {
  ns <- NS(id)
  uiOutput(outputId = ns("conditional_maleanatogramplot"))
}

maleAnatogramPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_maleanatogramplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "maleanatogram_plot"
        fieldnameMultiple <- "maleanatogram_plot_multiple"
        ImgFun <- "make_male_anatogram" # special case. original is make_colorful_anatogram
        query_type <- "gene"
        image_type <- "plot"
        SpinnerColor <- query_type
        PlotSizes <- paste0(plot_size_finder("make_male_anatogram"), "px")
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        }
      })
      output$maleanatogram_plot <- renderImage({
        shiny::validate(need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene."))
        list(src= load_image(input = data(), 
                             fun_name = "make_male_anatogram",  # special case. original is make_colorful_anatogram
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_male_anatogram")$plot_width,
             height = plot_size_finder("make_male_anatogram")$plot_height
        )
      }, deleteFile = FALSE)
      
      output$maleanatogram_plot_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene."))
        make_male_anatogram(input = data())
      })
    }
  )
}

femaleAnatogramPlot <- function(id) {
  ns <- NS(id)
  uiOutput(outputId = ns("conditional_femaleanatogramplot"))
}

femaleAnatogramPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_femaleanatogramplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "femaleanatogram_plot"
        fieldnameMultiple <- "femaleanatogram_plot_multiple"
        ImgFun <- "make_female_anatogram"
        query_type <- "gene"
        image_type <- "plot"
        SpinnerColor <- query_type
        PlotSizes <- paste0(plot_size_finder("make_female_anatogram"), "px") 
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        }
      })
      output$femaleanatogram_plot <- renderImage({
        shiny::validate(need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene."))
        list(src= load_image(input = data(), 
                             fun_name = "make_female_anatogram",  # special case. original is make_colorful_anatogram
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_female_anatogram")$plot_width,
             height = plot_size_finder("make_female_anatogram")$plot_height
        )
      }, deleteFile = FALSE)
      
      output$femaleanatogram_plot_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene."))
        make_female_anatogram(input = data())
      })
    }
  )
}

tissuePlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_tissue_exp_plot")))),
    fluidRow(plotOutput(outputId = ns("tissueplot"), height = "auto") %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
  )
  
}

tissuePlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_tissue_exp_plot <- renderText({paste0("Human tissue gene expression of ", str_c(data()$content, collapse = ", "))})
      output$tissueplot <- renderPlot({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene."))
        make_tissue(tissue, input = data())
      }, height = 1000)
    }
  )
}

# EXPRESSION PLOTS --------
cellGeneExpressionPlot <- function(id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cellgeneexpressionplot_title"))),
    uiOutput(ns("conditional_cellgeneexpressionplot")),
    tags$br()# ,
    # fluidRow(htmlOutput(outputId = ns("cellgeneexpressionplot_legend")))
  )
  
}

cellGeneExpressionPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellgeneexpressionplot_title <- renderText({paste0("Expression values for ", str_c(data()[[3]], collapse = ", "))})
      output$conditional_cellgeneexpressionplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cell_gene_plot"
        fieldnameMultiple <- "cell_gene_plot_multiple"
        ImgFun <- "make_cellexpression"
        query_type <- data()$type
        image_type <- "plot"
        SpinnerColor <- query_type
        PlotSizes <- paste0(plot_size_finder("make_cellexpression"), "px") 
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        }
      })
      output$cell_gene_plot <- renderImage({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(data()$content %in% gene) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(data()$content %in% cell_line) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
        }
        list(src= load_image(input = data(), 
                             fun_name = "make_cellexpression",  # special case. original is make_colorful_anatogram
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_cellexpression")$plot_width,
             height = plot_size_finder("make_cellexpression")$plot_height
        )
      }, deleteFile = FALSE)
      
      output$cell_gene_plot_multiple <- renderPlot({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
        }
        make_cellexpression(input = data(), var = "gene")
      })
    }
  )
}

cellProteinExpressionPlot <- function(id) {
  ns <- NS(id)
  tagList(
    h4(textOutput(outputId = ns("cellgeneproteinplot_title"))),
    uiOutput(ns("conditional_cellgeneproteinplot")),
    tags$br(),
    fluidRow(htmlOutput(outputId = ns("cellgeneproteinplot_legend")))
  )
  
}

cellProteinExpressionPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellgeneproteinplot_title <- renderText({paste0("Protein values for ", str_c(data()$content, collapse = ", "))})
      output$conditional_cellgeneproteinplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cell_gene_plot"
        fieldnameMultiple <- "cell_gene_plot_multiple"
        ImgFun <- "make_cellexpression"
        query_type <- data()$type
        image_type <- "plot"
        SpinnerColor <- "protein"
        PlotSizes <- paste0(plot_size_finder("make_cellexpression"), "px") 
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = SpinnerColor)
        }
      })
      output$cell_gene_plot <- renderImage({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
        }
        list(src= load_image(input = data(), 
                             fun_name = "make_cellexpression",  # special case. original is make_colorful_anatogram
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_cellexpression")$plot_width,
             height = plot_size_finder("make_cellexpression")$plot_height
        )
      }, deleteFile = FALSE)
      
      output$cell_gene_plot_multiple <- renderPlot({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
        }
        make_cellexpression(input = data(), var = "protein")
      })
      output$cellgeneproteinplot_legend <- renderText({paste0("<strong>Protein Values.</strong> Each point shows the ranked protein value across ", ifelse(data()$type == "gene", paste0(n_distinct(expression_long$X1)," cell lines."), paste0(n_distinct(expression_long$gene)," genes.")))})
    }
  )
}

cellGeneProteinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_gene_protein_plot")))),
    fluidRow(plotOutput(outputId = ns("cell_gene_protein_plot")) %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cellgeneprotein"))
  )
}

cellGeneProteinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_gene_protein_plot <- renderText({paste0("Gene vs. protein expression of ", str_c(data()$content, collapse = ", "))})
      output$cell_gene_protein_plot <- renderPlot({
        shiny::validate(
          need(expression_long %>% filter(gene %in% data()$content & !is.na(protein_expression)) %>%
                 nrow() > 0 && expression_long %>% filter(gene %in% data()$content & !is.na(gene_expression)) %>% nrow() > 0, "Cannot make the plot without all the data."))
        make_cellgeneprotein(input = data())
      })
    }
  )
}

cellLineGeneProteinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cellLine_gene_protein_plot")))),
    fluidRow(plotOutput(outputId = ns("cellLine_gene_protein_plot")) %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cellgeneprotein"))
  )
}

cellLineGeneProteinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cellLine_gene_protein_plot <- renderText({paste0("Gene vs. protein expression of ", str_c(data()$content, collapse = ", "))})
      output$cellLine_gene_protein_plot <- renderPlot({
        shiny::validate(
          need(expression_names$cell_line %in% data()$content, "No data found for this cell line."))
        # ggplotly(
        make_cellgeneprotein(input = data())#, tooltip = c("text")
        # )
      })
    }
  )
}

cellCoexpressionPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_coexp_plot")))),
    fluidRow(plotOutput(outputId = ns("cell_coexpression")) %>% 
               withSpinnerColor(plot_type = "cell") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cell_similarity"))
    )
}

cellCoexpressionPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$text_cell_coexp_plot <- renderText({paste0("Cell co-expression plot generated for ", str_c(data()$content, collapse = ", "))})
      output$cell_coexpression <- renderPlot({
        shiny::validate(
          need(data()$content %in% cell_line_exp_sim$cell1_name |
                 data()$content %in% cell_line_exp_sim$cell2_name, 
               "No data found for this cell line."))
        make_cell_similarity(input = data(), 
                             similarity = "expression")
      })      
    }
  )
}

# DEPENDENCY PLOTS -----
cellDependenciesPlot <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    # Scatterplot
    fluidRow(
      column(h4(textOutput(ns("text_cell_dep_plot"))), width = 12)), 
    fluidRow(column(textOutput(ns("essential_num")), width = 9),
             column(actionButton(ns("cell_dep_switch"), label = "Show dynamic plot"), width = 3)
    ),
    tags$br(),
    fluidRow(
      div(
        id = ns("cell_deps_static"),
        style = "padding-left:1%",
        uiOutput(outputId = ns("conditional_celldependenciesplot"))
      ),
      shinyjs::hidden(
        div(
          id = ns("cell_deps_dynamic"),
          style = "padding-left:1%",
          plotlyOutput(outputId = ns("cell_deps_dynamic_plot"))
        )
      )
    ),
    tags$br()
  )
}

cellDependenciesPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      observeEvent(input$cell_dep_switch, {
        
        if(input$cell_dep_switch %% 2 == 1){
          shinyjs::hide("cell_deps_static")
          shinyjs::show("cell_deps_dynamic")
          
          updateActionButton(session, "cell_dep_switch", label = "Show static plot")
          
        } else {
          shinyjs::hide("cell_deps_dynamic")
          shinyjs::show("cell_deps_static")
          
          updateActionButton(session, "cell_dep_switch", label = "Show dynamic plot")
        }
      })
      
      # Scatterplot
      output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", str_c(data()$content, collapse = ", "))})
      output$essential_num <- renderText({
        get_essential(input = data())
        #paste0("Essential in ", get_essential(input = data()), " cell lines")
      })
      output$conditional_celldependenciesplot <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cell_deps_static_plot"
        fieldnameMultiple <- "cell_deps_static_multiple_plot"
        ImgFun <- "make_celldeps"
        query_type <- data()$type
        image_type <- "plot"
        SpinnerColor <- query_type
        PlotSizes <- paste0(plot_size_finder(ImgFun), "px")
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = PlotSizes[2],
                      width = PlotSizes[1], inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotSizes[2],
                     width = PlotSizes[1]) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      #static
      output$cell_deps_static_plot <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        list(src= load_image(input = data(), 
                             fun_name = "make_celldeps", 
                             image_type = "plot"), #this calls the plot image, not the card
             width = plot_size_finder("make_celldeps")$plot_width,
             height = plot_size_finder("make_celldeps")$plot_height
        )
      }, deleteFile = FALSE)
      
      output$cell_deps_static_multiple_plot <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        withProgress(message = 'Almost there...', value = 1, {
          make_celldeps(input = data(), 
                        card = FALSE,
                        lineplot = TRUE,
                        scale = 0.3)
        })
      })
      #dynamic
      output$cell_deps_dynamic_plot <- renderPlotly({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        ggplotly(make_celldeps(input = data()), tooltip = "text" #,
                 # width = plot_size_finder("make_celldeps")$plot_width,
                 # height = plot_size_finder("make_celldeps")$plot_height
        )
      })
    }
  )
}

cellDependenciesDensityPlot <- function(id) {
  ns <- NS(id)
  tagList(
    h4(textOutput(outputId = ns("cell_dep_density_title"))),
    uiOutput(ns("density_plot_plot")),
    tags$br(),
    fluidRow(htmlOutput(outputId = ns("density_plot_legend")))
  )
}

cellDependenciesDensityPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_dep_density_title <- renderText({paste0("Dependency density plot for ", str_c(data()[[3]], collapse = ", "))})
      output$density_plot_plot <- renderUI({
        fluidRow(plotOutput(outputId = session$ns("cell_deps_density"),  height = "auto") %>% 
                   withSpinnerColor(plot_type = data()$type) #see shiny_helper.R
        )
        
      })
      output$cell_deps_density <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        make_cellbins(input = data())
      },
      height = function() length(data()$content) * 90 + 80)
      output$density_plot_legend <- renderText({paste0("<strong>Computed Densities.</strong> Kernel density estimate of dependency scores. Dependency scores across all ", ifelse(data()$type == "gene", "cell lines for queried genes", "genes for queried cell lines"), ", revealing overall influence of a gene on cellular fitness. The interval indicates the 95% quantile of the data, the dot indicates the median dependency score. The gray background highlights weak dependency values between -1 and 1.")})
    }
  )
}

cellDependenciesBarPlot <- function(id) {
  ns <- NS(id)
  tagList(
    h4(textOutput(outputId = ns("cell_dep_barplot_title"))),
    uiOutput(ns("bar_plot_plot")),
    tags$br(),
    fluidRow(htmlOutput(outputId = ns("bar_plot_legend")))
  )
}

cellDependenciesBarPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_dep_barplot_title <- renderText({paste0("Dependency barplot for ", str_c(data()$content, collapse = ", "))})
      output$bar_plot_plot <- renderUI({
        fluidRow(plotlyOutput(outputId = session$ns("cell_bar"),  height = "auto") %>% 
                   withSpinnerColor(plot_type = data()$type) #see shiny_helper.R
        )
      })
      output$cell_bar <- renderPlotly({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        ggplotly(make_cellbar(input = data()), tooltip = "text")
      })
      output$bar_plot_legend <- renderText({paste0("<strong>Dependency Bar Plot.</strong> Each bar shows the dependency scores of the queried ", ifelse(data()$type == "gene", "genes in a cell line.", "cell lines in a gene."), " Dependency scores less than -1 indicate a gene that is essential within a cell line. Dependency scores close to 0 mean no changes in fitness when the gene is knocked out. Dependency scores greater than 1 indicate gene knockouts lead to a gain in fitness.")})
    }
  )
}

cellDepsLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    # conditionalPanel(condition = paste0("output['", ns("click_deps_lin_high"), "'] == 1"),
    fluidRow(prettySwitch(inputId = ns("celllin_click_high"), "Show statistically significant")),
    # ),
    tags$br(),
    fluidRow(plotOutput(outputId = ns("cell_deps_lin"),  height = "auto") %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_lineage"))
  )
}

cellDepsLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      # output$click_deps_lin_high <- reactive({
      #   if(!is.null(data()$content)) {
      #     gene_len <- length(data()$content)
      #     return(gene_len)
      #   }
      # })
      # outputOptions(output, "click_deps_lin_high", suspendWhenHidden = FALSE)
      
      observeEvent(input$celllin_click_high, { #event to store the 'click'
      })
      output$cell_deps_lin <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No data found for this gene."))
        make_lineage(input = data(), 
                     highlight = input$celllin_click_high)
      },
      height = 550)
      observeEvent(input$sublin_click, { #event to store the 'click'
      })
    }
  )
}

cellDepsSubLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    # conditionalPanel(condition = paste0("output['", ns("click_deps_sublin_high"), "'] == 1"),
    fluidRow(prettySwitch(inputId = ns("cellsublin_click_high"), "Show statistically significant")),
    # ),
    tags$br(),
    fluidRow(plotOutput(outputId = ns("cell_deps_sublin"),  height = "auto") %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_sublineage"))
  )
}

cellDepsSubLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$cellsublin_click_high, { #event to store the 'click'
      })
      
      output$cell_deps_sublin <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No data found for this gene."))
        make_sublineage(input = data(),
                        highlight = input$cellsublin_click_high)
      },
      height = 1400)
    }
  )
}

cellDependenciesCorrPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_genecorr_plot")))),
    fluidRow(plotOutput(outputId = ns("gene_correlations")) %>% 
               withSpinnerColor(plot_type = "gene") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_correlation"))
  )
}

cellDependenciesCorrPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_genecorr_plot <- renderText({paste0("Gene correlation plot generated for ", str_c(data()$content, collapse = ", "))})
      output$gene_correlations <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_cor_nest$fav_gene, "No data found for this gene."))
        make_correlation(input = data())
      })      
    }
  )
}

compoundDependenciesPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compound_dep_plot")))),
    fluidRow(textOutput(ns("essential_num_compound"))),
    fluidRow(plotlyOutput(outputId = ns("cell_deps")) %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_celldeps"))
    )
}

compoundDependenciesPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_compound_dep_plot <- renderText({paste0("Viability plots generated for ", str_c(data()$content, collapse = ", "))})
      output$essential_num_compound <- renderText({
        paste0("Toxic in ", get_essential(input = data()), " cell lines")})
      output$cell_deps <- renderPlotly({
        shiny::validate(
          need(data()$content %in% prism_long$name, "No data found for this compound"))
        ggplotly(make_celldeps(input = data()), tooltip = "text")
      })      
    }
  )
}

compoundBinsPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("compound_bins"),  height = "auto") %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cellbins"))
  )
}

compoundBinsPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_bins <- renderPlot({
        shiny::validate(
          need(data()$content %in% prism_long$name, "")) #""left blank
        make_cellbins(input = data())
      },
      height = function() length(data()$content) * 90 + 80)
    }
  )
}

compoundLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("compound_lin"),  height = "auto") %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(tags$strong(ddh::make_legend("make_lineage"))
    )
  )
}

compoundLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_lin <- renderPlot({
        shiny::validate(
          need(data()$content %in% prism_long$name, "No data found for this gene."))
        make_lineage(input = data())
      },
      height = 550)
      observeEvent(input$compound_sublin_click, { #event to store the 'click'
      })
    }
  )
}

compoundSubLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(actionLink(inputId = ns("compound_sublin_click"), "View plot split into sublineages")),
    tags$br(),
    conditionalPanel(condition = paste0("input['", ns("compound_sublin_click"), "'] != 0"),
                     tags$br(),
                     fluidRow(plotOutput(outputId = ns("compound_sublin"),  height = "auto") %>% 
                                withSpinnerColor(plot_type = "compound") #see shiny_helper.R
                     ))
  )
}

compoundSubLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$compound_sublin_click, { #event to store the 'click'
      })
      output$compound_sublin <- renderPlot({
        shiny::validate(
          need(data()$content %in% prism_long$name, "No data found for this gene."))
        make_sublineage(input = data())
      },
      height = 1400)
    }
  )
}

compoundCorrelationPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compoundcorr_plot")))),
    fluidRow(plotOutput(outputId = ns("compound_correlations")) %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_correlation"))
    )
}

compoundCorrelationPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_compoundcorr_plot <- renderText({paste0("Compound correlation plot generated for ", str_c(data()$content, collapse = ", "))})
      output$compound_correlations <- renderPlot({
        shiny::validate(
          need(data()$content %in% prism_cor_nest$fav_drug, "No data found for this compound"))
        make_correlation(input = data())
      })      
    }
  )
}

cellCoessentialityPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_coess_plot")))),
    fluidRow(plotOutput(outputId = ns("cell_coessentiality")) %>% 
               withSpinnerColor(plot_type = "cell") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cell_similarity"))
    )
}

cellCoessentialityPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {

      output$text_cell_coess_plot <- renderText({paste0("Cell co-essentiality plot generated for ", str_c(data()$content, collapse = ", "))})
      output$cell_coessentiality <- renderPlot({
        shiny::validate(
          need(data()$content %in% cell_line_dep_sim$cell1_name |
                 data()$content %in% cell_line_dep_sim$cell2_name, 
               "No data found for this cell line."))
        make_cell_similarity(input = data())
      })      
    }
  )
}

# Exp v. Dep Plots -----
expdepPlot <- function(id) {
  ns <- NS(id)
  tagList(
    h4(textOutput(outputId = ns("exdep_plot_title"))),
    uiOutput(ns("exdep_plot_plot")),
    tags$br(),
    fluidRow(htmlOutput(outputId = ns("exdep_plot_legend")))
  )
}

expdepPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$exdep_plot_title <- renderText({paste0("Expression v. dependency plot for ", str_c(data()[[3]], collapse = ", "))})
      output$exdep_plot_plot <- renderUI({
        fluidRow(plotOutput(outputId = session$ns("cellexpdep"),  height = "500px") %>% 
                   withSpinnerColor(plot_type = data()$type) #see shiny_helper.R
        )
      })
      output$cellexpdep <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        make_expdep(input = data())
      })
      output$exdep_plot_legend <- renderText({paste0("<strong>Dependency versus Expression.</strong> Each point shows the dependency value compared to the expression value for ", ifelse(data()$type == "gene", "a cell line given a gene.", "a gene given a cell line."))})
    }
  )
}

# PUBMED PLOT --------------------------------------------------------
pubmedPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_pubmed_plot")))),
    fluidRow(textOutput(ns("text_pubmed_plot"))),
    uiOutput(ns("pubmed_plot_ui"))
  )
}

pubmedPlotServer <- function (id, data, session) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_pubmed_plot <- renderText({paste0("Publication history for ", str_c(data()$content, collapse = ", "))})
      output$text_pubmed_plot <- renderText({
        num <- make_pubmed_table(input = data())
        if(num[,1] != "No data available") {
          num <- nrow(num)
        } else {
          num <- 0
        }
        glue::glue('{num} annotated papers')
      })
      output$pubmed_plot_ui <- renderUI({
        fluidRow(plotOutput(outputId = session$ns("pubmed_plot"), height = 550) %>% 
                   withSpinnerColor(plot_type = data()$type) #see shiny_helper.R
        )
      })
      output$pubmed_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% pubmed$name, "No literature found"))
        make_pubmed(input = data())
      })      
    }
  )
}

pubmedCompoundPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_pubmed_compound_plot")))),
    fluidRow(textOutput(ns("text_pubmed_compound_plot"))),
    fluidRow(plotOutput(outputId = ns("pubmed_compound_plot")) %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    )
  )
}

pubmedCompoundPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_pubmed_compound_plot <- renderText({paste0("Publication history for ", str_c(data()$content, collapse = ", "))})
      output$text_pubmed_compound_plot <- renderText({
        num <- make_pubmed_table(input = data()) %>% count()
        glue::glue('{num} annotated papers')
      })
      output$pubmed_compound_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% pubmed$name, "no literature found"))
        make_pubmed(input = data())
      })      
    }
  )
}

pubmedCellLinePlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_pubmed_cell_line_plot")))),
    fluidRow(textOutput(ns("text_pubmed_cell_line_plot"))),
    fluidRow(plotOutput(outputId = ns("pubmed_cell_line_plot")) %>% 
               withSpinnerColor(plot_type = "cell") #see shiny_helper.R
    )
  )
}

pubmedCellLinePlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_pubmed_cell_line_plot <- renderText({paste0("Publication history for ", str_c(data()$content, collapse = ", "))})
      output$text_pubmed_cell_line_plot <- renderText({
        shiny::validate(
          need(data()$content %in% pubmed$name, ""))
        num <- make_pubmed_table(input = data()) %>% count()
        glue::glue('{num} annotated papers')
      })
      output$pubmed_cell_line_plot <- renderPlot({
        shiny::validate(
          need(data()$content %in% pubmed$name, "no literature found"))
        make_pubmed(input = data())
      })      
    }
  )
}
# CELL QUERY-----
## Cell Image Loader ----------
cellImage <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(imageOutput(outputId = ns("cell_image"), height = "auto") %>% 
               withSpinnerColor(plot_type = "cell") #see shiny_helper.R
    ),
    tags$br(),
    fluidRow(ddh::make_legend("make_cell_image"))
  )
}

cellImageServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_title <- renderText({paste0(str_c(data()$content[[1]], " Cell Line Images"))})
      output$text_cell_attribution <- renderText({
        if(file.exists(make_cell_image(input = data()))){
          atcc_id <- 
            cellosaurus %>% 
            dplyr::filter(name %in% data()$content[[1]]) %>% 
            dplyr::pull("ATCC")
          glue::glue('<a href="https://www.atcc.org/products/all/{atcc_id}.aspx#characteristics" target="_blank">Image from ATCC</a>')
        }
      })
        output$cell_image <- renderImage({
          shiny::validate(
            need(file.exists(make_cell_image(input = data())), "No cell image available"))
          list(src= make_cell_image(input = data()), 
               width = "100%")
        }, deleteFile = FALSE)
    }
      )
}

# FUNCTIONAL PLOT ----
cellFunctionalPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_functional_plot")))),
    fluidRow(prettySwitch(ns("remove_equal"), 
                          "Remove pathways represented for the same group of genes",
                          value = FALSE),
             numericInput(inputId = ns("num_genes"),
                          "Minumum number of genes required to represent a pathway",
                          value = 2),
             sliderInput(ns("num_terms"), 
                         "Number of pathways to show",
                         min = 2, max = 40, value = 10)
             ),
    fluidRow(plotlyOutput(outputId = ns("functional_plot"),  height = "auto") %>% 
      withSpinnerColor(plot_type = "cell")
      ), #see shiny_helper.R
    tags$br(),
    fluidRow(ddh::make_legend("make_functional_cell")),
    tags$br()
  )
}

cellFunctionalPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_functional_plot <- renderText({paste0("Differential pathway expression plot for ", str_c(data()$content, collapse = ", "))})
      output$functional_plot <- renderPlotly({
        ggplotly(
          make_functional_cell(input = data(),
                               num_pathways = input$num_terms,
                               num_genes = input$num_genes,
                               remove_equivalent_pathways = input$remove_equal),
          tooltip = "text"
        )
      })
    })
}

# CELL METADATA PLOT ----
cellMetadataPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("title_metadata_plot")))),
    fluidRow(selectizeInput(inputId = ns("metadata_factor"),
                            "Metadata factor:",
                            choices = c("Lineage" = "lineage",
                                        "Sublineage" = "sublineage")
                            ),
             sliderInput(ns("bonferroni_cutoff"), 
                         "Bonferroni cutoff",
                         min = 0.01, max = 1, value = 0.05)),
    fluidRow(plotOutput(outputId = ns("metadata_plot"),  height = "auto") %>% 
               withSpinnerColor(plot_type = "cell")
    ), #see shiny_helper.R
    tags$br(),
    fluidRow(ddh::make_legend("make_metadata_cell"))
  )
}

cellMetadataPlotServer <- function (id, data, type) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title_metadata_plot <- renderText({paste0("Lineage similarity plot for ", str_c(data()$content, collapse = ", "))})
      output$metadata_plot <- renderPlot({
        shiny::validate(
          need(make_cell_sim_table(similarity = type, 
                                   bonferroni_cutoff = input$bonferroni_cutoff,
                                   input = data()) %>% 
                 bind_rows() %>% 
                 nrow() > 0, 
               "No associations for this cell line.")
        )
        make_metadata_cell(input = data(),
                           cell_line_similarity = type,
                           metadata = input$metadata_factor,
                           bonferroni_cutoff = input$bonferroni_cutoff)
      },
      height = 550)
    })
}

# COMPOUND QUERY -----
## Compound structure plot --------------------------------------------------------
compoundStructure <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("compound_structure")) %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    )
  )
}

compoundStructureServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_structure <- renderPlot({
        shiny::validate(
          need(is.array(make_molecule_structure(input = data())), "No structure found for this compound."))
        make_molecule_structure(input = data())
      })      
    }
  )
}

