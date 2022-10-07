# SUMMARY CARDS -----
#shiny code to generate card modules to place on pages

card_contents_width = "200px"
card_contents_height = "300px"
card_title_height = "30px"
#card_text_height = "30px"
card_content_and_title_height = "330px"

card_cell_style <- paste0(
  "box-shadow: rgba(100, 100, 111, 0.2) 0px 7px 19px 0px; ", # card effect
  "text-align: center;", # horizontally center contents
  "padding: 4px;", # add spacing around content inside the card
  "margin-bottom: 4px;", # add spacing below the card
  "height: ", card_content_and_title_height
)

# Creates a flow layout where the children will have card styling
cardLayout <- function(...) {
  flowLayout(
    ...,
    cellArgs = list(style = card_cell_style)
  )
}

cardTitle <- function(...) {
  div(
    ...,
    style=paste0("height:", card_title_height, ";")
  )
}

divFlexAlignCenter <- function(title, ...) {
  div(
    cardTitle(title),
    ...,
    style="display: flex; align-items: center; flex-direction: column; height:100%;"
  )
}

# GENE -----

##gene barcode plot-----
barcodeDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    #really want to build a dynamic ahref here to point to the specific NFT
    #would also be cool to dynamically add "owned by 0x...", or "available to own" if not
    #probably will need OS API? 
    
    #barcode_url <- textOutput(ns("barcode_url"))

    imageOutput(outputId = ns("barcodedash"),
                height = card_contents_height,
                width = card_contents_width, inline = TRUE) %>%
      withSpinnerColor(plot_type = "gene") 
    )
}

barcodeDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      # #hope to make a dynamic a href to link to specific NFT, not just all NFTs
      # output$barcode_url <- renderText({
      #   glue::glue('https://www.datadrivenhypothesis.com/?show=gene&query_type=gene&symbol={data()$content}')
      # })
      output$barcodedash <- renderImage({
        shiny::validate(
          need(!is.null(load_image(input = data(), fun_name = "make_barcode", image_type = "card")), "No image found for this query.")
        )
        list(src = load_image(input = data(), fun_name = "make_barcode", image_type = "card"),
             width=card_contents_width,
             height = card_contents_width) #force to square
      }, deleteFile = FALSE)
    })
}

#load_image(input = list(type = "gene", content = c("ROCK1")), fun_name = "make_barcode", image_type = "plot")
#load_image(input = list(type = "gene", content = c("ROCK1", "ROCK2")), fun_name = "make_ideogram")


##ideogram plot-----
ideogramPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Gene",
    # not centering vertically this because image already has a top left margin
    uiOutput(outputId = ns("conditional_ideogramdash"))
  )
}

ideogramPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_ideogramdash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "ideogramdash"
        fieldnameMultiple <- "ideogramdash_multiple"
        ImgFun <- "make_ideogram"
        query_type <- "gene"
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      
      output$ideogramdash <- renderImage({
        shiny::validate(
          need(data()$content %in% gene_location$approved_symbol, "No data found for this gene.")
          #need(length(data()$content) == 1, "View the ideogram \nfor your \nmultigene query")
        )
        list(src = load_image(input = data(), fun_name = "make_ideogram"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$ideogramdash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% gene_location$approved_symbol, "No data found for this gene."))
        withProgress(message = 'Almost there...', value = 1, {
          make_ideogram(input = data(),
                        card = TRUE)
        })
      })
    })
}

#load_image(input = list(type = "gene", content = c("ROCK1")), fun_name = "make_ideogram")
#load_image(input = list(type = "gene", content = c("ROCK1", "ROCK2")), fun_name = "make_ideogram")

##Structure plots-----
structureDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Protein",
    uiOutput(outputId = ns("conditional_structuredash"))
  )
}

structureDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_structuredash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "structuredash"
        fieldnameMultiple <- "structuredash_multiple"
        ImgFun <- "make_structure"
        query_type <- "gene"
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$structuredash <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No structure found for this protein")
        )
        list(src= load_image(input = data(), fun_name = "make_structure"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$structuredash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No structure found for this protein"))
        withProgress(message = 'Almost there...', value = 1, {
          make_structure(input = data(), #raw function will return single card from multigene query
                         card = TRUE)
        })
      })
    })
}

## PUBMED plot-----
pubmedPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Literature",
    uiOutput(outputId = ns("conditional_pubmedgenedash"))
  )
}

pubmedPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_pubmedgenedash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "pubmed_gene_plot_dash"
        fieldnameMultiple <- "pubmed_gene_plot_dash_multiple"
        ImgFun <- "make_pubmed"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$pubmed_gene_plot_dash <- renderImage({
        shiny::validate(
          need(data()$content %in% pubmed$name, "No literature found"))
        list(src= load_image(input = data(), fun_name = "make_pubmed"), 
             width = card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$pubmed_gene_plot_dash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% pubmed$name, "No literature found"))
        withProgress(message = 'Almost there...', value = 1, {
          make_pubmed(input = data(), 
                      card = TRUE)
        })
      })
    }
  )
}

#EXPRESSION plots-----
##expression plot-----
cellExpressionPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Expression",
    uiOutput(outputId = ns("conditional_cellexpressionplotdash"))
  )
}

cellExpressionPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_cellexpressionplotdash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cellexpressiondash"
        fieldnameMultiple <- "cellexpressiondash_multiple"
        ImgFun <- "make_cellexpression"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$cellexpressiondash <- renderImage({
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
        list(src= load_image(input = data(), fun_name = "make_cellexpression"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$cellexpressiondash_multiple <- renderPlot({
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
        withProgress(message = 'Almost there...', value = 1, {
          make_cellexpression(input = data(), 
                              card = TRUE)
        })
      })
    })
}

##cell anatogram plot-----
cellAnatogramPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Subcellular Distribution",
    uiOutput(outputId = ns("conditional_cellanatogramplotdash"))
  )
}

cellAnatogramPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_cellanatogramplotdash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cellanatogramdash"
        fieldnameMultiple <- "cellanatogramdash_multiple"
        ImgFun <- "make_cellanatogram"
        query_type <- "gene"
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$cellanatogramdash <- renderImage({
        shiny::validate(
          need(data()$content %in% subcell$gene_name, "No subcellular location data for this gene.")
        )
        list(src= load_image(input = data(), fun_name = "make_cellanatogram"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$cellanatogramdash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% subcell$gene_name, "No subcellular location data for this gene."))
        withProgress(message = 'Almost there...', value = 1, {
          make_cellanatogram(input = data(), 
                             card = TRUE)
        })
      })
    })
}

##tissue anatogram plot-----
tissueAnatogramPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Tissue Distribution",
    uiOutput(outputId = ns("conditional_tissueanatogramplotdash"))
  )
}

tissueAnatogramPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_tissueanatogramplotdash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "tissueanatogramdash"
        fieldnameMultiple <- "tissueanatogramdash_multiple"
        ImgFun <- "make_female_anatogram"
        query_type <- "gene"
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$tissueanatogramdash <- renderImage({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No data found for this gene.")
          #need(length(data()$content) == 1, "View the anatogram \nfor your \nmultigene query")
        )
        list(src= load_image(input = data(), fun_name = "make_female_anatogram"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$tissueanatogramdash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No data found for this gene."))
        withProgress(message = 'Almost there...', value = 1, {
          make_female_anatogram(input = data(), 
                                  card = TRUE)
        })
      })
    })
}

#load_image(input = list(type = "gene", content = c("ROCK1")), fun_name = "make_female_anatogram")
#load_image(input = list(type = "gene", content = c("ROCK1", "ROCK2")), fun_name = "make_female_anatogram")

#DEPENDENCY plots-----
##cell dependencies plot-----
cellDependenciesPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Cell Dependencies",
    uiOutput(outputId = ns("conditional_celldependenciesplotdash"))
  )
}

cellDependenciesPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_celldependenciesplotdash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "cell_depdash"
        fieldnameMultiple <- "cell_depdash_multiple"
        ImgFun <- "make_celldeps"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$cell_depdash <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        list(src = load_image(input = data(), fun_name = "make_celldeps"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$cell_depdash_multiple <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        withProgress(message = 'Almost there...', value = 1, {
          make_celldeps(input = data(), 
                        card = TRUE,
                        scale = 0.3)
        })
      })
    }
  )
}

##cell dependencies table-----
cellDependenciesTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Co-essential genes",
    gt_output(outputId = ns("deptabledash"))
  )
}

cellDependenciesTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deptabledash <- render_gt({
        shiny::validate(
          need(nrow(make_top_table(input = data())) != 0, "No data found for this gene."))
        gt::gt(make_top_table(input = data()) %>% 
                 dplyr::mutate("Rank" = row_number()) %>% 
                 dplyr::select("Rank", "Gene" = "gene", "R^2" = "r2") %>% 
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}
##cell dependencies graph-----
cellDependenciesGraphDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependencies Graph",
    visNetworkOutput(outputId = ns("depgraphdash"))
  )
}

cellDependenciesGraphDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$depgraphdash <- renderVisNetwork({
        if(nrow(make_top_table(input = data())) == 0 &
           !(data()$content %in% unique(cell_line_dep_sim$cell1_name))) {
          make_empty_graph()
        } else {
          make_graph(input = data(), 
                     threshold = 10, deg = 2, corrType = "Positive", 
                     card = TRUE)
        }
      })
    })
}

##cell coexpression graph-----
cellExpressionGraphDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Co-Expression Graph",
    visNetworkOutput(outputId = ns("expgraphdash"))
  )
}

cellExpressionGraphDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$expgraphdash <- renderVisNetwork({
        if(nrow(make_top_table(input = data())) == 0 &
           !(data()$content %in% unique(cell_line_exp_sim$cell1_name))) {
          make_empty_graph()
        } else {
          make_graph(input = data(), 
                     threshold = 10, 
                     deg = 2, 
                     corrType = "Positive", 
                     cell_line_similarity = "expression",
                     card = TRUE)
        }
      })
    })
}

##cell functional plot-----
cellFunctionalPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Pathway Expression",
    plotOutput(outputId = ns("functional_plot"))
  )
}

cellFunctionalPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
    output$functional_plot <- renderPlot({
      make_functional_cell(input = data(),
                           card = TRUE,
                           num_pathways = 5,
                           remove_equivalent_pathways = TRUE)
      })
    })
}

#DRUG table-----
geneDrugsCorTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Correlated Drugs",
    gt_output(outputId = ns("genedrugscortabledash"))
  )
}

geneDrugsCorTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$genedrugscortabledash <- render_gt({
        shiny::validate(
          need(data()$content %in% gene_drugs_cor_table$fav_gene, "No gene data found for this gene."))
        gt::gt(make_gene_drugs_cor_table(input = data()) %>%
                 dplyr::mutate("Rank" = row_number()) %>%
                 dplyr::select(Rank, Drug = drug) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width)
    }
  )
}

cellDrugsTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Associated Drugs",
    gt_output(outputId = ns("celldrugstabledash"))
  )
}

cellDrugsTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$celldrugstabledash <- render_gt({
        shiny::validate(
          need(prism_long %>% 
                 dplyr::rename(X1 = 1) %>% 
                 left_join(expression_meta, by = "X1") %>% 
                 pull(cell_line) %in% data()$content,
               "No drug data found for this cell line."))
        gt::gt(make_cell_drugs_table(input = data()) %>% 
                 dplyr::mutate("Rank" = row_number()) %>%
                 dplyr::select(Rank, Drug = name) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width)
    }
  )
}

cellMetabolitesTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Metabolite Levels",
    gt_output(outputId = ns("cellmetabolitestabledash"))
  )
}

cellMetabolitesTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellmetabolitestabledash <- render_gt({
        shiny::validate(
          need(cell_metabolites %>% 
                 left_join(expression_names, by = c("DepMap_ID" = "X1")) %>% 
                 filter(cell_line %in% data()$content) %>% 
                 nrow() > 0, 
               "No metabolites found that associate with this cell line.")
          )
        gt::gt(make_metabolite_table(input = data()) %>% 
                 dplyr::mutate("Rank" = row_number()) %>%
                 dplyr::select(Rank, Metabolite = metabolite) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width)
    }
  )
}

# CELL -----
##cell image card-----
cellImageDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Cell",
    imageOutput(outputId = ns("cellimagedash"),
                height = card_contents_height,
                width = card_contents_width, inline = TRUE) %>%
      withSpinnerColor(plot_type = "cell") 
  )
}

cellImageDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellimagedash <- renderImage({
        shiny::validate(
          need(file.exists(make_cell_image(input = data())), "No image found for this cell")
        )
        list(src = load_image(input = data(), fun_name = "make_cell_image", image_type = "card"),
             width=card_contents_width,
             height = card_contents_width) #make square #card_contents_height
      }, deleteFile = FALSE)
    })
}

#cell dependencies table-----
compoundDependenciesTableDash <- function(id) {
  ns <- NS(id)
  gt_output(outputId = ns("depcompoundtabledash"))
}

compoundDependenciesTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$depcompoundtabledash <- render_gt({
        shiny::validate(
          need(nrow(make_compound_table(input = data())) != 0, "No data found for this drug"))
        gt::gt(make_compound_table(input = data()) %>% 
                 dplyr::mutate("Rank" = row_number()) %>% 
                 dplyr::select("Rank", "Drug" = "name", "R^2" = "r2") %>% 
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width)
    }
  )
}

# COMPOUND -----
#structure plot-----
compoundStructureDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Compound Structure",
    uiOutput(outputId = ns("conditional_compoundstructuredash"))
  )
}

compoundStructureDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_compoundstructuredash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "compounddash"
        fieldnameMultiple <- "compounddash_multiple"
        ImgFun <- "make_molecule_structure"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$compounddash <- renderImage({
        #shiny::validate(
        #  need(is.array(make_molecule_structure(compound = data()$content)), "No structure found for this compound."))
        # withProgress(message = 'Shiny molecule comin up...', value = 1, {
        # make_molecule_structure(compound = data()$content)
        list(src = load_image(input = data(), fun_name = "make_molecule_structure"), 
             width = card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$compounddash_multiple <- renderPlot({
        withProgress(message = 'Almost there...', value = 1, {
          make_molecule_structure(input = data(), 
                                  card = TRUE)
        })
      })
      
      # })
    })
}

#compound dependencies plot-----
compoundDependenciesPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Compound Dependencies",
    uiOutput(outputId = ns("conditional_compounddependenciesplotedash"))
  )
}

compoundDependenciesPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_compounddependenciesplotedash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "compound_depdash"
        fieldnameMultiple <- "compound_depdash_multiple"
        ImgFun <- "make_celldeps"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$compound_depdash <- renderImage({
        shiny::validate(need(data()$content %in% prism_long$name, "No data found for this gene."))
        list(src= load_image(input = data(), fun_name = "make_celldeps"), 
             width = card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$compound_depdash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% prism_long$name, "No data found for this gene."))
        withProgress(message = 'Almost there...', value = 1, {
          make_celldeps(input = data(), 
                        card = TRUE)
        })
      })
    }
  )
}

#compound dependencies graph-----
compoundDependenciesGraphDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependencies Graph",
    visNetworkOutput(outputId = ns("compoundgraphdash")) %>% 
      withSpinnerColor(plot_type = "compound") #see shiny_helper.R
  )
}

compoundDependenciesGraphDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compoundgraphdash <- renderVisNetwork({
        shiny::validate(
          need(nrow(make_compound_table(input = data())) != 0, "No data found."))
        make_graph(input = data(), threshold = 10, deg = 2, corrType = "Positive")
      })
    })
}

#pubmed compound plot-----
pubmedCompoundPlotDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Pubmed Compound",
    uiOutput(outputId = ns("conditional_pubmedcompoundplotedash"))
  )
}

pubmedCompoundPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_pubmedcompoundplotedash <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "pubmed_compound_plot_dash"
        fieldnameMultiple <- "pubmed_compound_plot_dash_multiple"
        ImgFun <- "make_pubmed"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$pubmed_compound_plot_dash <- renderImage({
        shiny::validate(
          need(data()$content %in% pubmed$name, "no literature found"))
        list(src= load_image(input = data(), fun_name = "make_pubmed"), 
             width = card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$pubmed_compound_plot_dash_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% pubmed$name, "no literature found"))
        withProgress(message = 'Almost there...', value = 1, {
          make_pubmed(input = data(), 
                      card = TRUE)
        })
      })
    }
  )
}

#DRUGS-----
drugGenesCorTableDash <- function(id) {
  ns <- NS(id)
  gt_output(outputId = ns("druggenescortabledash"))
}

drugGenesCorTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$druggenescortabledash <- render_gt({
        shiny::validate(
          need(data()$content %in% drug_genes_cor_table$fav_drug, "No gene data found for this compound."))
        gt::gt(make_drug_genes_cor_table(drug = data()$content) %>% 
                 dplyr::mutate("Rank" = row_number()) %>% 
                 dplyr::select(Rank, Gene = gene, 'R^2'=r2) %>% 
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

# TAB CARDS -----
## GENE -----
###DEPENDENCY tables-----
geneGoTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Go Table",
    gt_output(outputId = ns("genegotabletab"))
  )
}
geneGoTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$genegotabletab <- render_gt({
        shiny::validate(
          need(data()$content %in% unique(unlist(pathways$data, use.names = FALSE)), "GO Term pathways")
        )
        gt::gt(make_pathway_list(table_name = pathways, input = data()) %>% 
                 #dplyr::mutate(go = map_chr(go, internal_link))  %>% #from fun_helper.R
                 dplyr::select(Pathway = pathway, GO = go) %>%
                 dplyr::slice(1))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}


###PROTEIN plots-----
sizePlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Size",
    uiOutput(outputId = ns("conditional_sizeplottab"))
  )
}

sizePlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_sizeplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "sizetab"
        fieldnameMultiple <- "sizetab_multiple"
        ImgFun <- "make_proteinsize"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- "protein"
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$sizetab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No size found for this protein")
        )
        list(src= load_image(input = data(), fun_name = "make_proteinsize"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$sizetab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No size found for this protein"))
        withProgress(message = 'Almost there...', value = 1, {
          make_proteinsize(input = data(),
                           card = TRUE)
        })
      })
    })
}

sequencePlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Sequence",
    uiOutput(outputId = ns("conditional_sequenceplottab"))
  )
}

sequencePlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_sequenceplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "sequencetab"
        fieldnameMultiple <- "sequencetab_multiple"
        ImgFun <- "make_sequence"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- "protein"
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$sequencetab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No sequence found for this protein.")
        )
        list(src= load_image(input = data(), fun_name = "make_sequence"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$sequencetab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No sequence found for this protein."))
        withProgress(message = 'Almost there...', value = 1, {
          make_sequence(input = data(),
                        card = TRUE)
        })
      })
    })
}

signaturePlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Signature",
    uiOutput(outputId = ns("conditional_signatureplottab"))
  )
}

signaturePlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_signatureplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "signaturetab"
        fieldnameMultiple <- "signaturetab_multiple"
        ImgFun <- "make_radial"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- "protein"
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$signaturetab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No signature found for this gene.")
        )
        list(src= load_image(input = data(), fun_name = "make_radial"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$signaturetab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No signature found for this gene."))
        withProgress(message = 'Almost there...', value = 1, {
          make_radial(input = data(),
                      card = TRUE)
        })
      })
    })
}

structurePlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Structure",
    uiOutput(outputId = ns("conditional_structureplottab"))
  )
}

structurePlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_structureplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "structureplottab"
        fieldnameMultiple <- "structureplottab_multiple"
        ImgFun <- "make_structure"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- "protein"
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$structureplottab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No structure found for this protein")
        )
        list(src= load_image(input = data(), fun_name = "make_structure"),
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$structureplottab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No structure found for this protein"))
        withProgress(message = 'Almost there...', value = 1, {
          make_structure(input = data(),
                         card = TRUE)
        })
      })
    })
}

###EXPRESSION tabs-----
cellGeneExpressionTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Gene Expression",
    gt_output(outputId = ns("geneexptabletab"))
  )
}

cellGeneExpressionTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$geneexptabletab <- render_gt({
        if(data()$type == "gene") {
          shiny::validate(
            need(nrow(make_expression_table(input = data(), var = "gene")) > 0, 
                 "No expression data found for this gene."))
          gt::gt(make_expression_table(input = data(), var = "gene") %>% 
                   dplyr::slice(1:5) %>% 
                   dplyr::select(-Lineage, -Subtype))
        } else if (data()$type == "cell") {
          shiny::validate(
            need(nrow(make_expression_table(input = data(), var = "gene")) > 0, 
                 "No expression data found for this cell line."))
          gt::gt(make_expression_table(input = data(), var = "gene") %>% 
                   dplyr::slice(1:5) %>% 
                   dplyr::select(-`Gene Name`))
        }
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellProteinExpressionPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Protein Expression",
    uiOutput(outputId = ns("conditional_cellproteinexpressionplottab"))
  )
}

cellProteinExpressionPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_cellproteinexpressionplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "proteinexpressionplottab"
        fieldnameMultiple <- "proteinexpressionplottab_multiple"
        ImgFun <- "make_cellexpression"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- "protein" #manual override to protein
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$proteinexpressionplottab <- renderImage({
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

        list(src= load_image(input = data(), fun_name = "make_cellexpression"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$proteinexpressionplottab_multiple <- renderPlot({
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
        withProgress(message = 'Almost there...', value = 1, {
          make_cellexpression(input = data(),
                              var = "protein",
                              card = TRUE)
        })
      })
    })
}

cellProteinExpressionTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Protein Expression",
    gt_output(outputId = ns("proteinexptabletab"))
  )
}

cellProteinExpressionTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$proteinexptabletab <- render_gt({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this gene.")
          )
          gt::gt(make_expression_table(input = data(), var = "protein") %>% 
                   dplyr::slice(1:5) %>% 
                   dplyr::select(-Lineage, -Subtype))
          
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No protein data found for this cell line.")
          )
          gt::gt(make_expression_table(input = data(), var = "protein") %>% 
                   dplyr::slice(1:5) %>% 
                   dplyr::select(-`Gene Name`))
        }
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellGeneProteinPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Gene v. Protein",
    uiOutput(outputId = ns("conditional_cellgeneproteinplottab"))
  )
}
cellGeneProteinPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_cellgeneproteinplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "geneproteinplottab"
        fieldnameMultiple <- "geneproteinplottab_multiple"
        ImgFun <- "make_cellgeneprotein"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$geneproteinplottab <- renderImage({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No data found for this cell line.")
          )
        }
        list(src= load_image(input = data(), fun_name = "make_cellgeneprotein"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$geneproteinplottab_multiple <- renderPlot({
        if(data()$type == "gene") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   filter(gene %in% data()$content) %>% 
                   nrow() > 0, "No data found for this gene.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(expression_long %>% 
                   drop_na(protein_expression) %>% 
                   left_join(expression_names, by = "X1") %>% 
                   filter(cell_line %in% data()$content) %>% 
                   nrow() > 0, "No data found for this cell line.")
          )
        }
        withProgress(message = 'Almost there...', value = 1, {
          make_cellgeneprotein(input = data(),
                               card = TRUE)
        })
      })
    })
}

tissuePlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Tissue Plot",
    uiOutput(outputId = ns("conditional_tissueplottab"))
  )
}
tissuePlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_tissueplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "tissueplottab"
        fieldnameMultiple <- "tissueplottab_multiple"
        ImgFun <- "make_tissue"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$tissueplottab <- renderImage({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene.")
        )
        list(src= load_image(input = data(), fun_name = "make_tissue"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$tissueplottab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene.")
        )
        withProgress(message = 'Almost there...', value = 1, {
          make_tissue(input = data(),
                      card = TRUE)
        })
      })
    })
}

tissueTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Expression Values",
    gt_output(outputId = ns("tissuetabletab"))
  )
}
tissueTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tissuetabletab <- render_gt({
        shiny::validate(
          need(data()$content %in% tissue$gene_name, "No Tissue-specific expression data for this gene.")
        )
        gt::gt(make_humananatogram_table(input = data()) %>% 
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineExpressionPosTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Positive Associations",
    gt_output(outputId = ns("cell_exppostabletab"))
  )
}

cellLineExpressionPosTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_exppostabletab <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data(),
                                        similarity = "expression")$top_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data(),
                              similarity = "expression")$top_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineExpressionNegTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Negative Associations",
    gt_output(outputId = ns("cell_expnegtabletab"))
  )
}
cellLineExpressionNegTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_expnegtabletab <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data(),
                                        similarity = "expression")$bottom_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data(),
                              similarity = "expression")$bottom_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

###DEPENDENCY plots-----
cellDependenciesBarPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependencies Barplot",
    uiOutput(outputId = ns("conditional_celldependenciesbarplottab"))
  )
}

cellDependenciesBarPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_celldependenciesbarplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "depbardash"
        fieldnameMultiple <- "depbardash_multiple"
        ImgFun <- "make_cellbar"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            gene_symbols = data()$gene_symbols),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$gene_symbols) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      
      output$depbardash <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        list(src= load_image(input = data(), fun_name = "make_cellbar"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$depbardash_multiple <- renderPlot({
        shiny::validate(
          need(length(data()[[3]]) > 0, "No info found for this query.")
        )
        withProgress(message = 'Almost there...', value = 1, {
          make_cellbar(input = data(),
                       card = TRUE,
                       scale = 0.2)
        })
      })
    }
  )
}

cellDependenciesDensityPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependencies Density",
    uiOutput(outputId = ns("conditional_celldependenciesdensityplottab"))
  )
}

cellDependenciesDensityPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_celldependenciesdensityplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "densityplottab"
        fieldnameMultiple <- "densityplottab_multiple"
        ImgFun <- "make_cellbins"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$densityplottab <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        list(src= load_image(input = data(), fun_name = "make_cellbins"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      
      output$densityplottab_multiple <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        withProgress(message = 'Almost there...', value = 1, {
          make_cellbins(input = data(),
                        card = TRUE)
        })
      })
    })
}

cellDepsLinPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Lineage",
    uiOutput(outputId = ns("conditional_celldepslinplottab"))
  )
}
cellDepsLinPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_celldepslinplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "lineageplottab"
        fieldnameMultiple <- "lineageplottab_multiple"
        ImgFun <- "make_lineage"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$lineageplottab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No plot found for this protein")
        )
        list(src= load_image(input = data(), fun_name = "make_lineage"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$lineageplottab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No plot found for this protein"))
        withProgress(message = 'Almost there...', value = 1, {
          make_lineage(input = data(),
                       card = TRUE)
        })
      })
    })
}

cellDepsSubLinPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Sublineage",
    uiOutput(outputId = ns("conditional_celldepssublinplottab"))
  )
}
cellDepsSubLinPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_celldepssublinplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "sublineageplottab"
        fieldnameMultiple <- "sublineageplottab_multiple"
        ImgFun <- "make_sublineage"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$sublineageplottab <- renderImage({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No plot found for this protein")
        )
        list(src= load_image(input = data(), fun_name = "make_sublineage"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$sublineageplottab_multiple <- renderPlot({
        shiny::validate(
          need(data()$content %in% achilles_long$gene, "No plot found for this protein"))
        withProgress(message = 'Almost there...', value = 1, {
          make_sublineage(input = data(),
                          card = TRUE)
        })
      })
    })
}

expdepPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Expression v. Dependency",
    uiOutput(outputId = ns("conditional_expdepplottab"))
  )
}
expdepPlotTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$conditional_expdepplottab <- renderUI({
        #### CARD PARAMETERS
        fieldnameSingle <- "expdepplottab"
        fieldnameMultiple <- "expdepplottab_multiple"
        ImgFun <- "make_expdep"
        query_type <- data()$type
        image_type <- "card"
        SpinnerColor <- query_type
        ImgWidth <- card_contents_width
        PlotWidth <- card_contents_width
        ImgHeight <- card_contents_height
        PlotHeight <- card_contents_height
        ####
        img_path <- load_image(input = list(type = query_type,
                                            content = data()$content),
                               fun_name = ImgFun,
                               image_type = image_type)
        
        if((length(data()$content) == 1) & (!is.null(img_path))) {
          imageOutput(outputId = session$ns(fieldnameSingle),
                      height = ImgHeight,
                      width = ImgWidth, inline = TRUE) %>%
            withSpinnerColor(plot_type = query_type)
        } else {
          plotOutput(outputId = session$ns(fieldnameMultiple),
                     height = PlotHeight,
                     width = PlotWidth) %>%
            withSpinnerColor(plot_type = query_type)
        }
      })
      output$expdepplottab <- renderImage({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        list(src= load_image(input = data(), fun_name = "make_expdep"), 
             width=card_contents_width,
             height = card_contents_height) #defined above
      }, deleteFile = FALSE)
      output$expdepplottab_multiple <- renderPlot({
        shiny::validate(need(data()$content %in% achilles_long$gene |
                               data()$content %in% expression_meta$cell_line, 
                             "No data found for this query."))
        withProgress(message = 'Almost there...', value = 1, {
          make_expdep(input = data(),
                      card = TRUE)
        })
      })
    })
}

###DEPENDENCY tables-----
cellDependenciesTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependency Table",
    gt_output(outputId = ns("deptabletab"))
  )
}
cellDependenciesTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deptabletab <- render_gt({
        shiny::validate(
          need(nrow(make_dep_table(input = data())) > 0, "No data found for this gene."))
        gt::gt(make_dep_table(input = data()) %>% 
                 dplyr::select(`Cell Line`, contains(data()$content)) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellDependenciesPosTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Positive Correlations",
    gt_output(outputId = ns("deppostabletab"))
  )
}
cellDependenciesPosTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deppostabletab <- render_gt({
        shiny::validate(
          need(data()$content %in% master_top_table$fav_gene, "No data found for this gene."))
        gt::gt(make_top_table(input = data()) %>% 
                 # dplyr::mutate("Rank" = row_number()) %>% 
                 dplyr::select(c("Gene" = "gene", "Z-Score" = "z_score")) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellDependenciesPosPathwayTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Positive Enrichment",
    gt_output(outputId = ns("deppospathwaystab"))
  )
}
cellDependenciesPosPathwayTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deppospathwaystab <- render_gt({
        shiny::validate(
          need(data()$content %in% master_positive$fav_gene, "No data found for this gene."))
        gt::gt(make_enrichment_top(input = data()) %>% 
                 dplyr::select(`Gene List`) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellDependenciesNegTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Negative Correlations",
    gt_output(outputId = ns("depnegtabletab"))
  )
}
cellDependenciesNegTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$depnegtabletab <- render_gt({
        shiny::validate(
          need(data()$content %in% master_bottom_table$fav_gene, "No data found for this gene."))
        gt::gt(make_bottom_table(input = data()) %>% 
                 # dplyr::mutate("Rank" = row_number()) %>% 
                 dplyr::select(c("Gene" = "gene", "Z-Score" = "z_score")) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellDependenciesNegPathwayTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Negative Correlations",
    gt_output(outputId = ns("depnegpathwaystab"))
  )
}
cellDependenciesNegPathwayTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$depnegpathwaystab <- render_gt({
        shiny::validate(
          need(data()$content %in% master_negative$fav_gene, "No data found for this gene."))
        gt::gt(make_enrichment_bottom(input = data()) %>% 
                 dplyr::select(`Gene List`) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}


## CELL -----
# dash
cellSummaryTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Lineage Table",
    gt_output(outputId = ns("cell_summarydash"))
  )
}

cellSummaryTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summarydash <- render_gt({
        shiny::validate(
          need(data()$content %in% expression_meta$cell_line, 
          "No data found for this cell line.")
        )
        gt::gt(
          make_cell_line_table(input = data()) %>%
            dplyr::select(1, 3) %>%
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineDependenciesPosTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Co-essential Cell Lines",
    gt_output(outputId = ns("cell_deptabledash"))
  )
}

cellLineDependenciesPosTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_deptabledash <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data())$top_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data())$top_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineExpressionPosTableDash <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Co-expressed Cell Lines",
    gt_output(outputId = ns("cell_exppostabletab"))
  )
}

cellLineExpressionPosTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_exppostabletab <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data(),
                                        similarity = "expression")$top_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data(),
                              similarity = "expression")$top_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

# tabs
cellLineDependenciesPosTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Positive Associations",
    gt_output(outputId = ns("cell_deppostabletab"))
  )
}

cellLineDependenciesPosTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_deppostabletab <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data())$top_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data())$top_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineDependenciesNegTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Negative Associations",
    gt_output(outputId = ns("cell_depnegtabletab"))
  )
}
cellLineDependenciesNegTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_depnegtabletab <- render_gt({
        shiny::validate(
          need(nrow(make_cell_sim_table(input = data())$bottom_table) > 0, 
               "No data found for this cell line."))
        gt::gt(
          make_cell_sim_table(input = data())$bottom_table %>%
            dplyr::rename("Cell" = "cell2_name", "Bonferroni" = "bonferroni") %>% 
            dplyr::select(Cell, Bonferroni) %>%
            dplyr::arrange(Bonferroni) %>% 
            dplyr::mutate_if(is.numeric, ~ signif(., digits = 3)) %>% 
            dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineDependenciesTableTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Dependency Table",
    gt_output(outputId = ns("deptabletab"))
  )
}
cellLineDependenciesTableTabServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deptabletab <- render_gt({
        shiny::validate(
          need(nrow(make_dep_table(input = data())) > 0, "No data found for this cell line."))
        gt::gt(make_dep_table(input = data()) %>% 
                 dplyr::rename(Gene = gene) %>% 
                 dplyr::select(Gene, contains(data()$content)) %>%
                 dplyr::slice(1:5))
      },
      height = card_contents_height,
      width = card_contents_width
      )
    }
  )
}

cellLineMetadataPlotTab <- function(id) {
  ns <- NS(id)
  divFlexAlignCenter(
    "Lineage Similarity",
    plotOutput(outputId = ns("metadata_plot_tab"))
  )
}

cellLineMetadataPlotTabServer <- function (id, data, type) {
  moduleServer(
    id,
    function(input, output, session) {
      output$metadata_plot_tab <- renderPlot({
        shiny::validate(
          need(make_cell_sim_table(similarity = type, 
                                   bonferroni_cutoff = 0.05,
                                   input = data()) %>% 
                 bind_rows() %>% 
                 nrow() > 0, 
               "No associations for this cell line.")
          )
        withProgress(message = 'Almost there...', value = 1, {
          make_metadata_cell(input = data(),
                             cell_line_similarity = type,
                             card = TRUE
                             )
        })
      })
    })
}

## COMPOUND -----
