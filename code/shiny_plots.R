# Plots -----
cellDependenciesPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_dep_plot")))),
    fluidRow(plotlyOutput(outputId = ns("cell_deps"))),
    tags$br(),
    fluidRow(tags$strong(plot_celldeps_title), plot_celldeps_legend))
}

cellDependenciesPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cell_deps <- renderPlotly({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_celldeps(achilles, expression_names, data()$gene_symbols, mean_virtual_achilles), tooltip = "text")
        })
      })      
    }
  )
}

cellBinsPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("cell_bins"),  height = "auto")),
    tags$br(),
    fluidRow(tags$strong(plot_cellbins_title), plot_cellbins_legend)
  )
}

cellBinsPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_bins <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "")) #""left blank
        make_cellbins(achilles, expression_names, data()$gene_symbols)
      },
      height = function() length(data()$gene_symbols) * 90 + 80)
    }
  )
}

cellDepsLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("cell_deps_lin"),  height = "auto")),
    tags$br(),
    fluidRow(tags$strong(plot_celllin_title), plot_celllin_legend)
    )
}

cellDepsLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_deps_lin <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          make_lineage(achilles, expression_names, data()$gene_symbols)
        })
      },
      height = 550)
      observeEvent(input$sublin_click, { #event to store the 'click'
      })
      output$cell_deps_sublin <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          make_sublineage(achilles, expression_names, data()$gene_symbols)
        })
      },
      height = 1400)
    }
  )
}

cellDepsSubLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(actionLink(inputId = ns("sublin_click"), "View plot split into sublineages")),
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("sublin_click"), "'] != 0"), 
                     tags$br(),
                     fluidRow(plotOutput(outputId = ns("cell_deps_sublin"),  height = "auto")))
    )
}

cellDepsSubLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$sublin_click, { #event to store the 'click'
      })
      output$cell_deps_sublin <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          make_sublineage(achilles, expression_names, data()$gene_symbols)
        })
      },
      height = 1400)
    }
  )
}




cellAnatogramPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_subcell_exp_plot")))),
    fluidRow(plotOutput(outputId = ns("cellanatogram"))),
  )
    
}
  
cellAnatogramPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_subcell_exp_plot <- renderText({paste0("Sub-cellular expression of ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cellanatogram <- renderPlot({
        validate(
          need(data()$gene_symbols %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, data()$gene_symbols)
      })
    }
  )
}

cellExpressionPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_exp_plot")))),
    fluidRow(plotlyOutput(outputId = ns("cell_exp"))),
  )
}

cellExpressionPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_exp_plot <- renderText({paste0("Cell line expression of ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cell_exp <- renderPlotly({
        validate(
          need(data()$gene_symbols %in% colnames(expression), "")) #""left blank
        ggplotly(make_cellexpression(gene_symbol = data()$gene_symbols), tooltip = c("text"))
      })
    }
  )
}

geneCorrelationPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_genecorr_plot")))),
    fluidRow(plotOutput(outputId = ns("gene_correlations"))),
    tags$br(),
    fluidRow(tags$strong(plot_genecorrelations_title), plot_genecorrelations_legend))
}

geneCorrelationPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_genecorr_plot <- renderText({paste0("Gene correlation plot generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$gene_correlations <- renderPlot({
        validate(
          need(data()$gene_symbols %in% achilles_cor_nest$fav_gene, "No data found for this gene."))
        withProgress(message = 'Working hard under here...', value = 1, {
          make_correlation(achilles_cor_nest, data()$gene_symbols, mean_virtual_achilles, achilles_upper, achilles_lower)
        })
      })      
    }
  )
}

# Exp v. Dep Plots -----
cellexpdepPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cellexpdep_plot")))),
    fluidRow(plotlyOutput(outputId = ns("cellexpdep"))),
    tags$br(),
    fluidRow(tags$strong(plot_expdep_title), plot_expdep_legend))
}

cellexpdepPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cellexpdep_plot <- renderText({paste0("Dependency v. expression plots generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cellexpdep <- renderPlotly({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_expdep(expression_data = expression, celldeps_data = achilles, expression_join = expression_names, gene_symbol = data()$gene_symbols), tooltip = "text")
        })
      })      
    }
  )
}
