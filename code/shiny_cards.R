cellDependenciesPlotDash <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_depdash_plot")))),
    fluidRow(
      div(plotOutput(outputId = ns("cell_depdash"))), style = "box-shadow: rgba(100, 100, 111, 0.2) 0px 7px 29px 0px;") #https://getcssscan.com/css-box-shadow-examples
  )
}

cellDependenciesPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_depdash_plot <- renderText({paste0(str_c(data()$gene_symbols, collapse = ", "), " dependency plot")})
      output$cell_depdash <- renderPlot({
        validate(need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
         make_celldeps(achilles, expression_names, data()$gene_symbols, mean_virtual_achilles)
      })      
    }
  )
}