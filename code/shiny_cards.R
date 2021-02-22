box_shadow <- "box-shadow: rgba(100, 100, 111, 0.2) 0px 7px 29px 0px; height: 450px;"

#cell dependencies plot-----
cellDependenciesPlotDash <- function(id) {
  ns <- NS(id)
  div(plotOutput(outputId = ns("cell_depdash")), style = box_shadow) # https://getcssscan.com/css-box-shadow-examples
}

cellDependenciesPlotDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_depdash <- renderPlot({
        validate(need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        make_celldeps(achilles, expression_names, data()$gene_symbols, mean_virtual_achilles)
      })      
    }
  )
}
#cell dependencies table-----
cellDependenciesTableDash <- function(id) {
  ns <- NS(id)
  div(gt_output(outputId = ns("deptabledash")), style = box_shadow) #https://getcssscan.com/css-box-shadow-examples
}

cellDependenciesTableDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$deptabledash <- render_gt({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        gt::gt(make_top_table(master_top_table, data()$gene_symbols) %>% mutate("Rank" = row_number()) %>% select("Rank", "Gene", "R^2") %>% slice(1:8)) %>% 
          tab_header(
            title = "Co-essential genes"
          )
      })      
    }
  )
}
#cell dependencies graph-----
cellDependenciesGraphDash <- function(id) {
  ns <- NS(id)
  div(visNetworkOutput(outputId = ns("depgraphdash")), style = box_shadow) #https://getcssscan.com/css-box-shadow-examples
}

cellDependenciesGraphDashServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$depgraphdash <- renderVisNetwork({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found."))
        withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
          make_graph(master_top_table, master_bottom_table, data()$gene_symbols, threshold = 10, deg = 2, corrType = "Positive")
        })
      })
    })
}