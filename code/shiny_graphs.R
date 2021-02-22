geneNetworkGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    sidebarLayout(
      sidebarPanel(sliderInput(inputId = ns("deg"),
                               label = "Filter connections (<)",
                               value = 2, min = 1, max = 10),
                   sliderInput(inputId = ns("threshold"),
                               label = "# Related genes",
                               value =10, min = 10, max = 20),
                   selectInput(inputId = ns("corrType"),
                                 label = "Correlations",
                               choices = c("Positive", "Negative", "Positive and Negative"),
                               selected = "Positive"),
                   actionButton(inputId = ns("update"), 
                                label = "Update", 
                                width = "100%"),
                   tags$br(),
                   width = 3),
      mainPanel(visNetworkOutput(outputId = ns("graph"), height = "70vh"),# 70vh corresponds to 70% the size of the viewing port
                width = 9)
    )
  )
}
  
geneNetworkGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_graph <- renderText({paste0("Interactive Network graph for ", str_c(data()$gene_symbols, collapse = ", "))})
      #establish reactive value
      rv <- reactiveValues(degree = 2, 
                           threshold = 10,
                           corrType = "Positive" )
      
      #update value upon call
      observeEvent(input$update, {
        rv$degree <- input$deg
        rv$threshold <- input$threshold
        rv$corrType <- input$corrType
      })
      
      networkGraph <- reactive({
        make_graph(master_top_table, master_bottom_table, data()$gene_symbols, threshold = rv$threshold, deg = rv$degree, corrType = rv$corrType, displayHeight = '80vh', displayWidth = '100%', tooltipLink = TRUE) %>% 
          visLegend(position = "right", width = .5, zoom = F)
      })
      
      output$graph <- renderVisNetwork({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found."))
        withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
         networkGraph()
        })
      })
    }
      
  )
}

