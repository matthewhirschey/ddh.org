geneNetworkGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    sidebarLayout(
      sidebarPanel(sliderInput(inputId = ns("deg"),
                               label = "Filter \nConnections (<)",
                               value = 2, min = 1, max = 10),
                   sliderInput(inputId = ns("threshold"),
                               label = "n related \nGenes",
                               value =10, min = 10, max = 20),
                   selectInput(inputId = ns("corrType"),
                                 label = "Types of Correlations to Include:",
                               choices = c("Positive and Negative", "Positive", "Negative"),
                               selected = "Positive and Negative"),
                   actionButton(inputId = ns("update"), 
                                label = "Update", 
                                width = "100%"),
                   tags$br(),
                   downloadLink(outputId = ns("downloadNetwork"), 'Download network'),
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
                           corrType = "Positive and Negative" )
      
      #update value upon call
      observeEvent(input$update, {
        rv$degree <- input$deg
        rv$threshold <- input$threshold
        rv$corrType <- input$corrType
      })
      
      networkGraph <- reactive({
        make_graph(master_top_table, master_bottom_table, data()$gene_symbols, threshold = rv$threshold, deg = rv$degree, corrType = rv$corrType, displayHeight = '80vh', displayWidth = '100%', tooltipLink = TRUE) %>% 
          visLegend(position = "right", width = .25, zoom = F)
      })
      
      output$graph <- renderVisNetwork({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found."))
        withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
         networkGraph()
        })
      })
      
      output$downloadNetwork <- downloadHandler(
        filename = function() {
          if(length(data()$gene_symbols) == 1){         
            paste('network_', data()$gene_symbols, '_', Sys.Date(), '.html', sep='')
          }
          else{
            exportName <- ""
            for(geneName in data()$gene_symbols){
              exportName <- paste0(geneName, "_", exportName)
            }
           paste('custom_network_', exportName, Sys.Date(), '.html', sep='')
          }
        },
        content = function(file) {
          make_graph(master_top_table, master_bottom_table, data()$gene_symbols, threshold = rv$threshold, deg = rv$degree, corrType = rv$corrType, displayHeight = '80vh', displayWidth = '100%', tooltipLink = FALSE) %>% 
            visLegend(position = "right", width = .25, zoom = F) %>%
            visSave(file)
        }
      )

      
    }
      
  )
}

