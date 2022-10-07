geneNetworkGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    uiOutput(ns("network_graph"))
  )
}

geneNetworkGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_graph <- renderText({paste0("Co-essentiality network graph for ", str_c(data()$content, collapse = ", "))})
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
        make_graph(input = data(), 
                   threshold = rv$threshold, 
                   deg = rv$degree, 
                   corrType = rv$corrType, 
                   cell_line_similarity = "dependency",
                   displayHeight = '80vh', 
                   displayWidth = '100%', 
                   tooltipLink = TRUE) %>% 
          visLegend(position = "right", width = .25, zoom = F) #.5 fixes cut-off, but makes it too wide
      })
      
      output$network_graph <- renderUI({
        sidebarLayout(
          sidebarPanel(sliderInput(inputId = session$ns("deg"),
                                   label = "Filter connections (<)",
                                   value = 2, min = 1, max = 10),
                       sliderInput(inputId = session$ns("threshold"),
                                   label = paste0("# Related ", 
                                                  ifelse(data()$type == "gene",
                                                         "gene", "cell line")),
                                   value =10, min = 10, max = 20),
                       selectInput(inputId = session$ns("corrType"),
                                   label = "Associations",
                                   choices = c("Positive", "Negative", "Positive and Negative"),
                                   selected = "Positive"),
                       actionButton(inputId = session$ns("update"), 
                                    label = "Update", 
                                    width = "100%"),
                       tags$br(),
                       width = 3),
          
          mainPanel(visNetworkOutput(outputId = session$ns("graph"), height = "70vh") %>% # 70vh corresponds to 70% the size of the viewing port
                      withSpinnerColor(plot_type = data()$type), #see shiny_helper.R
                    width = 9)
        )
      })
      
      output$graph <- renderVisNetwork({
        if(data()$type == "gene") {
          shiny::validate(
            need(data()$content %in% master_top_table$fav_gene | 
                   data()$content %in% master_bottom_table$fav_gene,
                 "Not enough data to make a graph.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(data()$content %in% unique(cell_line_dep_sim$cell1_name),
                 "Not enough data to make a graph.")
          )
        }
        networkGraph()
        })
    }
  )
}

geneNetworkGraphExp <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    uiOutput(ns("network_graph"))
  )
}

geneNetworkGraphExpServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_graph <- renderText({paste0("Co-expression network graph for ", str_c(data()$content, collapse = ", "))})
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
        make_graph(input = data(), 
                   threshold = rv$threshold, 
                   deg = rv$degree, 
                   corrType = rv$corrType, 
                   cell_line_similarity = "expression",
                   displayHeight = '80vh', 
                   displayWidth = '100%', 
                   tooltipLink = TRUE) %>% 
          visLegend(position = "right", width = .25, zoom = F) #.5 fixes cut-off, but makes it too wide
      })
      
      output$network_graph <- renderUI({
        sidebarLayout(
          sidebarPanel(sliderInput(inputId = session$ns("deg"),
                                   label = "Filter connections (<)",
                                   value = 2, min = 1, max = 10),
                       sliderInput(inputId = session$ns("threshold"),
                                   label = paste0("# Related ", 
                                                  ifelse(data()$type == "gene",
                                                         "gene", "cell line")),
                                   value =10, min = 10, max = 20),
                       selectInput(inputId = session$ns("corrType"),
                                   label = "Associations",
                                   choices = c("Positive", "Negative", "Positive and Negative"),
                                   selected = "Positive"),
                       actionButton(inputId = session$ns("update"), 
                                    label = "Update", 
                                    width = "100%"),
                       tags$br(),
                       width = 3),
          
          mainPanel(visNetworkOutput(outputId = session$ns("graph"), height = "70vh") %>% # 70vh corresponds to 70% the size of the viewing port
                      withSpinnerColor(plot_type = data()$type), #see shiny_helper.R
                    width = 9)
        )
      })
      
      output$graph <- renderVisNetwork({
        if(data()$type == "gene") {
          shiny::validate(
            need(data()$content %in% master_top_table$fav_gene | 
                   data()$content %in% master_bottom_table$fav_gene,
                 "Not enough data to make a graph.")
          )
        } else if(data()$type == "cell") {
          shiny::validate(
            need(data()$content %in% unique(cell_line_exp_sim$cell1_name),
                 "Not enough data to make a graph.")
          )
        }
        networkGraph()
      })
    }
  )
}

####

compoundNetworkGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_compound_graph")))),
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
      mainPanel(visNetworkOutput(outputId = ns("compound_graph"), height = "70vh") %>% # 70vh corresponds to 70% the size of the viewing port
                  withSpinnerColor(plot_type = "compound"), #see shiny_helper.R
                width = 9)
    )
  )
}

compoundNetworkGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_compound_graph <- renderText({paste0("Network graph for ", str_c(data()$content, collapse = ", "))})
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
      
      networkCompoundGraph <- reactive({
        make_graph(input = data(), 
                   threshold = rv$threshold, 
                   deg = rv$degree, 
                   corrType = rv$corrType, 
                   displayHeight = '80vh', 
                   displayWidth = '100%', 
                   tooltipLink = TRUE) %>% 
          visLegend(position = "right", width = .25, zoom = F) #.5 fixes cut-off, but makes it too wide
      })
      
      output$compound_graph <- renderVisNetwork({
        shiny::validate(
          need(data()$content %in% prism_cor_nest$fav_drug, "No data found."))
        networkCompoundGraph()
      })
    }
    
  )
}

geneBipartiteGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    sidebarLayout(
      sidebarPanel(sliderInput(inputId = ns("threshold"),
                               label = "# Related genes",
                               value =10, min = 10, max = 20),
                   selectInput(inputId = ns("corrType"),
                               label = "Correlations",
                               choices = c("Positive", "Negative", "Positive and Negative"),
                               selected = "Positive"),
                   checkboxInput(inputId = ns("simplify"), 
                                 label = "Simplify",
                                 value = TRUE),
                   actionLink(inputId = ns("censor_click"), "Censor"), # conditional panel for censor list
                   conditionalPanel(condition = paste0("input['", ns("censor_click"), "'] != 0"), 
                                    checkboxGroupInput(
                                      inputId = ns("censor"), 
                                      label = NULL, 
                                      choices = c("Water", "Adenosine triphosphate", "ADP", "NAD", "NADH", "Oxygen"), 
                                      selected = c("Water"))),
                   actionButton(inputId = ns("update"), 
                                label = "Update", 
                                width = "100%"),
                   tags$br(),
                   width = 3),
      mainPanel(visNetworkOutput(outputId = ns("graph")) %>%
                  withSpinnerColor(plot_type = "gene"), #see shiny_helper.R
                width = 9)
    )
  )
}

geneBipartiteGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_graph <- renderText({paste0("Bipartite network graph for ", str_c(data()$content, collapse = ", "))})
      
      #event to store the censor 'click'
      observeEvent(input$censor_click, { 
      })
      #establish reactive value
      rv <- reactiveValues(threshold = 10,
                           corrType = "Positive", 
                           collapsed = TRUE)
      
      #update value upon call
      observeEvent(input$update, {
        rv$threshold <- input$threshold
        rv$corrType <- input$corrType
        rv$collapsed <- input$simplify
        rv$censor <- input$censor
      })
      
      bipartiteNetworkGraph <- reactive({
        make_bipartite_graph(input = data(), 
                             threshold = rv$threshold, 
                             corrType = rv$corrType, 
                             collapsed = rv$collapsed, 
                             censor = rv$censor) %>% 
          visLegend(position = "right", width = .2, zoom = F) #.5 fixes cut-off, but makes it too wide
      })
      
      output$graph <- renderVisNetwork({
        shiny::validate(
          need(data()$content %in% hmdb_proteins$fav_gene, "No data found."))
        bipartiteNetworkGraph()
      })
    }
    
  )
}

compoundBipartiteGraph <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_graph")))),
    fluidRow(visNetworkOutput(outputId = ns("graph"))  %>% 
               withSpinnerColor(plot_type = "compound") #see shiny_helper.R
    )
  )
}

compoundBipartiteGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_graph <- renderText({paste0("Bipartite network graph for ", str_c(data()$content, collapse = ", "))})
      
      bipartiteNetworkGraph <- reactive({
        make_bipartite_graph(input = data()) %>% 
          visLegend(position = "right", width = .2, zoom = F) #.5 fixes cut-off, but makes it too wide
      })
      
      output$graph <- renderVisNetwork({
        shiny::validate(
          need(data()$content %in% hmdb_metabolites$fav_metabolite, "No data found."))
        bipartiteNetworkGraph()
      })
    }
    
  )
}
