library(tictoc)
library(pryr)
library(shiny)
source(here::here("code", "install_libraries.R"))

#LOAD DATA-----
#read params
source(here::here("code", "app_params.R"))
source(here::here("code", "generate_colors.R"))
source(here::here("code", "app_data.R"))

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_helper.R"), local = TRUE)
source(here::here("code", "fun_tables.R"), local = TRUE)
source(here::here("code", "fun_plots.R"), local = TRUE)
source(here::here("code", "fun_graphs.R"), local = TRUE)
source(here::here("code", "fun_reports.R"), local = TRUE)
source(here::here("code", "fun_text.R"), local = TRUE)
source(here::here("code", "fun_structures.R"), local = TRUE)

#SHINY FUNCTIONS-----
source(here::here("code", "shiny_tables.R"), local = TRUE)
source(here::here("code", "shiny_plots.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)
source(here::here("code", "shiny_reports.R"), local = TRUE)
source(here::here("code", "shiny_text.R"), local = TRUE)
source(here::here("code", "shiny_cards.R"), local = TRUE)
source(here::here("code", "shiny_structures.R"), local = TRUE)
source(here::here("code", "shiny_download.R"), local = TRUE)

# data function to pass to the module
my_gene_data <- function() {
    list(
        type="gene",
        subtype="gene",
        id="ROCK2",
        gene_symbols=c("ROCK2")
    )
}

testPage <- function (id) {
    ns <- NS(id)
    # module render function to test
    div(
        geneTitle(ns("title_var"))
    )
}

testPageServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            # module server function to test
            geneTitleServer("title_var", my_gene_data)
        }
    )
}

ui <- fluidPage(
    titlePanel("Test App"),
    uiOutput("pageUI")
)

server <- function(input, output) {
    output$pageUI <- renderUI({
        testPage("test")
    })
    testPageServer("test")
}

# Run the application 
shinyApp(ui = ui, server = server)
