library(shiny)
library(tidyverse)
library(corrr)
library(here)
library(lubridate)
library(feather)
library(rmarkdown)

#LOAD DATA-----
read_gene_summary_into_environment <- function(tmp.env) {
  # Read gene_summary saved as RData using: save(gene_summary, file=here::here("data", "gene_summary.RData"))
  load(here::here("data", "gene_summary.RData"), envir=tmp.env)
}

#read in datasets saved from depmap_generate_pathways.Rmd
load(file=here::here("data", "master_bottom_table.RData"))
load(file=here::here("data", "master_top_table.RData"))
load(file=here::here("data", "master_plots.RData"))
load(file=here::here("data", "master_positive.RData"))
load(file=here::here("data", "master_negative.RData"))


#FUNCTIONS-----
gene_summary_details <- function(gene_summary) {
  title <- paste0(gene_summary$approved_symbol, ": ", gene_summary$approved_name)
  tagList(
    h3(title),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene"), tags$dd(gene_summary$approved_symbol),
      tags$dt("Name"), tags$dd(gene_summary$approved_name),
      tags$dt("aka"), tags$dd(gene_summary$aka),
      tags$dt("Entrez ID"), tags$dd(gene_summary$ncbi_gene_id), 
      tags$dt("Gene Summary"), tags$dd(gene_summary$entrez_summary)
    )
  )
}

# renders 'not found' or details about gene
gene_summary_ui <- function(gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    tmp.env <- environment()
    read_gene_summary_into_environment(tmp.env)
    gene_summary_row <- tmp.env$gene_summary %>%
      filter(approved_symbol == gene_symbol)
    if (dim(gene_summary_row)[1] == 0) {
      result <- tagList(
        h4(paste0("Gene symbol \"", gene_symbol, "\" not found. Please make sure this is the 'Official' gene symbol and not an alias."))
      )
    } else {
      title <- paste0(gene_summary_row$approved_symbol, ": ", gene_summary_row$approved_name)
      result <- gene_summary_details(gene_summary_row)
    }
  }
  result
}

make_top_table <- function(gene_symbol) {
  master_top_table %>% 
    dplyr::filter(fav_gene == gene_symbol) %>% 
    unnest(data) %>% 
    select(-fav_gene) %>% 
    arrange(desc(r2))
}

make_bottom_table <- function(gene_symbol) {
  master_bottom_table %>% 
    dplyr::filter(fav_gene == gene_symbol) %>% 
    unnest(data) %>% 
    select(-fav_gene) %>% 
    arrange(r2)
}
make_enrichment_table <- function(table, gene_symbol) {
  table %>% 
    filter(fav_gene == gene_symbol) %>% 
    unnest(data) %>% 
    select(-fav_gene) %>% 
    arrange(Adjusted.P.value)
}

make_plot1 <- function(gene_symbol) {
  master_plots %>% 
    filter(fav_gene == gene_symbol) %>% 
    pluck("plot1")
}
make_plot2 <- function(gene_symbol) {
  master_plots %>% 
    filter(fav_gene == gene_symbol) %>% 
    pluck("plot2")
}
render_depmap_report_to_file <- function(file, gene_symbol) {
  tmp.env <- environment()
  read_gene_summary_into_environment(tmp.env)
  render_dummy_report(file, gene_symbol, tmp.env)
}

render_dummy_report <- function (file, gene_symbol, tmp.env) {
  fav_gene_summary <- tmp.env$gene_summary %>%
    filter(approved_symbol == gene_symbol)
  rmarkdown::render("report_dummy_depmap.rmd", output_file = file)
}


#UI------
ui <- fluidPage(
  navbarPage(
    title = "Data-Driven Hypothesis", 
    tabPanel("Home", 
             "text",
             hr(),
             textInput(inputId = "gene_symbol", label = "Enter gene symbol", value ='SETD1B'), 
             actionButton(inputId = "go", label = "Generate"), 
             hr(), 
             uiOutput("gene_summary")),
    tabPanel("Cell Dependencies", 
             fluidRow(splitLayout(cellWidths = c("50%", "50%"),
              plotOutput(outputId = "plot2"), 
              plotOutput(outputId = "plot1"))),
             fluidRow(splitLayout(cellWidths = c("50%", "50%"),
              "text", 
              "text"))),
    tabPanel("Similar",
             fluidRow("text"),
             fluidRow(dataTableOutput(outputId = "dep_top")), 
             fluidRow(dataTableOutput(outputId = "pos_enrich"))),
    tabPanel("Dissimilar", 
             fluidRow(dataTableOutput(outputId = "dep_bottom")), 
             fluidRow(dataTableOutput(outputId = "neg_enrich"))),
    tabPanel("Methods", 
             h3("Methods"), 
             "description"),
    tabPanel("Generate Report",
             h2("Report Generator"), 
             "To generate a report, click on the button below",
             br(),
             downloadButton("report", "Generate report"))
  )
)

#SERVER-----
server <- function(input, output, session) {
  data <- eventReactive(input$go, {input$gene_symbol})
  output$gene_summary <- renderUI({
    # render details about the gene symbol user entered
    gene_summary_ui(data())
  })
  output$dep_top <- renderDataTable(
    make_top_table(data())
  )
  output$dep_bottom <- renderDataTable(
    make_bottom_table(data())
  )
  output$plot1 <- renderPlot(
    make_plot1(data())
  )
  output$plot2 <- renderPlot(
    make_plot2(data())
  )
  output$pos_enrich <- renderDataTable(
    make_enrichment_table(master_positive, data())
  )
  output$neg_enrich <- renderDataTable(
    make_enrichment_table(master_negative, data())
  )
  output$report <- downloadHandler(
    # create pdf depmap report
    filename = paste0(data(), "_depmap.pdf"),
    content = function(file) {
      render_depmap_report_to_file(file, data())
    }
  )
}

shinyApp(ui, server)
