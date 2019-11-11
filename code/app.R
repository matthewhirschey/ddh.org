library(shiny)
library(tidyverse)
library(here)
library(lubridate)
library(feather)
library(rmarkdown)


read_gene_summary_into_environment <- function(tmp.env) {
  # Read gene_summary saved as RData using: save(gene_summary, file=here::here("data", "gene_summary.RData"))
  load(here::here("data", "gene_summary.RData"), envir=tmp.env)
}

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
    ),
    downloadButton("report", "Generate report"),
  )
}

# renders 'not found' or details about gene and download button
gene_summary_ui <- function(gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    tmp.env <- environment()
    read_gene_summary_into_environment(tmp.env)
    gene_summary_row <- tmp.env$gene_summary %>%
      filter(approved_symbol == gene_symbol)
    if (dim(gene_summary_row)[1] == 0) {
      result <- tagList(
        h4(paste0("Gene symbol \"", gene_symbol, "\" not found."))
      )
    } else {
      title <- paste0(gene_summary_row$approved_symbol, ": ", gene_summary_row$approved_name)
      result <- gene_summary_details(gene_summary_row)
    }
  }
  result
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

ui <- fluidPage(
    titlePanel("Depmap"),
    sidebarLayout(
        sidebarPanel(
            textInput("gene_symbol", "Enter gene symbol", "", placeholder='BRCA1')
        ),
        mainPanel(
            uiOutput("gene_summary")
        )
    )
)

server <- function(input, output, session) {
    output$gene_summary <- renderUI({
        # render details about the gene symbol user entered
        gene_summary_ui(input$gene_symbol)
    })
    output$report <- downloadHandler(
        # create pdf depmap report
        filename = paste0(input$gene_symbol, "_depmap.pdf"),
        content = function(file) {
          render_depmap_report_to_file(file, input$gene_symbol)
        }
    )
}

shinyApp(ui, server)
