library(shiny)
library(tidyverse)
library(corrr)
library(here)
library(lubridate)
library(feather)
library(rmarkdown)

#LOAD DATA-----~20sec
read_gene_summary_into_environment <- function(tmp.env) {
  # Read gene_summary saved as RData using: save(gene_summary, file=here::here("data", "gene_summary.RData"))
  load(here::here("data", "gene_summary.RData"), envir=tmp.env)
}

# Read achilles saved as RData using: save(achilles, file=here::here("data", "achilles.RData"))
achilles <- load(here::here("data", "achilles.RData"))
#achilles_long <- achilles %>% pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

# Read achilles_cor saved as RData using: save(achilles_cor, file=here::here("data", "achilles_cor.RData"))
achilles_cor <- load(here::here("data", "achilles_cor.RData"))

# Read expression_id saved as RData using: save(expression_id, file=here::here("data", "expression_id.RData")), need this to get cell line names
#expression_id <- load(here::here("data", "expression_id.Rdata"))
#expression_join <- expression_id %>% 
#  rename(X1 = dep_map_id) %>% 
#  select(X1, stripped_cell_line_name, lineage)


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

#UI------
ui <- fluidPage(
  navbarPage(
    title = "Data-Driven Hypothesis", 
    tabPanel("Home", 
             textInput(inputId = "gene_symbol", label = "Enter gene symbol", value ='BRCA1'), 
             actionButton(inputId = "go", label = "Generate"), 
             hr(), 
             uiOutput("gene_summary")),
    tabPanel("Cell Dependencies", 
             splitLayout(
             plotOutput(outputId = "plot2"), 
             plotOutput(outputId = "plot1"))),
    tabPanel("Similar", 
             dataTableOutput(outputId = "dep_top")),
    tabPanel("Dissimilar", 
             dataTableOutput(outputId = "dep_bottom")),
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
    achilles_cor %>% 
    focus(data()) %>% 
    arrange(desc(.[[2]])) %>% #use column index
    filter(.[[2]] > 0.1963) %>% #formerly top_n(20), but changed to mean +/- 3sd
    rename(approved_symbol = rowname) %>% 
    left_join(gene_summary, by = "approved_symbol") %>% 
    select(approved_symbol, approved_name, data()) %>% 
    rename(gene = approved_symbol, name = approved_name, r2 = data())
  )
  output$dep_bottom <- renderDataTable(
    achilles_cor %>% 
      focus(data()) %>% 
      arrange(.[[2]]) %>% #use column index
      filter(.[[2]] < -0.1899) %>% #formerly top_n(20), but changed to mean +/- 3sd
      rename(approved_symbol = rowname) %>% 
      left_join(gene_summary, by = "approved_symbol") %>% 
      select(approved_symbol, approved_name, data()) %>% 
      rename(gene = approved_symbol, name = approved_name, r2 = data())
  )
  output$plot1 <- renderPlot(
    achilles_long %>% 
      filter(gene == data()) %>% 
      left_join(expression_join, by = "X1") %>% 
      select(stripped_cell_line_name, lineage, dep_score) %>% 
      ggplot() +
      geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "lightgray") +
      labs(x = "Dependency Score (binned)") + 
      theme_light()
  )
  output$plot2 <- renderPlot(
    achilles_long %>% 
      filter(gene == data()) %>% 
      left_join(expression_join, by = "X1") %>% 
      select(stripped_cell_line_name, lineage, dep_score) %>% 
      ggplot() +
      geom_point(aes(x = fct_rev(fct_reorder(target_achilles$stripped_cell_line_name, target_achilles$dep_score, .desc = TRUE)), y = dep_score)) +
      labs(x = "Cell Lines", y = "Dependency Score") +
      geom_hline(yintercept = mean_virtual_achilles) +
      geom_hline(yintercept = achilles_upper, linetype="dashed") +
      geom_hline(yintercept = achilles_lower, linetype="dashed") +
      geom_hline(yintercept = 0) +
      theme_light() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
      NULL
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
