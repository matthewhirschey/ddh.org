library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(networkD3)
library(corrr)
library(here)
library(lubridate)
library(feather)
library(rmarkdown)
library(markdown)
library(tidygraph)
library(ggraph)
library(viridis)
library(cowplot)
library(plotly)
library(DT)
library(future)
library(promises)

render_report_in_background <- FALSE
if (supportsMulticore()) {
  plan(multicore)
  render_report_in_background <- TRUE
}

#LOAD DATA-----
#read current release information
source(here::here("code", "current_release.R"))

#read data from create_gene_summary.R
gene_summary <- readRDS(here::here("data", "gene_summary.Rds"))
     
#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here("data", paste0(release, "_achilles.Rds")))
expression_join <- readRDS(file=here::here("data", paste0(release, "_expression_join.Rds")))

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here("data", "sd_threshold.Rds"))
achilles_lower <- readRDS(file = here::here("data", "achilles_lower.Rds"))
achilles_upper <- readRDS(file = here::here("data", "achilles_upper.Rds"))
mean_virtual_achilles <- readRDS(file = here::here("data", "mean_virtual_achilles.Rds"))
sd_virtual_achilles <- readRDS(file = here::here("data", "sd_virtual_achilles.Rds"))

#read data from generate_depmap_pathways.R
master_bottom_table <- readRDS(file=here::here("data", "master_bottom_table.Rds"))
master_top_table <- readRDS(file=here::here("data", "master_top_table.Rds"))
master_positive <- readRDS(file=here::here("data", "master_positive.Rds"))
master_negative <- readRDS(file=here::here("data", "master_negative.Rds"))

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_tables.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "fun_graphs.R"))

#shiny functions
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

gene_search_results <- function(gene_summary_row) {
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    div(
      tags$a(title, href=paste0("?show=gene&symbol=", gene_summary_row["approved_symbol"]))
      ),
    div(tags$strong("Aka:"), gene_summary_row["aka"]),
    div(tags$strong("Entrez ID:"), gene_summary_row["ncbi_gene_id"]),
    hr()
  )
}

gene_summary_not_found <- function(gene_symbol) {
  if(sum(str_detect(gene_summary$aka, paste0("(?<![:alnum:])", gene_symbol, "(?![:alnum:]|\\-)"))) > 0) {
    result <- tagList(h4(paste0("Gene symbol ", gene_symbol, " not found. Make sure to use the 'Official' gene symbol. Did you mean any of these: ", str_c(pull(gene_summary[str_which(gene_summary$aka, gene_symbol), 2]), collapse = ", "), "?")))
  } else {
    result <- tagList(h4(paste0("Gene symbol ", gene_symbol, " not found. Please make sure this is the 'Official' gene symbol and not an alias.")))
  }  
  return(result)
}

# renders 'not found' or details about gene
gene_summary_ui <- function(gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    gene_summary_row <- gene_summary %>%
      filter(approved_symbol == gene_symbol)
    if (dim(gene_summary_row)[1] == 0) {
      result <- gene_summary_not_found(gene_symbol)
    } else {
      title <- paste0(gene_summary_row$approved_symbol, ": ", gene_summary_row$approved_name)
      result <- gene_summary_details(gene_summary_row)
    }
  }
  result
}

render_report_to_file <- function(file, gene_symbol) {
  if (gene_symbol %in% colnames(achilles)) {
    src <- normalizePath('report_depmap_app.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_depmap_app.Rmd', overwrite = TRUE)
    out <- render_complete_report(file, gene_symbol)
    file.rename(out, file)
  } else {
    src <- normalizePath('report_dummy_depmap.Rmd')
    
    # temporarily switch to the temp dir, in case you do not have write
    # permission to the current working directory
    
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    
    file.copy(src, 'report_dummy_depmap.Rmd', overwrite = TRUE)
    out <- render_dummy_report(file, gene_symbol)
    file.rename(out, file)
  }
}

render_complete_report <- function (file, gene_symbol) {
  fav_gene_summary <- gene_summary %>%
    filter(approved_symbol == gene_symbol)
  p1 <- make_celldeps(achilles, expression_join, gene_symbol, mean_virtual_achilles)
  p2 <- make_cellbins(achilles, expression_join, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles, expression_join, gene_symbol) %>% slice(1:10)
  target_achilles_top <- make_achilles_table(achilles, expression_join, gene_symbol) %>% dplyr::arrange(desc(`Dependency Score`)) %>%  slice(1:10)
  dep_top <- make_top_table(master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_table(master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_table(master_negative, gene_symbol)
  graph_report <- make_graph_report(master_top_table, master_bottom_table, gene_symbol)
  rmarkdown::render("report_depmap_app.Rmd", output_file = file)

}
render_dummy_report <- function (file, gene_symbol) {
  fav_gene_summary <- gene_summary %>%
    filter(approved_symbol == gene_symbol)
  rmarkdown::render("report_dummy_depmap.Rmd", output_file = file)
}

# Define the fields we want to save from the form
fields <- c("first_name", "last_name", "email")

save_data <- function(input) {
  # put variables in a data frame
  data <- data.frame(matrix(nrow=1,ncol=0))
  for (x in fields) {
    data[[x]] <- input[[x]]
  }
  data$submit_time <- date()
  
  # Create a unique file name
  file_name <- sprintf(
    "%s_%s.csv",
    as.integer(Sys.time()), 
    digest::digest(data) #gives it a unique name
  )
  
  #create dir
  directory_path <- here::here("user-data")
  dir.create(file.path(directory_path), showWarnings = FALSE)
  # Write the file to the local system as csv without column headers for ease of use
  write_csv(data, path=file.path(directory_path, file_name), col_names=FALSE)
}

search_panel <- function() {
  searchInput(
    inputId = "gene_or_pathway",
    placeholder = "TP53",
    btnSearch = icon("search")
  )
}

head_tags <- tags$head(includeHTML("gtag.html"),includeScript("returnClick.js"))

main_title <- "Data-Driven Hypothesis"

search_tab_panel <- div(
  search_panel()
)

home_page <- tagList(
  head_tags,
  navbarPage(title = main_title),
  h4("Enter gene symbol"),
  search_panel(),
)

search_page <- tagList(
  head_tags,
  navbarPage(title = main_title),
  div(search_panel(), style="float: right"),
  h3(textOutput("search_title")),
  div(div(h3("Gene", class="panel-title"), class="panel-heading"),
      div(uiOutput("genes_search_result"), class="panel-body"),
      class="bg-info panel panel-default"
  )
)


#UI------
gene_page <- fluidPage(
  head_tags,
  navbarPage(
    title = main_title,
    tabPanel("Home",
             "Data-driven hypothesis is a resource developed by the", a(href = "http://www.hirscheylab.org", "Hirschey Lab"), "to functionally map human genes. A typical use case starts by querying a gene, identifying genes that share similar patterns or behaviors across several measures, in order to discover novel genes in established processes or new functions for well-studied genes.",
             hr(),
             uiOutput("gene_summary"),
             div(search_panel(), style="float: right")
    ),
    navbarMenu(title = "Cell Dependencies",
               tabPanel("Plots",
                         fluidRow(h4(textOutput("text_cell_dep_plot"))),
                         fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                              plotlyOutput(outputId = "cell_deps"),
                                              plotlyOutput(outputId = "cell_bins"))),
                         fluidRow(htmlOutput(outputId = "plot_text"))),
               tabPanel("Table",
                        fluidRow(h4(textOutput("text_cell_dep_table"))),
                        fluidRow(dataTableOutput(outputId = "target_achilles")))
    ),
    navbarMenu(title = "Similar",
               tabPanel("Genes",
<<<<<<< HEAD
                        fluidRow(h4(textOutput("text_dep_top"))),
                        fluidRow(dataTableOutput(outputId = "dep_top"))),
=======
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of similar genes"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_dep_top"))),
                                         fluidRow(checkboxGroupInput(inputId = "vars_dep_top", 
                                                                     "Select columns:",
                                                                     c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                                     selected = c("Z Score", "Co-publication Count"), 
                                                                     inline = TRUE)),
                                         fluidRow(dataTableOutput(outputId = "dep_top")))),
>>>>>>> master
               tabPanel("Pathways",
                        fluidRow(h4(textOutput("text_pos_enrich"))),
                        fluidRow(dataTableOutput(outputId = "pos_enrich")))),
    navbarMenu(title = "Dissimilar",
               tabPanel("Genes",
<<<<<<< HEAD
                        fluidRow(h4(textOutput("text_dep_bottom"))),
                        fluidRow(dataTableOutput(outputId = "dep_bottom"))),
=======
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of dissimilar genes"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_dep_bottom"))),
                                         fluidRow(checkboxGroupInput(inputId = "vars_dep_bottom", 
                                                                     "Select columns:",
                                                                     c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                                     selected = c("Z Score", "Co-publication Count"), 
                                                                     inline = TRUE)),
                                         fluidRow(dataTableOutput(outputId = "dep_bottom")))),
>>>>>>> master
               tabPanel("Pathways",
                        fluidRow(h4(textOutput("text_neg_enrich"))),
                        fluidRow(dataTableOutput(outputId = "neg_enrich")))),
    tabPanel("Graph",
             #sidebarLayout( #UNCOMMENT THIS SECTION WHEN input$deg is working
             #sidebarPanel(numericInput(inputId = "deg",
             #label = "Minimum # of Connections",
             #value = 2, min = 1, max = 50)),
             mainPanel(forceNetworkOutput(outputId = "graph"))
             #uncomment this parenthesis)
    ),
    tabPanel("Methods",
             includeHTML(here::here("code","methods.html"))),
    tabPanel("Download Report",
              h2("Report Generator"),
              conditionalPanel(condition = 'input.submit == 0',
                               "Please enter your name and email address to download a report", 
                               textInput("first_name", "First Name", ""), 
                               textInput("last_name", "Last Name", ""), 
                               textInput("email", "Email Address", ""), 
                               actionButton(inputId = "submit", 
                                            label = "Enter")),
              conditionalPanel(condition = 'input.submit != 0', 
                               "To generate a report, click on the button below",
                               br(),
                               downloadButton(outputId = "report", label = "Download report")))
  )
)

#SERVER-----

gene_callback <- function(input, output, session) {
  data <- reactive({
    gene_symbol <- getQueryString()$symbol
    if_else(str_detect(gene_symbol, "orf"), gene_symbol, str_to_upper(gene_symbol))
  })

  # When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    save_data(input)
  })

  output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", data())})
  output$text_cell_dep_table <- renderText({paste0("Dependency table generated for ", data())})
  output$text_dep_top <- renderText({paste0("Genes with similar dependencies as ", data())})
  output$text_pos_enrich <- renderText({paste0("Pathways of genes with similar dependencies as ", data())})
  output$text_dep_bottom <- renderText({paste0("Genes with inverse dependencies as ", data())})
  output$text_neg_enrich <- renderText({paste0("Pathways of genes with inverse dependencies as ", data())})

  output$gene_summary <- renderUI({
    # render details about the gene symbol user entered
    gene_summary_ui(data())
  })
  output$dep_top <- DT::renderDataTable({
    validate(
      need(data() %in% master_top_table$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_top_table(master_top_table, data()) %>% dplyr::select("Gene", "Name", input$vars_dep_top), 
      options = list(pageLength = 25))
  })
  output$dep_bottom <- DT::renderDataTable({
    validate(
      need(data() %in% master_bottom_table$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_bottom_table(master_bottom_table, data()) %>% dplyr::select("Gene", "Name", input$vars_dep_bottom),
    options = list(pageLength = 25))
  })
  output$cell_deps <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    withProgress(message = 'Wait for it...', value = 1, {
      ggplotly(make_celldeps(achilles, expression_join, data(), mean_virtual_achilles), tooltip = "text")
    })
  })
  output$cell_bins <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "")) #""left blank
    ggplotly(make_cellbins(achilles, expression_join, data()))
  })
  output$plot_text <- renderUI({
    validate(
      need(data() %in% colnames(achilles), "")) #""left blank
    HTML("<p></p><p><b>Left:</b> Cell Line Dependency Curve. Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essentail within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out. <b>Right:</b> Histogram of Binned Dependency Scores. Dependency scores across all cell lines for queried gene partitioned into 0.25 unit bins. Shape of the histogram curve reveals overall influence of queried gene on cellular fitness</p>")
    })
  output$target_achilles <- DT::renderDataTable({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    make_achilles_table(achilles, expression_join, data())
  })
  output$pos_enrich <- DT::renderDataTable({
    validate(
      need(data() %in% master_positive$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_enrichment_table(master_positive, data()),
    options = list(pageLength = 25))
  })
  output$neg_enrich <- DT::renderDataTable({
    validate(
      need(data() %in% master_negative$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_enrichment_table(master_negative, data()),
      options = list(pageLength = 25))
  })
  output$graph <- renderForceNetwork({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
    make_graph(master_top_table, master_bottom_table, data())
    })
  })
  output$report <- downloadHandler(
    # create pdf report
    filename = function() {paste0(data(), "_ddh.pdf")},
    content = function(file) {
      gene_symbol <- data() # reactive data must be read outside of a future
      progress_bar <- Progress$new()
      progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
      if (render_report_in_background) {
        result <- future({
          render_report_to_file(file, gene_symbol)
        })
        finally(result, function(){
          progress_bar$close()
        })
      } else {
        render_report_to_file(file, gene_symbol)
        progress_bar$close()
      }
    }
  )
}

# Creat output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageContent")
  )
)

pages_ui <- list(
  home=home_page,
  search=search_page,
  gene=gene_page
)

search_callback <- function(input, output, session) {
  output$search_title <- renderText({
    query <- getQueryString()
    paste("Search results for", query$query)
  })
  output$genes_search_result <- renderUI({
    query <- getQueryString()
    selected_genes <- gene_summary %>%
      filter(str_detect(approved_symbol, paste0("^", query$query))) %>%
      head(10)
    apply(selected_genes, 1, gene_search_results)
  })
}

server <- shinyServer(function(input, output, session) {
  output$pageContent <- renderUI({
    query <- getQueryString()
    show_page <- query$show
    if (is.null(show_page)) {
      show_page <- "home"
    }
    pages_ui[show_page]
  })
  
  observeEvent(input$gene_or_pathway_search, { 
    updateQueryString(paste0("?show=search&query=", input$gene_or_pathway), mode="push")
  })
  
  # the below item is problematic for the back button
  search_callback(input, output, session)
  gene_callback(input, output, session)
})

shinyApp(ui, server)
