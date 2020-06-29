library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(networkD3)
library(corrr)
library(here)
library(lubridate)
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
library(gganatogram)

render_report_in_background <- FALSE
if (supportsMulticore()) {
  plan(multicore)
  render_report_in_background <- TRUE
}

#LOAD DATA-----
#read current release information
source(here::here("code", "current_release.R"))

#read data from create_*.R
gene_summary <- readRDS(here::here("data", paste0(release, "_gene_summary.Rds")))
pathways <- readRDS(here::here("data", paste0(release, "_pathways.Rds")))

#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here("data", paste0(release, "_achilles.Rds")))
expression_join <- readRDS(file=here::here("data", paste0(release, "_expression_join.Rds")))

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here("data", paste0(release, "_sd_threshold.Rds")))
achilles_lower <- readRDS(file = here::here("data", paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here("data", paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here("data", paste0(release, "_sd_virtual_achilles.Rds")))

#read data from generate_depmap_tables & pathways.R
master_bottom_table <- readRDS(file=here::here("data", paste0(release, "_master_bottom_table.Rds")))
master_top_table <- readRDS(file=here::here("data", paste0(release, "_master_top_table.Rds")))
master_positive <- readRDS(file=here::here("data", paste0(release, "_master_positive.Rds")))
master_negative <- readRDS(file=here::here("data", paste0(release, "_master_negative.Rds")))
surprise_genes <- readRDS(file=here::here("data", paste0(release, "_surprise_genes.Rds")))
censor_genes <- readRDS(file=here::here("data", paste0(release, "_censor_genes.Rds")))

#read data from generate_subcell_data.R
subcell <- readRDS(file=here::here("data", paste0(release, "_subcell.Rds")))

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_tables.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "fun_graphs.R"))
source(here::here("code", "fun_reports.R"))

#shiny functions
search_panel <- function() {
  searchInput(
    inputId = "gene_or_pathway",
    placeholder = "genes, pathways, or GO number", #search
    btnSearch = icon("search")
  )
}

surprise <- function(table = master_top_table, surprise_vec = surprise_genes) {
  gene_symbol <- sample(surprise_vec, 1)
  gene_symbol_url <- paste0("?show=detail&content=gene&symbol=", gene_symbol)
  return(gene_symbol_url)
}

query_result_row <- function(row) {
  if (row$contents == 'gene') {
    gene_query_result_row(row)
  } else if (row$contents == 'pathway') {
      pathway_query_result_row(row)
  } else {
    gene_list_query_result_row(row)
  }
}

gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
   h4(
     tags$strong("Gene:"),
     tags$a(title, href=paste0("?show=detail&content=gene&symbol=", gene_summary_row["approved_symbol"]))
   ),
   div(tags$strong("Aka:"), gene_summary_row["aka"]),
   div(tags$strong("Entrez ID:"), gene_summary_row["ncbi_gene_id"]),
   hr()
  )
}

pathway_query_result_row <- function(row) {
  pathways_row <- row$data
  gene_symbols <- lapply(pathways_row$data, function(x) { paste(x$gene, collapse=', ') })
  title <- paste0(pathways_row$pathway, " (GO:", pathways_row$go, ")")
  list(
    h4(
      tags$strong("Pathway:"),
      tags$a(title, href=paste0("?show=detail&content=pathway&go=", pathways_row$go))
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
    ),
    hr()
  )
}

gene_list_query_result_row <- function(row) {
  gene_summary_rows <- row$data
  title <- row$key
  
  known_gene_symbols <- gene_summary_rows %>% 
    filter(known == TRUE) %>%
    pull(approved_symbol)
  has_known_gene_symbols <- !is_empty(known_gene_symbols)
  
  unknown_gene_symbols <- gene_summary_rows %>% 
    filter(known == FALSE) %>%
    pull(approved_symbol)
  has_unknown_gene_symbols <- !is_empty(unknown_gene_symbols)
  
  known_gene_symbols_tags <- NULL
  if (has_known_gene_symbols) {
    gene_query_param <- paste0("custom_gene_list=", paste(known_gene_symbols, collapse=","))
    href <- paste0("?show=detail&content=custom_gene_list&", gene_query_param)
    known_gene_symbols_tags <- list(
      tags$h6("Known Gene Symbols"),
      tags$a(paste(known_gene_symbols, collapse=", "), href=href)
    )
  }

  unknown_gene_symbols_tags <- NULL
  if (has_unknown_gene_symbols) {
    unknown_gene_symbols_tags <- list(
      tags$h6("Unknown Gene Symbols"),
      tags$div(paste(unknown_gene_symbols, collapse=", "))
    )
  }

  list(
    h4(
      tags$strong("Custom Gene List")#,
      #tags$span(title)
    ),
    known_gene_symbols_tags,
    unknown_gene_symbols_tags,
    hr()
  )
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
      tags$dt("Gene Summary"), tags$dd(gene_summary$entrez_summary)
    ), 
    hr(),
    plotOutput(outputId = "cellanatogram"), 
    hr(), 
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

pathway_summary_details <- function(pathways_row) {
  gene_symbols <- lapply(pathways_row$data, function(x) { paste(x$gene, collapse=', ') })
  title <- paste0("Pathway: ", pathways_row$pathway, " (GO:", pathways_row$go, ")")
  list(
    h4(
      tags$strong(title),
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
      tags$dt("Pathway Description"), 
      tags$dd(pathways_row$def)
    ),
    hr(),
    plotOutput(outputId = "cellanatogram"), 
    hr(), 
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

gene_list_summary_details <- function(custom_gene_list) {
  gene_symbols <- paste(custom_gene_list, collapse=', ')
  title <- paste0("Custom Gene List") #, gene_symbols
  list(
    h4(
      tags$strong(title),
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
    ),
    hr(),
    plotOutput(outputId = "cellanatogram"),
    hr(),
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

# renders details about gene
gene_summary_ui <- function(gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    gene_summary_row <- gene_summary %>%
      filter(approved_symbol == gene_symbol)
      result <- gene_summary_details(gene_summary_row)
  }
  result
}

pathway_summary_ui <- function(pathway_go) {
  pathway_row <- pathways %>%
    filter(go == pathway_go)
  pathway_summary_details(pathway_row)
}

gene_list_summary_ui <- function(custom_gene_list) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- gene_summary %>%
    filter(approved_symbol %in% custom_gene_list) %>%
    pull(approved_symbol)
  gene_list_summary_details(valid_gene_symbols)
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

#UI------
### HEAD
head_tags <- tags$head(includeHTML("gtag.html"),includeScript("returnClick.js"))

### universal elements
main_title <- HTML('<a href="." style="color:black;">Data-Driven Hypothesis</a>')

search_tab_panel <- div(
  search_panel()
)

### HOME (landing) PAGE
home_page <- tagList(
  head_tags,
  HTML('<center><img src="hex_ddh.png"></center>'),
  tags$div(
    tags$br(),
    HTML('<center>Data-driven hypothesis is a resource developed by the <a href="http://www.hirscheylab.org" style="color:black;">Hirschey Lab</a> for predicting functional relationships for thousands of genes across the human genome.</center>'), 
    tags$br(),
    tags$br()),
  HTML("<center>"),
  search_panel(), 
  actionLink(inputId = "example_click", "See example searches"), 
  ", ", 
  actionLink(inputId = "pathway_click", "browse the pathways"), 
  ", or",
  htmlOutput("get_lucky", inline = TRUE),
  HTML("</center>"),
  conditionalPanel(condition = 'input.example_click == 0',
                   ""),
  conditionalPanel(condition = 'input.example_click != 0', 
                   tags$br(),
                   h4("Examples"),
                   HTML('<h5>Search for</h5>
                        <ul>
                        <li>A single gene, such as <a href="?show=detail&content=gene&symbol=TP53">TP53</a> or <a href="?show=detail&content=gene&symbol=BRCA1">BRCA1</a></li>
                        <li>A pathway name, such as <a href="?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="?show=detail&content=pathway&go=0006695">Cholesterol Biosynthetic Process</a></li>
                        <li>The Gene Ontology biological process identifier, such as <a href="?show=search&query=1901989">1901989</a>, which will find <a href="?show=detail&content=pathway&go=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                        <li>A custom list of genes (separated by commas), such as <a href="?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="?show=detail&content=custom_gene_list&custom_gene_list=BRCA1,BRCA2">a custom gene list</a></li>
                       </ul>')
                   ),
  conditionalPanel(condition = 'input.pathway_click == 0',
                   ""),
  conditionalPanel(condition = 'input.pathway_click != 0', 
                   tags$br(),
                   h4("GO Biological Processes"),
                   dataTableOutput(outputId = "pathway_table")), 
)

### SEARCH PAGE
search_page <- tagList(
  head_tags,
  navbarPage(title = main_title, windowTitle = "Data-Driven Hypothesis | A Hirschey Lab Resource"),
  div(search_panel(), style="float: right"),
  h3(textOutput("search_title")),
  div(div(h3("Results", class="panel-title"), class="panel-heading"),
      div(uiOutput("genes_search_result"), class="panel-body"),
      class="bg-info panel panel-default"
  )
)

### DETAILS PAGE: shows either gene or pathway data
detail_page <- fluidPage(
  head_tags,
  navbarPage(title = main_title, windowTitle = "Data-Driven Hypothesis | A Hirschey Lab Resource", 
    tabPanel("Home",
             div(search_panel(), style="float: right"),
             uiOutput("detail_summary")
    ),
    navbarMenu(title = "Cell Dependencies",
               tabPanel("Plots",
                        fluidRow(h4(textOutput("text_cell_dep_plot"))),
                        fluidRow(plotlyOutput(outputId = "cell_deps")),
                        tags$br(),
                        fluidRow(tags$strong(plot_celldeps_title), plot_celldeps_legend),
                        tags$hr(),
                        fluidRow(plotlyOutput(outputId = "cell_bins")),
                        tags$br(),
                        fluidRow(tags$strong(plot_cellbins_title), plot_cellbins_legend),
                        tags$hr(),
                        fluidRow(plotlyOutput(outputId = "cell_deps_lin")),
                        tags$br(),
                        fluidRow(tags$strong(plot_celllin_title), plot_celllin_legend, 
                        actionLink(inputId = "sublin_click", "View sublineage plot below")), #add conditional panel for subtype
                        tags$br(), 
                        conditionalPanel(condition = 'input.sublin_click == 0',
                                         ""),
                        conditionalPanel(condition = 'input.sublin_click != 0', 
                                         tags$br(),
                                         fluidRow(plotlyOutput(outputId = "cell_deps_sublin")))
                        ),
               tabPanel("Table",
                        fluidRow(h4(textOutput("text_cell_dep_table"))),
                        fluidRow(dataTableOutput(outputId = "target_achilles")))
    ),
    navbarMenu(title = "Similar",
               tabPanel("Genes",
                        fluidRow(h4(textOutput("text_dep_top"))),
                        fluidRow(
                          column(6, checkboxGroupInput(inputId = "vars_dep_top", 
                                                    "Select columns:",
                                                    c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                    selected = c("Z Score", "Co-publication Count"), 
                                                    inline = TRUE)),
                          column(6, fluidRow(sliderInput(inputId = "num_sim_genes",
                                                           "Censor genes with more than n associations:",
                                                           min = 100,
                                                           max = 1000,
                                                           value = 1000,
                                                           step = 100)), 
                                      fluidRow(column(3, actionButton(inputId = "censor", "Submit")),
                                               column(3, actionButton(inputId = "reset", "Reset"))))
                          ),
                        hr(),
                        fluidRow(dataTableOutput(outputId = "dep_top"))),
               tabPanel("Pathways",
                        fluidRow(h4(textOutput("text_pos_enrich"))),
                        fluidRow(dataTableOutput(outputId = "pos_enrich")))),
    navbarMenu(title = "Dissimilar",
               tabPanel("Genes",
                        fluidRow(h4(textOutput("text_dep_bottom"))),
                        fluidRow(checkboxGroupInput(inputId = "vars_dep_bottom", 
                                                    "Select columns:",
                                                    c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                    selected = c("Z Score", "Co-publication Count"), 
                                                    inline = TRUE)),
                        fluidRow(dataTableOutput(outputId = "dep_bottom"))),
               tabPanel("Pathways",
                        fluidRow(h4(textOutput("text_neg_enrich"))),
                        fluidRow(dataTableOutput(outputId = "neg_enrich")))),
    tabPanel("Graph",
             sidebarLayout(
               sidebarPanel(sliderInput(inputId = "deg",
                                        label = "Filter \nConnections (<)",
                                        value = 2, min = 1, max = 10),
                            sliderInput(inputId = "threshold",
                                        label = "n related \nGenes",
                                        value =10, min = 10, max = 20),
               actionButton(inputId = "update", 
                            label = "Update"),
               width = 3),
             mainPanel(forceNetworkOutput(outputId = "graph"),
                       width = 9)
    )
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

detail_page_visible <- function() {
  getQueryString()$show == 'detail'
}

gene_callback <- function(input, output, session) {
  data <- reactive({
    if (detail_page_visible()) {
      content <- getQueryString()$content
      if (content == 'gene') {
        gene_symbol <- getQueryString()$symbol
      } else if (content == 'pathway') {
        pathway_go <- getQueryString()$go
        pathway_row <- pathways %>%
          filter(go == pathway_go)
        pathway_row$data[[1]]$gene
      } else {
        custom_gene_list <- getQueryString()$custom_gene_list
        c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
      }
    }
  })

  #reactive 'get lucky' to dynamically update within a single session
  output$get_lucky <- renderUI({
    tags$a(href = surprise(), "get lucky")
  })

  #When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    save_data(input)
  })

  #examples & pathways
  observeEvent(input$example_click, { #event to store the 'click'
  })
  observeEvent(input$pathway_click, { #event to store the 'click'
  })
  
  output$pathway_table <- DT::renderDataTable({
    DT::datatable(make_pathway_table(pathways) %>% dplyr::rename(Pathway = pathway, GO = go, Genes = genes), 
                  options = list(pageLength = 10))
  })
  
  output$cellanatogram <- renderPlot({
    validate(
      need(data() %in% subcell$gene_name, "No subcellular location data for this gene."))
    make_cellanatogram(subcell, data())
  })
  
  output$cellanatogram_table <- DT::renderDataTable({
    validate(
      need(data() %in% subcell$gene_name, ""))
    DT::datatable(make_cellanatogram_table(subcell, data()), 
                  options = list(pageLength = 10))
  })
    
  output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", str_c(data(), collapse = ", "))})
  output$text_cell_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data(), collapse = ", "))})
  output$text_dep_top <- renderText({paste0("Genes with similar dependencies as ", str_c(data(), collapse = ", "))})
  output$text_pos_enrich <- renderText({paste0("Pathways of genes with similar dependencies as ", str_c(data(), collapse = ", "))})
  output$text_dep_bottom <- renderText({paste0("Genes with inverse dependencies as ", str_c(data(), collapse = ", "))})
  output$text_neg_enrich <- renderText({paste0("Pathways of genes with inverse dependencies as ", str_c(data(), collapse = ", "))})

  output$detail_summary <- renderUI({
    if (detail_page_visible()) {
      content <- getQueryString()$content
      if (content == 'gene') {
        gene_summary_ui(data())
      } else if (content == 'pathway') {
        pathway_summary_ui(getQueryString()$go)
      } else {
        custom_gene_list <- getQueryString()$custom_gene_list
        gene_list_summary_ui(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
      }
    }
    # render details about the gene symbol or pathway user chose
  })
  #censor reactive values
  censor_status <- reactiveValues(choice = FALSE, 
                                  num_sim_genes = 1000)
  
  observeEvent(input$censor, {
    censor_status$choice <- TRUE
    censor_status$num_sim_genes <- input$num_sim_genes
  })
  
  observeEvent(input$reset, {
    censor_status$choice <- FALSE
    censor_status$num_sim_genes <- 1000
  })
  
  output$dep_top <- DT::renderDataTable({
    validate(
      need(data() %in% master_top_table$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_top_table(master_top_table, data()) %>% 
        dplyr::mutate(link = paste0("<center><a href='?show=detail&content=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
        dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_top) %>%
        censor(censor_genes, censor_status$choice, censor_status$num_sim_genes),
      escape = FALSE,
      options = list(pageLength = 25))
  })
  output$dep_bottom <- DT::renderDataTable({
    validate(
      need(data() %in% master_bottom_table$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_bottom_table(master_bottom_table, data()) %>% 
        dplyr::mutate(link = paste0("<center><a href='?show=detail&content=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
        dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_bottom),
      escape = FALSE,
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
    ggplotly(make_cellbins(achilles, expression_join, data()), tooltip = c("text"))
  })
  output$cell_deps_lin <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    withProgress(message = 'Wait for it...', value = 1, {
      ggplotly(make_lineage(achilles, expression_join, data()), tooltip = "text")
    })
  })
  observeEvent(input$sublin_click, { #event to store the 'click'
  })
  output$cell_deps_sublin <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    withProgress(message = 'Wait for it...', value = 1, {
      ggplotly(make_sublineage(achilles, expression_join, data()), height = 1100, tooltip = "text")
    })
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
      make_enrichment_top(master_positive, data()),
    options = list(pageLength = 25))
  })
  output$neg_enrich <- DT::renderDataTable({
    validate(
      need(data() %in% master_negative$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_enrichment_bottom(master_negative, data()),
      options = list(pageLength = 25))
  })

  #establish reactive value
  rv <- reactiveValues(degree = 2, 
                       threshold = 10)

  #update value upon call
  observeEvent(input$update, {
    rv$degree <- input$deg
    rv$threshold <- input$threshold
  })
  
  output$graph <- renderForceNetwork({
    validate(
      need(data() %in% colnames(achilles), "No data found."))
    withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
    make_graph(master_top_table, master_bottom_table, data(), threshold = rv$threshold, deg = rv$degree)
    })
  })
  output$report <- downloadHandler(
    # create pdf report
    filename = function() {
      type = getQueryString()$content
      if(type == "gene"){
      paste0(data(), "_ddh.pdf")
      } else if (type == "pathway") {
        go <- getQueryString()$go
      paste0("go_", go, "_ddh.pdf")
      } else {
        paste0("custom_", data(), "_ddh.pdf")
      }
        },
    content = function(file) {
      gene_symbol <- data() # reactive data must be read outside of a future
      progress_bar <- Progress$new()
      progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
      if (render_report_in_background) {
        result <- future({
          render_report_to_file(file, 
                                gene_symbol, 
                                type = getQueryString()$content,
                                summary1 = gene_summary, 
                                summary2 = pathways,
                                cellbins_data = achilles, 
                                expression_data = expression_join, 
                                celldeps_data = achilles,
                                mean = mean_virtual_achilles,
                                cellanatogram_data = subcell,
                                toptable_data = master_top_table, 
                                bottomtable_data = master_bottom_table,
                                enrichmenttop_data = master_positive, 
                                enrichmentbottom_data = master_negative, 
                                achilles_data = achilles)
        })
        finally(result, function(){
          progress_bar$close()
        })
      } else {
        render_report_to_file(file, 
                              gene_symbol, 
                              type = getQueryString()$content,
                              summary1 = gene_summary, 
                              summary2 = pathways,
                              cellbins_data = achilles, 
                              expression_data = expression_join, 
                              celldeps_data = achilles,
                              mean = mean_virtual_achilles,
                              cellanatogram_data = subcell,
                              toptable_data = master_top_table, 
                              bottomtable_data = master_bottom_table,
                              enrichmenttop_data = master_positive, 
                              enrichmentbottom_data = master_negative, 
                              achilles_data = achilles)
        progress_bar$close()
      }
    }
  )
}

# Create output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageContent")
  )
)

pages_ui <- list(
  home=home_page,
  search=search_page,
  detail=detail_page
)

search_callback <- function(input, output, session) {
  output$search_title <- renderText({
    query <- getQueryString()
    paste0("Search results for '", query$query, "'")
  })
  output$genes_search_result <- renderUI({
    query <- getQueryString()
    if (grepl(',', query$query)) {
      query_results_table <- gene_list_query_results_table(gene_summary, query$query)
    } else {
      query_results_table <- gene_or_pathway_query_results_table(gene_summary, pathways, query$query)
    }
    if (nrow(query_results_table) > 0) {
      apply(query_results_table, 1, query_result_row)
    }
    else {
      "No results found."
    }
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
  
  search_callback(input, output, session)
  gene_callback(input, output, session)
})

shinyApp(ui, server)
