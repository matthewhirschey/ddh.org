library(shiny)
library(tidyverse)
library(plotly)
library(networkD3)
library(corrr)
library(here)
library(lubridate)
library(feather)
library(rmarkdown)
library(markdown)

#LOAD DATA-----
read_gene_summary_into_environment <- function(tmp.env) {
  # Read gene_summary saved as RData using: save(gene_summary, file=here::here("data", "gene_summary.RData"))
  load(here::here("data", "gene_summary.RData"), envir=tmp.env)
}
load(file=here::here("data", "achilles.RData"))
load(file=here::here("data", "expression_join.RData"))

#read in datasets saved from depmap_generate_pathways.Rmd
load(file=here::here("data", "master_bottom_table.RData"))
load(file=here::here("data", "master_top_table.RData"))
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

make_enrichment_table <- function(table, gene_symbol) { #master_positive, master_negative
  table %>% 
    dplyr::filter(fav_gene == gene_symbol) %>% 
    unnest(data) %>% 
    select(enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>% 
    arrange(Adjusted.P.value)
}

setup_plot <- function(gene_symbol) {
  #plot setup
  target_achilles <- achilles %>% 
    select(X1, gene_symbol) %>% 
    left_join(expression_join, by = "X1") %>% 
    rename(dep_score = gene_symbol) %>% 
    select(cell_line, lineage, dep_score) 
  
  target_achilles_top <- target_achilles %>% 
    arrange(desc(dep_score)) %>% 
    top_frac(dep_score, n = 10)
  
  target_achilles_bottom <- target_achilles %>% 
    arrange(dep_score) %>% 
    top_frac(dep_score, n = 10)
}
  
make_plot <- function(gene_symbol, op) { #op = cell_bins, cell_deps
  #plot setup
  setup_plot(gene_symbol)
  target_achilles$cell_line <- fct_reorder(target_achilles$cell_line, target_achilles$dep_score, .desc = FALSE)
  
  switch(op, 
  #plot1:cell_bins
  cell_bins = ggplot(target_achilles) +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "lightgray") +
    labs(x = "Dependency Score (binned)") + 
    theme_light(),
    #plot2:cell_deps
  cell_deps = ggplot(target_achilles) +
    geom_point(aes(x = cell_line, y = dep_score), alpha = 0.2) +
    labs(x = "Cell Lines", y = "Dependency Score") +
    geom_hline(yintercept = mean_virtual_achilles) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    geom_point(data = target_achilles_top, aes(x = cell_line, y = dep_score), color = "red") +
    geom_point(data = target_achilles_bottom, aes(x = cell_line, y = dep_score), color = "red") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
    NULL, 
  stop("Specify your plot!")
  )
}

setup_graph <- function(gene_symbol, threshold = 10) {
  dep_network <- tibble()
  
 #find top and bottom correlations for fav_gene
  dep_top <- make_top_table(gene_symbol) %>% 
    slice(1:threshold)
  
  dep_bottom <- make_bottom_table(gene_symbol) %>% 
    slice(1:threshold) #limit for visualization?
  
  #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
  related_genes <- dep_top %>% 
    bind_rows(dep_bottom) %>% 
    dplyr::pull(gene)
  
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in related_genes){
    message("Getting correlations related to ", gene_symbol, ", including ", i)
    dep_top_related <- achilles_cor %>% 
      focus(i) %>% 
      arrange(desc(.[[2]])) %>% #use column index
      filter(.[[2]] > achilles_upper) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i) %>% 
      rename(y = rowname, r2 = i) %>% 
      select(x, y, r2) %>% 
      slice(1:threshold) #limit for visualization?
    
    dep_bottom_related <- achilles_cor %>% 
      focus(i) %>% 
      arrange(.[[2]]) %>% #use column index
      filter(.[[2]] < achilles_lower) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i) %>% 
      rename(y = rowname, r2 = i) %>% 
      select(x, y, r2) %>% 
      slice(1:threshold) #limit for visualization?
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>% 
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>% 
      bind_rows(dep_related)
  }
  return(dep_network)
}

make_graph <- function(dep_network = dep_network, deg = 2) { #change to make_graph
  #setup graph 
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  nodes <-  as_tibble(graph_network) %>% 
    rowid_to_column("id") %>% 
    mutate(degree = degree(graph_network), 
           group = 1) %>% 
    arrange(desc(degree))
  
  links <- graph_network %>% 
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>% 
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>% 
    as.data.frame
  
  # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1
  
  forceNetwork(Links = links_filtered, Nodes = nodes_filtered, Source = "from", Target ="to", NodeID = "name", Group = "group", zoom = TRUE, bounded = TRUE, opacityNoHover = 100)
}

render_depmap_report_to_file <- function(file, gene_symbol) {
  tmp.env <- environment()
  read_gene_summary_into_environment(tmp.env)
  render_complete_report(file, gene_symbol, tmp.env)
}

render_complete_report <- function (file, gene_symbol, tmp.env) {
  fav_gene_summary <- tmp.env$gene_summary %>%
    filter(approved_symbol == gene_symbol)
  rmarkdown::render("report_depmap_app.rmd", output_file = paste0(gene_symbol, '_report.pdf'))
  
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
             textInput(inputId = "gene_symbol", label = "Enter gene symbol", value ='TP53'), 
             actionButton(inputId = "go", label = "Generate"), 
             hr(), 
             uiOutput("gene_summary")),
    navbarMenu(title = "Cell Dependencies",
               tabPanel("Plots",
                        fluidRow("text"),
                        fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                             plotOutput(outputId = "plotdeps"), 
                                             plotOutput(outputId = "plotbins"))),
                        fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                             "text", 
                                             "text"))),
               tabPanel("Table",
                        fluidRow("text"),
                        fluidRow(dataTableOutput(outputId = "target_achilles_top")))
               ),
    navbarMenu(title = "Similar",
               tabPanel("Genes",
                        fluidRow("text"),
                        fluidRow(dataTableOutput(outputId = "dep_top"))),
               tabPanel("Pathways",
                        fluidRow("text"),
                        fluidRow(dataTableOutput(outputId = "pos_enrich")))),
    navbarMenu(title = "Dissimilar",
               tabPanel("Genes",
                        fluidRow("text"),
                        fluidRow(dataTableOutput(outputId = "dep_bottom"))),
               tabPanel("Pathways",
                        fluidRow("text"),
                        fluidRow(dataTableOutput(outputId = "neg_enrich")))),
    tabPanel("Graph", 
             #sidebarLayout( #UNCOMMENT THIS SECTION WHEN input$deg is working
               #sidebarPanel(numericInput(inputId = "deg", 
                                        #label = "Minimum # of Connections", 
                                        #value = 2, min = 1, max = 50)),
               #conditionalPanel(condition = "input.gene_symbol == NULL", 
               #                 p("Enter a gene symbol to generate a graph")),
               mainPanel(forceNetworkOutput(outputId = "graph"))#)
    ),
    tabPanel("Methods", 
             includeMarkdown(here::here("code", "methods.md"))
             ),
    tabPanel("Download Report",
             h2("Report Generator"), 
             "To generate a report, click on the button below",
             br(),
             downloadButton(outputId = "report", label = "Download report"))
  )
)

#SERVER-----
server <- function(input, output, session) {
  data <- eventReactive(input$go, {
    str_to_upper(input$gene_symbol)})
  
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
  output$plotdeps <- renderPlot(
    make_plot(data(), "cell_deps")
  )
  output$plotbins <- renderPlot(
    make_plot(data(), "cell_bins")
  )
  output$target_achilles_top <- renderDataTable(
    setup_plot(data()),
    target_achilles_top, 
    options = list(pageLength = 12)
  )
  output$pos_enrich <- renderDataTable(
    make_enrichment_table(master_positive, data())
  )
  output$neg_enrich <- renderDataTable(
    make_enrichment_table(master_negative, data())
  )
  output$graph <- renderForceNetwork(
    make_graph(setup_graph(data()))
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
