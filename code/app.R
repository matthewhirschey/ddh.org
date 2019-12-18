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
library(tidygraph)
library(ggraph)
library(viridis)
library(cowplot)
library(plotly)

#LOAD DATA-----
#read current release information
source(here::here("code", "current_release.R"))

#read data from creat_gene_summary.R
read_gene_summary_into_environment <- function(tmp.env) {
  # Read gene_summary saved as RData using: save(gene_summary, file=here::here("data", "gene_summary.RData"))
  load(here::here("data", "gene_summary.RData"), envir=tmp.env)
}
#read data from generate_depmap_data.R
load(file=here::here("data", paste0(release, "_achilles.RData")))
load(file=here::here("data", paste0(release, "_achilles_cor.RData")))
load(file=here::here("data", paste0(release, "_expression_join.RData")))

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here("data", "sd_threshold.Rds"))
achilles_lower <- readRDS(file = here::here("data", "achilles_lower.Rds"))
achilles_upper <- readRDS(file = here::here("data", "achilles_upper.Rds"))
mean_virtual_achilles <- readRDS(file = here::here("data", "mean_virtual_achilles.Rds"))
sd_virtual_achilles <- readRDS(file = here::here("data", "sd_virtual_achilles.Rds"))

#read data from generate_depmap_pathways.R
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
      ifelse(sum(str_detect(gene_summary$aka, paste0("^", gene_symbol))) > 0,
             result <- tagList(h4(paste0("Found ", sum(str_detect(gene_summary$aka, paste0("^", gene_symbol))), " entries. Make sure to use the 'Official' gene symbol."))),
             result <- tagList(h4(paste0("Gene symbol ", gene_symbol, " not found. Please make sure this is the 'Official' gene symbol and not an alias."))))
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
    arrange(desc(r2)) %>%
    rename("Gene" = "gene", "Name" = "name", "R^2" = "r2")
}

make_bottom_table <- function(gene_symbol) {
  master_bottom_table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    unnest(data) %>%
    select(-fav_gene) %>%
    arrange(r2) %>%
    rename("Gene" = "gene", "Name" = "name", "R^2" = "r2")
}

make_enrichment_table <- function(table, gene_symbol) { #master_positive, master_negative
  table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    unnest(data) %>%
    select(enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    arrange(Adjusted.P.value) %>%
    rename("Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(gene_symbol) {
  target_achilles <- achilles %>%
    select(X1, gene_symbol) %>%
    left_join(expression_join, by = "X1") %>%
    rename(dep_score = gene_symbol) %>%
    select(cell_line, lineage, dep_score) %>%
    arrange(dep_score) %>%
    rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Dependency Score" = "dep_score")
  return(target_achilles)
}

make_cellbins <- function(gene_symbol) {
  achilles %>% #plot setup
    select(X1, gene_symbol) %>%
    left_join(expression_join, by = "X1") %>%
    rename(dep_score = gene_symbol) %>%
    select(cell_line, lineage, dep_score) %>%
    arrange(dep_score) %>%
    ggplot() +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "lightgray") +
    labs(x = "Dependency Score (binned)") +
    theme_cowplot()
}

make_celldeps <- function(gene_symbol) {
  achilles %>% #plot setup
    select(X1, gene_symbol) %>%
    left_join(expression_join, by = "X1") %>%
    rename(dep_score = gene_symbol) %>%
    select(cell_line, lineage, dep_score) %>%
    arrange(dep_score) %>%
    ggplot() +
    geom_point(aes(x = fct_reorder(cell_line, dep_score, .desc = FALSE), y = dep_score, text = paste0("Cell Line: ", cell_line)), alpha = 0.2) +
    labs(x = "Cell Lines", y = "Dependency Score") +
    geom_hline(yintercept = mean_virtual_achilles) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    theme_cowplot() +
    #geom_point(data = target_achilles_top, aes(x = cell_line, y = dep_score), color = "red") +
    #geom_point(data = target_achilles_bottom, aes(x = cell_line, y = dep_score), color = "red") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
    NULL
}

make_graph <- function(gene_symbol, threshold = 10, deg = 2) {
  dep_network <- tibble()

 #find top and bottom correlations for fav_gene
  dep_top <- make_top_table(gene_symbol) %>%
    slice(1:threshold)

  dep_bottom <- make_bottom_table(gene_symbol) %>%
    slice(1:threshold) #limit for visualization?

  #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
  related_genes <- dep_top %>%
    bind_rows(dep_bottom) %>%
    dplyr::pull("Gene")

  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in related_genes){
    message("Getting correlations from ", i, " related to ", gene_symbol)
    dep_top_related <- achilles_cor %>%
      focus(i) %>%
      arrange(desc(.[[2]])) %>% #use column index
      filter(.[[2]] > achilles_upper) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i, origin = "pos") %>%
      rename(y = rowname, r2 = i) %>%
      select(x, y, r2, origin) %>%
      slice(1:threshold) #limit for visualization?

    dep_bottom_related <- achilles_cor %>%
      focus(i) %>%
      arrange(.[[2]]) %>% #use column index
      filter(.[[2]] < achilles_lower) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i, origin = "neg") %>%
      rename(y = rowname, r2 = i) %>%
      select(x, y, r2, origin) %>%
      slice(1:threshold) #limit for visualization?

    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)

    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  nodes <-  as_tibble(graph_network) %>%
    rowid_to_column("id") %>%
    mutate(degree = igraph::degree(graph_network),
           group = case_when(name %in% gene_symbol == TRUE ~ "query",
                             name %in% dep_top$Gene ~ "pos",
                             name %in% dep_bottom$Gene ~ "neg",
                             TRUE ~ "connected")) %>%
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

  #Pos #74D055, Neg #3A568C, Q #FDE825, C #450D53
  node_color <- 'd3.scaleOrdinal(["#74D055", "#3A568C", "#FDE825", "#450D53"])'

  forceNetwork(Links = links_filtered, Nodes = nodes_filtered, Source = "from", Target ="to", NodeID = "name", Group = "group", zoom = TRUE, bounded = TRUE, opacityNoHover = 100, Nodesize = "degree", colourScale = node_color)
  }

make_graph_report <- function(gene_symbol, threshold = 10, deg = 2) {
  dep_network <- tibble()
  
  #find top and bottom correlations for fav_gene
  dep_top <- make_top_table(gene_symbol) %>%
    slice(1:threshold)
  
  dep_bottom <- make_bottom_table(gene_symbol) %>%
    slice(1:threshold) #limit for visualization?
  
  #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
  related_genes <- dep_top %>%
    bind_rows(dep_bottom) %>%
    dplyr::pull("Gene")
  
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in related_genes){
    message("Getting correlations from ", i, " related to ", gene_symbol)
    dep_top_related <- achilles_cor %>%
      focus(i) %>%
      arrange(desc(.[[2]])) %>% #use column index
      filter(.[[2]] > achilles_upper) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i, origin = "pos") %>%
      rename(y = rowname, r2 = i) %>%
      select(x, y, r2, origin) %>%
      slice(1:threshold) #limit for visualization?
    
    dep_bottom_related <- achilles_cor %>%
      focus(i) %>%
      arrange(.[[2]]) %>% #use column index
      filter(.[[2]] < achilles_lower) %>% #formerly top_n(20), but changed to mean +/- 3sd
      mutate(x = i, origin = "neg") %>%
      rename(y = rowname, r2 = i) %>%
      select(x, y, r2, origin) %>%
      slice(1:threshold) #limit for visualization?
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  nodes <-  as_tibble(graph_network) %>%
    rowid_to_column("id") %>%
    mutate(degree = igraph::degree(graph_network),
           group = case_when(name %in% gene_symbol == TRUE ~ "query",
                             name %in% dep_top$Gene ~ "pos",
                             name %in% dep_bottom$Gene ~ "neg",
                             TRUE ~ "connected")) %>%
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

  links_filtered$from <- match(links_filtered$from, nodes_filtered$id)
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id)

  graph_network_ggraph <- tidygraph::tbl_graph(nodes = nodes_filtered, edges = links_filtered)

  graph_network_ggraph %>%
    ggraph::ggraph(layout = "auto") +
    geom_edge_fan(aes(edge_width = abs(r2)), alpha = 0.3) +
    geom_node_point(aes(size = degree, color = group), alpha = 0.7) +
    geom_node_label(aes(label = name), repel = TRUE) +
    scale_colour_viridis(discrete = TRUE, name = "Group", labels = c("Query", "Positive", "Negative", "Connected")) +
    theme_graph(base_family = 'Helvetica')
}

render_report_to_file <- function(file, gene_symbol) {
  tmp.env <- environment()
  read_gene_summary_into_environment(tmp.env)
  src <- normalizePath('report_depmap_app.Rmd')

  # temporarily switch to the temp dir, in case you do not have write
  # permission to the current working directory

  owd <- setwd(tempdir())
  on.exit(setwd(owd))

  file.copy(src, 'report_depmap_app.Rmd', overwrite = TRUE)
  out <- render_complete_report(file, gene_symbol, tmp.env)
  file.rename(out, file)
}

render_complete_report <- function (file, gene_symbol, tmp.env) { #how to leverage tmp.env?
  fav_gene_summary <- tmp.env$gene_summary %>%
    filter(approved_symbol == gene_symbol)
  p1 <- make_celldeps(gene_symbol)
  p2 <- make_cellbins(gene_symbol)
  target_achilles_bottom <- make_achilles_table(gene_symbol) %>% slice(1:10)
  target_achilles_top <- make_achilles_table(gene_symbol) %>% dplyr::arrange(desc(`Dependency Score`)) %>%  slice(1:10)
  dep_top <- make_top_table(gene_symbol)
  flat_top_complete <- make_enrichment_table(master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(gene_symbol)
  flat_bottom_complete <- make_enrichment_table(master_negative, gene_symbol)
  graph_report <- make_graph_report(gene_symbol)
  rmarkdown::render("report_depmap_app.rmd", output_file = file)

}
render_dummy_report <- function (file, gene_symbol, tmp.env) {
  fav_gene_summary <- tmp.env$gene_summary %>%
    filter(approved_symbol == gene_symbol)
  rmarkdown::render("report_dummy_depmap.rmd", output_file = file)
}


#UI------
ui <- fluidPage(
  tags$head(includeHTML("gtag.html")),
  navbarPage(
    title = "Data-Driven Hypothesis",
    tabPanel("Home",
             "Data-driven hypothesis is a resource developed by the", a(href = "http://www.hirscheylab.org", "Hirschey Lab"), "to functionally map human genes. A typical use case starts by querying a gene, identifying genes that share similar patterns or behaviors across several measures, in order to discover new novel genes in established processes or new funtions for well-studied genes.",
             hr(),
             textInput(inputId = "gene_symbol", label = "Enter gene symbol", value ='TP53'),
             actionButton(inputId = "go", label = "Generate"),
             hr(),
             uiOutput("gene_summary")),
    navbarMenu(title = "Cell Dependencies",
               tabPanel("Plots",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate dependency plots"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_cell_dep_plot"))),
                                         fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                                                              plotlyOutput(outputId = "cell_deps"),
                                                              plotlyOutput(outputId = "cell_bins"))),
                                         HTML("<p></p><p><b>Left:</b> Cell Line Dependency Curve. Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essentail within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out. <b>Right:</b> Histogram of Binned Dependency Scores. Dependency scores across all cell lines for queried gene partitioned into 0.25 unit bins. Shape of the histogram curve reveals overall influence of queried gene on cellular fitness</p>"))),
               tabPanel("Table",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of dependent cell line"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_cell_dep_table"))),
                                         fluidRow(dataTableOutput(outputId = "target_achilles"))))
    ),
    navbarMenu(title = "Similar",
               tabPanel("Genes",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of similar genes"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_dep_top"))),
                                         fluidRow(dataTableOutput(outputId = "dep_top")))),
               tabPanel("Pathways",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of similar pathways"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_pos_enrich"))),
                                         fluidRow(dataTableOutput(outputId = "pos_enrich"))))),
    navbarMenu(title = "Dissimilar",
               tabPanel("Genes",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of dissimilar genes"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_dep_bottom"))),
                                         fluidRow(dataTableOutput(outputId = "dep_bottom")))),
               tabPanel("Pathways",
                        conditionalPanel(condition = 'input.go == 0',
                                         "Enter a gene symbol to generate a table of dissimilar pathways"),
                        conditionalPanel(condition = 'input.go != 0', 
                                         fluidRow(h4(textOutput("text_neg_enrich"))),
                                         fluidRow(dataTableOutput(outputId = "neg_enrich"))))),
    tabPanel("Graph",
             #sidebarLayout( #UNCOMMENT THIS SECTION WHEN input$deg is working
             #sidebarPanel(numericInput(inputId = "deg",
             #label = "Minimum # of Connections",
             #value = 2, min = 1, max = 50)),
             conditionalPanel(condition = 'input.go == 0',
                              "Enter a gene symbol to generate a graph"),
             conditionalPanel(condition = 'input.go != 0', 
                              mainPanel(forceNetworkOutput(outputId = "graph")))
             #uncomment this parenthesis)
    ),
    tabPanel("Methods",
             includeMarkdown("methods.md")
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
  output$dep_top <- renderDataTable(
    make_top_table(data())
  )
  output$dep_bottom <- renderDataTable(
    make_bottom_table(data())
  )
  output$cell_deps <- renderPlotly(
    ggplotly(make_celldeps(data()), tooltip = "text")
  )
  output$cell_bins <- renderPlotly(
    ggplotly(make_cellbins(data()))
  )
  output$target_achilles <- renderDataTable(
    make_achilles_table(data()),
    options = list(pageLength = 12)
  )
  output$pos_enrich <- renderDataTable(
    make_enrichment_table(master_positive, data())
  )
  output$neg_enrich <- renderDataTable(
    make_enrichment_table(master_negative, data())
  )
  output$graph <- renderForceNetwork({
    withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
    make_graph(data())
    })
  })
  output$report <- downloadHandler(
    # create pdf report
    filename = function() {paste0(data(), "_ddh.pdf")},
    content = function(file) {
      render_report_to_file(file, data())
    }
  )
}

shinyApp(ui, server)
