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

make_top_table <- function(gene_symbol) {
  master_top_table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-fav_gene) %>%
    dplyr::mutate(r2 = round(r2, 3)) %>% 
    dplyr::arrange(desc(r2)) %>%
    dplyr::rename("Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_bottom_table <- function(gene_symbol) {
  master_bottom_table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-fav_gene) %>%
    dplyr::mutate(r2 = round(r2, 3)) %>% 
    dplyr::arrange(r2) %>%
    dplyr::rename("Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_enrichment_table <- function(table, gene_symbol) { #master_positive, master_negative
  table %>%
    dplyr::filter(fav_gene == gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(gene_symbol) {
  target_achilles <- achilles %>%
    dplyr::select(X1, gene_symbol) %>%
    dplyr::left_join(expression_join, by = "X1") %>%
    dplyr::rename(dep_score = gene_symbol) %>%
    dplyr::select(cell_line, lineage, dep_score) %>%
    dplyr::mutate(dep_score = round(dep_score, 3)) %>% 
    dplyr::arrange(dep_score) %>%
    dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Dependency Score" = "dep_score")
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
    geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "gray", fill = "#02224C") +
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
    geom_point(aes(x = fct_reorder(cell_line, dep_score, .desc = FALSE), 
                   y = dep_score, 
                   text = paste0("Cell Line: ", cell_line)), 
               alpha = 0.2, color = "#02224C") +
    labs(x = "Cell Lines", y = "Dependency Score") +
    geom_hline(yintercept = mean_virtual_achilles) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    theme_cowplot() +
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
    dep_top_related <- make_top_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    dep_bottom_related <- make_bottom_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)

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
    dep_top_related <- make_top_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    dep_bottom_related <- make_bottom_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
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
  p1 <- make_celldeps(gene_symbol)
  p2 <- make_cellbins(gene_symbol)
  target_achilles_bottom <- make_achilles_table(gene_symbol) %>% slice(1:10)
  target_achilles_top <- make_achilles_table(gene_symbol) %>% dplyr::arrange(desc(`Dependency Score`)) %>%  slice(1:10)
  dep_top <- make_top_table(gene_symbol)
  flat_top_complete <- make_enrichment_table(master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(gene_symbol)
  flat_bottom_complete <- make_enrichment_table(master_negative, gene_symbol)
  graph_report <- make_graph_report(gene_symbol)
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



#UI------
ui <- fluidPage(
  tags$head(includeHTML("gtag.html"), 
            includeScript("returnClick.js")),
  navbarPage(
    title = "Data-Driven Hypothesis",
    tabPanel("Home",
             "Data-driven hypothesis is a resource developed by the", a(href = "http://www.hirscheylab.org", "Hirschey Lab"), "to functionally map human genes. A typical use case starts by querying a gene, identifying genes that share similar patterns or behaviors across several measures, in order to discover novel genes in established processes or new functions for well-studied genes.",
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
                                         fluidRow(htmlOutput(outputId = "plot_text"))
                                         )),
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
                                         fluidRow(checkboxGroupInput(inputId = "vars_dep_top", 
                                                                     "Select columns:",
                                                                     c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                                     selected = c("Z Score", "Co-publication Count"), 
                                                                     inline = TRUE)),
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
                                         fluidRow(checkboxGroupInput(inputId = "vars_dep_bottom", 
                                                                     "Select columns:",
                                                                     c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                                                     selected = c("Z Score", "Co-publication Count"), 
                                                                     inline = TRUE)),
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
             conditionalPanel(condition = 'input.go == 0',
                              "Enter a gene symbol to generate reports"),
             conditionalPanel(condition = 'input.go != 0',
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
                                               downloadButton(outputId = "report", label = "Download report"))))
  )
)

#SERVER-----
server <- function(input, output, session) {
  data <- eventReactive(input$go, {
    if_else(str_detect(input$gene_symbol, "orf"), input$gene_symbol, str_to_upper(input$gene_symbol))})
  
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
      make_top_table(data()) %>% dplyr::select("Gene", "Name", input$vars_dep_top), 
      options = list(pageLength = 25))
  })
  output$dep_bottom <- DT::renderDataTable({
    validate(
      need(data() %in% master_bottom_table$fav_gene, "No data found for this gene."))
    DT::datatable(
      make_bottom_table(data()) %>% dplyr::select("Gene", "Name", input$vars_dep_bottom),
    options = list(pageLength = 25))
  })
  output$cell_deps <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    withProgress(message = 'Wait for it...', value = 1, {
      ggplotly(make_celldeps(data()), tooltip = "text")
    })
  })
  output$cell_bins <- renderPlotly({
    validate(
      need(data() %in% colnames(achilles), "")) #""left blank
    ggplotly(make_cellbins(data()))
  })
  output$plot_text <- renderUI({
    validate(
      need(data() %in% colnames(achilles), "")) #""left blank
    HTML("<p></p><p><b>Left:</b> Cell Line Dependency Curve. Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essentail within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out. <b>Right:</b> Histogram of Binned Dependency Scores. Dependency scores across all cell lines for queried gene partitioned into 0.25 unit bins. Shape of the histogram curve reveals overall influence of queried gene on cellular fitness</p>")
    })
  output$target_achilles <- DT::renderDataTable({
    validate(
      need(data() %in% colnames(achilles), "No data found for this gene."))
    make_achilles_table(data())
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
    make_graph(data())
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

shinyApp(ui, server)
