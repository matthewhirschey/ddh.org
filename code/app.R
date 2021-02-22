library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(visNetwork)
library(corrr)
library(here)
library(lubridate)
library(rmarkdown)
library(markdown)
library(tidygraph)
library(ggraph)
library(cowplot)
library(plotly)
library(DT)
library(future)
library(promises)
library(gganatogram)
library(ggdist)
library(showtext)
library(htmltools)
library(gt)

#LOAD DATA-----
#read params
source(here::here("code", "app_params.R"))

#read data from create_*.R
gene_summary <- readRDS(here::here(app_data_dir, paste0(release, "_gene_summary.Rds")))
pathways <- readRDS(here::here(app_data_dir, paste0(release, "_pathways.Rds")))

#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles.Rds")))
# achilles_cor <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles_cor.Rds")))
achilles_cor_nest <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles_cor_nest.Rds")))
expression <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression.Rds")))
expression_meta <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_meta.Rds")))
expression_names <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_names.Rds")))

#read drug names for search
prism_names <- readRDS(here::here(app_data_dir, paste0(release, "_prism_names.Rds"))) 

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here(app_data_dir, paste0(release, "_sd_threshold.Rds")))
achilles_lower <- readRDS(file = here::here(app_data_dir, paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here(app_data_dir, paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here(app_data_dir, paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here(app_data_dir, paste0(release, "_sd_virtual_achilles.Rds")))

expression_upper <- readRDS(file = here::here(app_data_dir, paste0(release, "_expression_upper.Rds")))
expression_lower <- readRDS(file = here::here(app_data_dir, paste0(release, "_expression_lower.Rds")))
mean_virtual_expression <- readRDS(file = here::here(app_data_dir, paste0(release, "_mean_virtual_expression.Rds")))
sd_virtual_expression <- readRDS(file = here::here(app_data_dir, paste0(release, "_sd_virtual_expression.Rds")))

#read data from generate_depmap_tables & pathways.R
master_bottom_table <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_bottom_table.Rds")))
master_top_table <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_top_table.Rds")))
master_positive <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_positive.Rds")))
master_negative <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_negative.Rds")))
surprise_genes <- readRDS(file=here::here(app_data_dir, paste0(release, "_surprise_genes.Rds")))
censor_genes <- readRDS(file=here::here(app_data_dir, paste0(release, "_censor_genes.Rds")))

#read data from generate_subcell_data.R
subcell <- readRDS(file=here::here(app_data_dir, paste0(release, "_subcell.Rds")))

#read data from generate_proteins_data.R
proteins <- readRDS(file=here::here(app_data_dir, paste0(release, "_proteins.Rds")))

#PRIVATE DATA CURRENTLY IN BETA
if(private == TRUE) {
  #load these data
} 

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_tables.R"), local = TRUE)
source(here::here("code", "fun_plots.R"), local = TRUE)
source(here::here("code", "fun_graphs.R"), local = TRUE)
source(here::here("code", "fun_reports.R"), local = TRUE)
source(here::here("code", "fun_text.R"), local = TRUE)

#SHINY FUNCTIONS-----
source(here::here("code", "shiny_tables.R"), local = TRUE)
source(here::here("code", "shiny_plots.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)
source(here::here("code", "shiny_reports.R"), local = TRUE)
source(here::here("code", "shiny_text.R"), local = TRUE)
source(here::here("code", "shiny_cards.R"), local = TRUE)

#MISC FUNCTIONS-----
font_add_google("Nunito Sans", "Nunito Sans")
font_add_google("Roboto Slab", "Roboto Slab")

render_report_in_background <- FALSE
if (supportsMulticore()) {
  plan(multicore)
  render_report_in_background <- TRUE
}

### HEAD
head_tags <- tags$head(includeHTML("gtag.html"), includeCSS("styles.css"))

### universal elements
main_title <- HTML('<a href="." style="color:black;">DATA-DRIVEN HYPOTHESIS</a>')
window_title <- "Data-Driven Hypothesis | A Hirschey Lab Resource"

ddhNavbarPage <- function(..., formContent = NULL, id = NULL) {
  navbarPageWithForm(title = main_title, id = id, windowTitle = window_title, formContent=formContent, ...)
}

navbarPageWithForm <- function (title, ..., id = NULL, windowTitle = title, formContent = NULL, formClass = "navbar-form navbar-right")
{
  # code below is based on shiny::navbarPage
  pageTitle <- title
  navbarClass <- "navbar navbar-default navbar-static-top"
  tabs <- list(...)
  tabset <- shiny:::buildTabset(tabs, "nav navbar-nav", NULL, id)
  containerDiv <- div(class = "container-fluid",
                      div(class = "navbar-header", span(class = "navbar-brand", pageTitle)),
                      div(class = formClass, formContent),
                      tabset$navList
                      )
  bootstrapPage(title = windowTitle,
                tags$nav(class = navbarClass, role = "navigation", containerDiv),
                div(class = "container-fluid", tabset$content))
}

### list of all pages rendered by this app
page_names <- list(
  home="home",
  search="search",
  gene="gene",
  pathway="pathway",
  gene_list="gene_list", 
  cell="cell", 
  lineage="lineage",
  cell_list="cell_list",
  compound="compound",
  moa="moa",
  compound_list="compound_list"
)

### HOME Modules ----

#SEARCH BOX----
# module to input search term and navigate to the search screen
querySearchInput <- function(id) {
  ns <- NS(id)
  searchInput(
    inputId = ns("gene_or_pathway"),
    placeholder = "genes, pathways, or GO number", #search
    btnSearch = icon("search")
  )
}

querySearchServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$gene_or_pathway, {
        if (input$gene_or_pathway != '') {
          updateQueryString(paste0("?show=search&query=", input$gene_or_pathway), mode="push")
        }
      })
    }
  )
}

#LUCKY GENE----
# module to display a random interesting gene and navigate to the detail screen for that gene
getLuckyLink <- function(id) {
  ns <- NS(id)
  htmlOutput(ns("get_lucky"), inline = TRUE)
}

surprise <- function(surprise_vec) {
  gene_symbol <- sample(surprise_vec, 1)
  gene_symbol_url <- paste0("?show=gene&query_type=gene&symbol=", gene_symbol)
  return(gene_symbol_url)
}

getLuckyServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$get_lucky <- renderUI({
        tags$a(href = surprise(surprise_genes), "get lucky")
      })
    }
  )
}

#EXAMPLE SEARCHES----
# module to display a random interesting gene and navigate to the detail screen for that gene

exampleSearchesLink <- function(id) {
  ns <- NS(id)
  actionLink(inputId = ns("example_click"), "See example searches")
}

exampleSearchesLinkServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$example_click, {}) # event to store the 'click'
    }
  )
}

exampleSearchesPanel <- function(id) {
  ns <- NS(id)
  notZeroConditionPanel(ns("example_click"), 
                        tags$br(),
                        h4("Examples"),
                        HTML('<h5>Search for</h5>
                        <ul>
                        <li>A single gene, such as <a href="?show=gene&query_type=gene&symbol=TP53">TP53</a> or <a href="?show=gene&query_type=gene&symbol=BRCA1">BRCA1</a></li>
                        <li>A pathway name, such as <a href="?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="?show=pathway&query_type=pathway&go=0006695">Cholesterol Biosynthetic Process</a></li>
                        <li>The Gene Ontology biological process identifier, such as <a href="?show=search&query=1901989">1901989</a>, which will find <a href="?show=pathway&query_type=pathway&go=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                        <li>A custom list of genes (separated by commas), such as <a href="?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="?show=gene_list&query_type=custom_gene_list&custom_gene_list=BRCA1,BRCA2">a custom gene list</a></li>
                       </ul>')
  )
}

exampleSearchesPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
    }
  )
}

#HOME(landing) PAGE----
homePage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    HTML('<center><br><br><img src="ddh_logo.png", width = "338" ></center>'),
    tags$div(
      tags$br(),
      HTML('<center>Data-driven hypothesis is a resource developed by the <a href="http://www.hirscheylab.org" style="color:black;">Hirschey Lab</a> for predicting functional relationships for thousands of genes across the human genome.</center>'), 
      tags$br(),
      tags$br()),
    HTML("<center>"),
    querySearchInput(ns("search")), 
    exampleSearchesLink(ns("examples")), 
    ", ", 
    browsePathwaysLink(ns("pathways")),
    ", or",
    getLuckyLink(ns("lucky")),
    HTML("</center>"),
    exampleSearchesPanel(ns("examples")),
    browsePathwaysPanel(ns("pathways")) 
  )
}

homePageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      querySearchServer("search")
      getLuckyServer("lucky")
      exampleSearchesLinkServer("examples")
      exampleSearchesPanelServer("examples")
      browsePathwaysLinkServer("pathways")
      browsePathwaysPanelServer("pathways")
    }
  )
}

### SEARCH PAGE ----
searchPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage(formContent=querySearchInput(ns("search"))),
    h3(textOutput("search_title")),
    div(div(h3("Results", class="panel-title"), class="panel-heading"),
        div(uiOutput(ns("genes_search_result")), class="panel-body"),
        class="bg-info panel panel-default"
    )
  )
}

query_result_row <- function(row) {
  func <- query_type_to_query_result_row[[row$query_type]]
  func(row)
}

searchPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      querySearchServer("search")
      output$search_title <- renderText({
        query <- getQueryString()
        paste0("Search results for '", query$query, "'")
      })
      output$genes_search_result <- renderUI({
        query <- getQueryString()
        query_results_table <- search_tables(gene_summary, pathways, expression_names, prism_names, query$query)
        if (nrow(query_results_table) > 0) {
          apply(query_results_table, 1, query_result_row)
        }
        else {
          "No results found."
        }
      })      
    }
  )
}

### SEARCH RESULT ROWS ----

gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    h4(
      tags$strong("Gene:"),
      tags$a(title, href=paste0("?show=gene&query_type=gene&symbol=", gene_summary_row["approved_symbol"]))
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
      tags$a(title, href=paste0("?show=pathway&query_type=pathway&go=", pathways_row$go))
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
    href <- paste0("?show=gene_list&query_type=custom_gene_list&", gene_query_param)
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
      tags$strong("Custom Gene List"),
      tags$span(title)
    ),
    known_gene_symbols_tags,
    unknown_gene_symbols_tags,
    hr()
  )
}

cell_query_result_row <- function(row) {
  expression_names_row <- row$data
  title <- paste0(expression_names_row["cell_line"])
  list(
    h4(
      tags$strong("Cell:"),
      tags$a(title, href=paste0("?show=cell&query_type=cell&cell_line=", expression_names_row["cell_line"]))
    ),
    div(tags$strong("Lineage:"), expression_names_row["lineage"]),
    div(tags$strong("Sublineage:"), expression_names_row["lineage_subtype"]),
    hr()
  )
}

lineage_query_result_row <- function(row) {
  expression_names_row <- row$data$data[[1]]
  cell_lines <- paste0(expression_names_row$cell_line, collapse=", ")
  list(
    h4(
      tags$strong("Lineage:"),
      tags$a(row$title, href=paste0("?show=lineage&query_type=lineage&lineage=", row$key))
    ),
    div(tags$strong("Cells:"), cell_lines),
    hr()
  )
}

lineage_subtype_query_result_row <- function(row) {
  expression_names_row <- row$data$data[[1]]
  cell_lines <- paste0(expression_names_row$cell_line, collapse=", ")
  list(
    h4(
      tags$strong("Sublineage:"),
      tags$a(row$title, href=paste0("?show=lineage_subtype&query_type=lineage_subtype&lineage_subtype=", row$key))
    ),
    div(tags$strong("Cells:"), cell_lines),
    hr()
  )
}

cell_list_query_result_row <- function(row) {
  expression_names_rows <- row$data
  title <- row$key

  known_expression_names <- expression_names_rows %>%
    filter(known == TRUE) %>%
    pull(cell_line)
  has_known_expression_names <- !is_empty(known_expression_names)

  unknown_expression_names <- expression_names_rows %>%
    filter(known == FALSE) %>%
    pull(cell_line)
  has_unknown_expression_names <- !is_empty(unknown_expression_names)

  known_cell_line_tags <- NULL
  if (has_known_expression_names) {
    cell_list_param <- paste0("custom_cell_list=", paste(known_expression_names, collapse=","))
    href <- paste0("?show=cell_list&query_type=custom_cell_list&", cell_list_param)
    known_cell_line_tags <- list(
      tags$h6("Known Cell Lines"),
      tags$a(paste(known_expression_names, collapse=", "), href=href)
    )
  }

  unknown_cell_line_tags <- NULL
  if (has_unknown_expression_names) {
    unknown_cell_line_tags <- list(
      tags$h6("Unknown Cell Lines"),
      tags$div(paste(unknown_expression_names, collapse=", "))
    )
  }

  list(
    h4(
      tags$strong("Custom Cell Line List"),
      tags$span(title)
    ),
    known_cell_line_tags,
    unknown_cell_line_tags,
    hr()
  )
}

compound_query_result_row <- function(row) {
  prism_name_row <- row$data
  list(
    h4(
      tags$strong("Compound:"),
      tags$a(prism_name_row$name, href=paste0("?show=compound&query_type=compound&compound=", prism_name_row$name))
    ),
    div(tags$strong("Mechanism of Action:"), prism_name_row$moa),
    div(tags$strong("CID:"), prism_name_row$cid),
    hr()
  )
}

moa_query_result_row <- function(row) {
  prism_name_row <- row$data$data[[1]]
  compounds <- paste0(prism_name_row$name, collapse=", ")
  list(
    h4(
      tags$strong("Compound Mechanism of Action:"),
      tags$a(row$title, href=paste0("?show=moa&query_type=moa&moa=", row$key))
    ),
    div(tags$strong("Compounds:"), compounds),
    hr()
  )
}

compound_list_query_result_row <- function(row) {
  prism_names_rows <- row$data
  title <- row$key

  known_compound_names <- prism_names_rows %>%
    filter(known == TRUE) %>%
    pull(name)
  has_known_compound_names <- !is_empty(known_compound_names)

  unknown_compound_names <- prism_names_rows %>%
    filter(known == FALSE) %>%
    pull(name)
  has_unknown_compound_names <- !is_empty(unknown_compound_names)

  known_compound_tags <- NULL
  if (has_known_compound_names) {
    compound_list_param <- paste0("custom_compound_list=", paste(known_compound_names, collapse=","))
    href <- paste0("?show=compound_list&query_type=compound_list&", compound_list_param)
    known_compound_tags <- list(
      tags$h6("Known Compounds"),
      tags$a(paste(known_compound_names, collapse=", "), href=href)
    )
  }

  unknown_compound_tags <- NULL
  if (has_unknown_compound_names) {
    unknown_compound_tags <- list(
      tags$h6("Unknown Compounds"),
      tags$div(paste(unknown_compound_names, collapse=", "))
    )
  }

  list(
    h4(
      tags$strong("Custom Compound List"),
      tags$span(title)
    ),
    known_compound_tags,
    unknown_compound_tags,
    hr()
  )
}

# specifies how to render the results for a specific query_type
# functions that generate rows in fun_tables.R eg. gene_list_query_results_table()
query_type_to_query_result_row = list(
  gene=gene_query_result_row,
  pathway=pathway_query_result_row,
  gene_list=gene_list_query_result_row,
  cell=cell_query_result_row,
  lineage=lineage_query_result_row,
  lineage_subtype=lineage_subtype_query_result_row,
  cell_list=cell_list_query_result_row,
  compound=compound_query_result_row,
  moa=moa_query_result_row,
  compound_list=compound_list_query_result_row
)

# PAGE MODULES-----
if(private == TRUE) {
  source(here::here("code", "page_gene_private.R"), local = TRUE) ### PRIVATE GENE PAGE ----
  source(here::here("code", "page_cell_private.R"), local = TRUE) ### PRIVATE CELL PAGE ----
  source(here::here("code", "page_compound_private.R"), local = TRUE) ### PRIVATE COMPOUND PAGE ----

} else {
  source(here::here("code", "page_gene.R"), local = TRUE) ### GENE PAGE ----
  source(here::here("code", "page_cell.R"), local = TRUE) ### CELL PAGE ----
  source(here::here("code", "page_compound.R"), local = TRUE) ### COMPOUND PAGE ----
} 

# Create output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageUI")
  )
)

pages <- list(
  home=homePage(page_names$home),
  search=searchPage(page_names$search),
  gene=genePage(page_names$gene, type = "gene"), #put var here
  pathway=genePage(page_names$pathway, type = "pathway"),
  gene_list=genePage(page_names$gene_list, type = "gene_list"),
  cell=cellPage(page_names$cell, type = "cell"),
  lineage=cellPage(page_names$lineage, type = "lineage"),
  lineage_subtype=cellPage(page_names$lineage_subtype, type = "lineage_subtype"),
  cell_list=cellPage(page_names$cell_list, type = "cell_list"),
  compound=compoundPage(page_names$compound, type = "compound"),
  moa=compoundPage(page_names$moa, type = "moa"),
  compound_list=compoundPage(page_names$compound_list, type = "compound_list")
)

server <- shinyServer(function(input, output, session) {
  output$pageUI <- renderUI({
    query <- getQueryString()
    show_page <- query$show
    if (is.null(show_page)) {
      show_page <- page_names$home
    }
    pages[show_page]
  })
  homePageServer(page_names$home)
  searchPageServer(page_names$search)
  genePageServer(page_names$gene, type = "gene")
  genePageServer(page_names$pathway, type = "pathway")
  genePageServer(page_names$gene_list, type = "gene_list")
  cellPageServer(page_names$cell, type = "cell")
  cellPageServer(page_names$lineage, type = "lineage")
  cellPageServer(page_names$lineage_subtype, type = "lineage_subtype")
  cellPageServer(page_names$cell_list, type = "cell_list")
  compoundPageServer(page_names$compound, type = "compound")
  compoundPageServer(page_names$moa, type = "moa")
  compoundPageServer(page_names$compound_list, type = "compound_list")
  session$onSessionEnded(function() {
    delete_tmp_zip_directory(session)
  })
})

shinyApp(ui, server)

