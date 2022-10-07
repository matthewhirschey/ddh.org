#DDH PARAMS-----
source(here::here("code", "app_params.R"), local = TRUE)

#LOAD LIBRARIES-----
source(here::here("code", "install_libraries.R"))

#ESTABLISH PRIVATE-----
source(here::here("code", "private.R"))

#DOWNLOAD/LOAD DATA-----
ddh::download_ddh_data(app_data_dir = app_data_dir,
                       test = testMode,
                       overwrite = FALSE)

ddh::load_ddh_data(app_data_dir = app_data_dir)

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_search.R"), local = TRUE)

#SHINY FUNCTIONS-----
source(here::here("code", "shiny_helper.R"), local = TRUE)
source(here::here("code", "shiny_tables.R"), local = TRUE)
source(here::here("code", "shiny_plots.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)
source(here::here("code", "shiny_reports.R"), local = TRUE)
source(here::here("code", "shiny_text.R"), local = TRUE)
source(here::here("code", "shiny_cards.R"), local = TRUE)
source(here::here("code", "shiny_download.R"), local = TRUE)

# HEAD----
head_tags <- tags$head(includeHTML("gtag.html"), includeCSS("styles.css"))

### universal elements
main_title <- HTML('<a href="." style="color:black;">DATA-DRIVEN HYPOTHESIS</a>')
window_title <- "Data-Driven Hypothesis | A Hirschey Lab Resource"

ddhNavbarPage <- function(..., formContent = NULL, id = NULL) {
  title_with_form <- tagList(
    main_title,
    tags$div(class="ddh-nav-form", formContent)
  )
  navbarPage(title = title_with_form, id = id, windowTitle = window_title, ...)
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
  metabolite="metabolite",
  compound_list="compound_list"
)

#HOME PAGE----
homePage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    HTML('<center><br><br><img src="ddh_logo.png", width = "338" ></center>'),
    tags$div(
      tags$br(),
      HTML("<center>Data-driven hypothesis is a resource developed by the <a href='http://www.hirscheylab.org' style='color:black;'>Hirschey Lab</a> for predicting functional relationships for thousands of genes across the human genome.</center>"), 
      tags$br(),
      tags$br()),
    HTML("<center>"),
    querySearchInput(ns("search")), 
    exampleSearchesLink(ns("examples")), 
    ", ", 
    HTML("<a href='/methods/start-here.html'>read the manual</a>"),
    ", or",
    getLuckyLink(ns("lucky")),
    HTML("</center>"),
    exampleSearchesPanel(ns("examples"))
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
    }
  )
}

##Search----
# module to input search term and navigate to the search screen
querySearchInput <- function(id) {
  ns <- NS(id)
  searchInput(
    inputId = ns("gene_or_pathway"),
    placeholder = "genes, pathways, or GO number", #change to genes, cells, or compounds
    btnSearch = icon("magnifying-glass")
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

##EXAMPLE SEARCHES----
# module to display a list of example searches

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
  source(here::here("code", "examples.R")) #pull out so methods can use it
  notZeroConditionalPanel(ns("example_click"), #toggle?
                          tagList(
                            tags$br(),
                            h3("Examples"),
                            HTML(examples), 
                            browsePathwaysLink(ns("pathways")),
                            browsePathwaysPanel(ns("pathways")) 
                          )
  )
}

exampleSearchesPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      browsePathwaysLinkServer("pathways")
      browsePathwaysPanelServer("pathways")
    }
  )
}

##LUCKY GENE----
# module to display a random interesting gene and navigate to the detail screen for that gene
getLuckyLink <- function(id) {
  ns <- NS(id)
  htmlOutput(ns("get_lucky"), inline = TRUE)
}

surprise <- function(surprise_vec) {
  gene_symbol <- sample(surprise_vec, 1)
  gene_symbol_url <- paste0("?show=gene&query=", gene_symbol)
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

# SEARCH PAGE ----
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
  func <- subtype_to_query_result_row[[row$subtype]]
  func(row)
}

searchPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      querySearchServer("search")
      output$search_title <- renderText({
        query_string <- getQueryString()
        paste0("Search results for '", query_string$query, "'")
      })
      output$genes_search_result <- renderUI({
        query_string <- getQueryString()
        query_results_table <- search_tables(gene_summary, pathways, expression_names, prism_names, hmdb_names, query_string$query)
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

# SEARCH RESULT ROWS ----

gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    h4(
      tags$strong("Gene:"),
      tags$a(title, href=paste0("?show=gene&query=", gene_summary_row["approved_symbol"]))
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
      tags$a(title, href=paste0("?show=pathway&query=", pathways_row$go))
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
    pull(approved_symbol) %>% 
    unique(.)
  has_known_gene_symbols <- !is_empty(known_gene_symbols)
  
  unknown_gene_symbols <- gene_summary_rows %>% 
    filter(known == FALSE) %>%
    pull(approved_symbol)
  has_unknown_gene_symbols <- !is_empty(unknown_gene_symbols)
  
  known_gene_symbols_tags <- NULL
  if (has_known_gene_symbols) {
    gene_query_param <- paste0("query=", paste(known_gene_symbols, collapse=","))
    href <- paste0("?show=gene_list&", gene_query_param)
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
      tags$a(title, href=paste0("?show=cell&query=", expression_names_row["cell_line"]))
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
      tags$a(row$title, href=paste0("?show=lineage&query=", row$key))
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
      tags$a(row$title, href=paste0("?show=lineage_subtype&query=", row$key))
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
    cell_list_param <- paste0("query=", paste(known_expression_names, collapse=","))
    href <- paste0("?show=cell_list&", cell_list_param)
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
      tags$a(prism_name_row$name, href=paste0("?show=compound&query=", prism_name_row$name))
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
      tags$a(row$title, href=paste0("?show=moa&query=", row$key))
    ),
    div(tags$strong("Compounds:"), compounds),
    hr()
  )
}

metabolite_query_result_row <- function(row) {
  hmdb_name_row <- row$data
  compounds <- paste0(hmdb_name_row$name, collapse=", ")
  list(
    h4(
      tags$strong("Metabolite:"),
      tags$a(hmdb_name_row$name, href=paste0("?show=metabolite&query=", row$name))
    ),
    div(tags$strong("Class:"), hmdb_name_row$class),
    div(tags$strong("CID:"), hmdb_name_row$cid),
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
    compound_list_param <- paste0("query=", paste(known_compound_names, collapse=","))
    href <- paste0("?show=compound_list&", compound_list_param)
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

# specifies how to render the results for a specific subtype
# functions that generate rows in fun_tables.R eg. gene_list_query_results_table()
subtype_to_query_result_row = list(
  gene=gene_query_result_row,
  pathway=pathway_query_result_row,
  gene_list=gene_list_query_result_row,
  cell=cell_query_result_row,
  lineage=lineage_query_result_row,
  lineage_subtype=lineage_subtype_query_result_row,
  cell_list=cell_list_query_result_row,
  compound=compound_query_result_row,
  moa=moa_query_result_row,
  metabolite=metabolite_query_result_row,
  compound_list=compound_list_query_result_row
)

# PAGE MODULES-----
source(here::here("code", "page_gene.R"), local = TRUE) ### GENE PAGE ----
source(here::here("code", "page_cell.R"), local = TRUE) ### CELL PAGE ----
source(here::here("code", "page_compound.R"), local = TRUE) ### COMPOUND PAGE ----

# Create output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageUI")
  )
)

pages <- list(
  home=homePage(page_names$home),
  search=searchPage(page_names$search),
  gene=genePage(page_names$gene, subtype = "gene"),
  pathway=genePage(page_names$pathway, subtype = "pathway"),
  gene_list=genePage(page_names$gene_list, subtype = "gene_list"),
  cell=cellPage(page_names$cell, subtype = "cell"),
  lineage=cellPage(page_names$lineage, subtype = "lineage"),
  lineage_subtype=cellPage(page_names$lineage_subtype, subtype = "lineage_subtype"),
  cell_list=cellPage(page_names$cell_list, subtype = "cell_list"),
  compound=compoundPage(page_names$compound, subtype = "compound"),
  moa=compoundPage(page_names$moa, subtype = "moa"),
  metabolite=compoundPage(page_names$metabolite, subtype ="metabolite"),
  compound_list=compoundPage(page_names$compound_list, subtype = "compound_list")
)

server <- shinyServer(function(input, output, session) {
  options(shiny.usecairo=TRUE) # ensure high quality images
  output$pageUI <- renderUI({
    query_string <- getQueryString()
    show_page <- query_string$show
    if (is.null(show_page)) {
      show_page <- page_names$home
    }
    pages[show_page]
  })
  homePageServer(page_names$home)
  searchPageServer(page_names$search)
  genePageServer(page_names$gene, subtype = "gene")
  genePageServer(page_names$pathway, subtype = "pathway")
  genePageServer(page_names$gene_list, subtype = "gene_list")
  cellPageServer(page_names$cell, subtype = "cell")
  cellPageServer(page_names$lineage, subtype = "lineage")
  cellPageServer(page_names$lineage_subtype, subtype = "lineage_subtype")
  cellPageServer(page_names$cell_list, subtype = "cell_list")
  compoundPageServer(page_names$compound, subtype = "compound")
  compoundPageServer(page_names$metabolite, subtype = "metabolite")
  compoundPageServer(page_names$moa, subtype = "moa")
  compoundPageServer(page_names$compound_list, subtype = "compound_list")
  # session$onSessionEnded(function() {
  #   delete_tmp_zip_directory(session)
  # })
})

shinyApp(ui, server)

