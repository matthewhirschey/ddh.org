cellPage <- function (id, type) {
  ns <- NS(id)

  summary_var <- ""
  if (type == "cell") {
    summary_var <- cellSummaryText(ns("summary"))
  }
  if (type == "lineage") {
    summary_var <- lineageSummaryText(ns("summary"))
  }
  if (type == "lineage_subtype") {
    summary_var <- lineageSubtypeSummaryText(ns("summary"))
  }
  if (type == "cell_list") {
    summary_var <- cellListSummaryText(ns("summary"))
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      tabPanel("SUMMARY",
               summary_var),
      formContent=querySearchInput(ns("search")))
  )
}

cellPageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {
      if(type == "cell") {
        data <- reactive({
          cell_line <- getQueryString()$cell_line
          list(
            type=type,
            id=cell_line,
            cell_line=cell_line
          )
        })
      }
      if(type == "lineage") {
        data <- reactive({
          lineage_str <- getQueryString()$lineage
          if (!is.null(lineage_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage == lineage_str)
            list(
              type=type,
              id=lineage_str,
              cell_line=expression_name_row$cell_line
            )
          }
        })
      }
      if(type == "lineage_subtype") {
        data <- reactive({
          lineage_subtype_str <- getQueryString()$lineage_subtype
          if (!is.null(lineage_subtype_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage_subtype == lineage_subtype_str)
            list(
              type=type,
              id=lineage_subtype_str,
              cell_line=expression_name_row$cell_line
            )
          }
        })
      }
      if(type == "cell_list") {
        data <- reactive({
          custom_cell_list <- getQueryString()$custom_cell_list
          cell_lines <- c(str_split(custom_cell_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_cell_list,
            cell_line=cell_lines
          )
        })
      }
      cellSummaryTextServer("summary", data)
      lineageSummaryTextServer("summary", data)
      lineageSubtypeSummaryTextServer("summary", data)
      cellListSummaryTextServer("summary", data)
      querySearchServer("search")
    }
  )
}
