compoundPage <- function (id, type) {
  ns <- NS(id)

  summary_var <- ""
  if (type == "compound") {
    summary_var <- compoundSummaryText(ns("summary"))
  }
  if (type == "moa") {
    summary_var <- moaSummaryText(ns("summary"))
  }
  if (type == "compound_list") {
    summary_var <- compoundListSummaryText(ns("summary"))
  }

  tagList(
    head_tags,
    ddhNavbarPage(
      tabPanel("SUMMARY",
               summary_var),
      formContent=querySearchInput(ns("search")))
  )
}

compoundPageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {
      if(type == "compound") {
        data <- reactive({
          compound_name <- getQueryString()$compound
          list(
            type=type,
            id=compound_name,
            compound=compound_name
          )
        })
      }
      if(type == "moa") {
        data <- reactive({
          moa_str <- getQueryString()$moa
          if (!is.null(moa_str)) {
            prism_rows <- prism_names %>%
              filter(moa == moa_str)
            list(
              type=type,
              id=moa_str,
              compound=prism_rows$name
            )
          }
        })
      }
      if(type == "compound_list") {
        data <- reactive({
          custom_compound_list <- getQueryString()$custom_compound_list
          compounds_list <- c(str_split(custom_compound_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_compound_list,
            compound=compounds_list
          )
        })
      }
      compoundSummaryTextServer("summary", data)
      moaSummaryTextServer("summary", data)
      compoundListSummaryTextServer("summary", data)
      querySearchServer("search")
    }
  )
}
