downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Report Generator"),
    textOutput(ns("help_message")),
    br(),
    conditionalPanel(condition = paste0("input['", ns("generate_report"), "'] == 0"),
      actionButton(ns("generate_report"), "Generate report")
    ),
    conditionalPanel(condition = paste0("output['", ns("report_zip_path"), "'] != ''"),
                     downloadButton(outputId = ns("report"), label = "Download report",
                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
    )
  )
}

downloadReportPanelServer <- function(id, type, data, private) {
  moduleServer(
    id,
    function(input, output, session) {
      # reactive variable to hold path to the download zip
      report_zip_path <- reactiveVal("")
      # make zip path a reactive output so it can be used in a conditionalPanel
      output$report_zip_path <- renderText({
        report_zip_path()
      })
      # populate the report_zip_path variable even though it isn't displayed
      outputOptions(output, "report_zip_path", suspendWhenHidden = FALSE)

      # update help message based on status of report generation
      output$help_message <- renderText({
        if (input$generate_report == 0) {
          "To generate a report, click on the button below."
        } else {
          if (report_zip_path() != "") {
            "Report complete. Click on the button below to download."
          }
        }
      })

      # user clicks generate report save zip into temp_zip_dir
      observeEvent(input$generate_report, {
        if(type == "gene"){
          filename <- paste0(data()$id, "_ddh.zip")
        } else if (type == "pathway") {
          filename <- paste0("go_", data()$id, "_ddh.zip")
        } else {
          filename <- paste0("custom_", paste(data()$gene_symbols, collapse="_"), "_ddh.zip")
        }
        zip_path <- file.path(get_or_create_tmp_zip_directory(session), filename)

        data_values <- data() # reactive data must be read outside of a future
        progress_bar <- Progress$new()
        progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
        if (render_report_in_background) {
          result <- future({
            render_report_to_file(data_values=data_values, file=zip_path, private)
          })
          finally(result, function(){
            report_zip_path(zip_path)
            progress_bar$close()
          })
        } else {
          render_report_to_file(data_values=data_values, file=zip_path, private)
          report_zip_path(zip_path)
          progress_bar$close()
        }
      })

      output$report <- downloadHandler(
        # create pdf report
        filename = function() {
          basename(report_zip_path())
        },
        content = function(file) {
          source_zip_path <- report_zip_path()
          message(file)
          message(source_zip_path)
          file.copy(source_zip_path, file)
        }
      )
    }
  )
} 

get_or_create_tmp_zip_directory <- function(session) {
  if (is.null(session$userData$temp_zip_dir)) {
    session$userData$temp_zip_dir <- tempfile(pattern="tmpzip", tmpdir=here::here("report"))
    dir.create(session$userData$temp_zip_dir)
  }
  session$userData$temp_zip_dir
}

delete_tmp_zip_directory <- function(session) {
  temp_zip_dir <- session$userData$temp_zip_dir
  if (!is.null(temp_zip_dir)) {
    unlink(temp_zip_dir, recursive = TRUE)
    session$userData$temp_zip_dir <- NULL
  }
}
