downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    h2("Report Generator"),
    br(),
    conditionalPanel(condition = paste0("input['", ns("send_report_msg"), "'] == 0"),
                     "Enter your name and email address to receive a report", 
                     br(),
                     br(),
                     textInput(ns("first_name"), "First Name", ""), 
                     textInput(ns("last_name"), "Last Name", ""), 
                     textInput(ns("email_address"), "Email Address", ""), 
                     actionButton(inputId = ns("send_report_msg"), 
                                  label = "Generate")), 
    conditionalPanel(condition = paste0("input['", ns("send_report_msg"), "'] != 0"),
                     textOutput(ns("confirmation_message")))
  )
}

downloadReportPanelServer <- function(id, data, privateMode) {
  moduleServer(
    id,
    function(input, output, session) {
      #disable send_report_msg if email is empty, from https://deanattali.com/shinyjs/overview
      observe({
        if (is.null(input$email_address) || input$email_address == "") {
          shinyjs::disable("send_report_msg")
        } else {
          shinyjs::enable("send_report_msg")
        }
      })
      
      # user clicks generate report send message to sqs
      observeEvent(input$send_report_msg, {
        send_report_message(first_name = input$first_name,
                            last_name = input$last_name,
                            email_address = input$email_address,
                            input = data(),
                            private = privateMode)
      })
      
      # give message
      output$confirmation_message <- renderText({
        glue::glue("The {(data()$query)} report will be generated and emailed to {input$email_address}")
      })
      
      # alert
      # observeEvent(input$send_report_msg, {
      #   shinyjs::alert("While you're waiting, check out some other queries")
      # })
    })
}

