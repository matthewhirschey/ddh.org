##SHINY DOWNLOADS
downloadTab <- function (id) {
    ns <- NS(id)
    
    tabPanel("DOWNLOADS",
             column(3, 
                    wellPanel(
                      fluidRow(
                        h4("ABOUT")
                      ),
                      fluidRow(#make sure the methods file is in www dir
                        HTML('<br><center><a href=\"methods/index.html\" target=\"_blank\"><img src=\"noun_Document.png\", width=\"50\", height=\"50\"></a><br>'), 
                        "METHODS", 
                        HTML('</center>'),
                      ),
                      fluidRow(
                        HTML('<br><center><a href=\"https://www.biorxiv.org/content/10.1101/2020.07.17.208751v1\" target=\"_blank\"><img src=\"noun_Document.png\", width=\"50\", height=\"50\"></a><br>'), 
                        "PREPRINT", 
                        HTML('</center>'),
                      ),
                      fluidRow(
                        HTML('<br><center><a href=\"https://github.com/matthewhirschey/ddh\" target=\"_blank\"><img src=\"GitHub-Mark-64px.png\", width=\"50\", height=\"50\"></a><br>'), 
                        "CODE", 
                        HTML('</center>'),
                      ),
                      fluidRow(
                        HTML('<br><center><a href=\"methods.html\" target=\"_blank\"><img src=\"noun_Document.png\", width=\"50\", height=\"50\"></a><br>'), 
                        "FAQ", 
                        HTML('</center>')
                      ), 
                      fluidRow(
                        HTML('<br><center><a href=\"https://www.twitter.com/ddhypothesis\" target=\"_blank\"><img src=\"2021 Twitter logo - black.png\", width=\"50\", height=\"50\"></a><br>'), 
                        "TWITTER", 
                        HTML('</center>'),
                      )
                    )),
             column(9, 
                    downloadReportPanel(ns("download"))
                    )
    )
}

downloadTabServer <- function (id, data, privateMode) {
  moduleServer(
    id,
    function(input, output, session) {
      downloadReportPanelServer("download", data, privateMode)
      
    }
  )
}
