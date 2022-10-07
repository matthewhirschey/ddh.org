private <- function(x, publicValue=NULL) {
  if (privateMode == TRUE) { x } else { publicValue }
}

private_message_html <- "This analysis is currently in Beta, and is only available to premium subscribers. If you want to upgrade to support this project and access these data, subscribe at ', a(href = 'http://www.datadrivenhypothesis.com', 'datadrivenhypothesis.com')"
private_message_md <- "This analysis is currently in Beta, and is only available to premium subscribers. If you want to upgrade to support this project and access these data, subscribe at [datadrivenhypothesis.com](http://www.datadrivenhypothesis.com)"

private_msg <- function() {
  if (privateMode == FALSE) {
    htmltools::tagList(private_message_html) 
  }
}
