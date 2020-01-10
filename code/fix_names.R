fix_names <- function(wrong_name) {
  var <- str_which(gene_summary$aka, paste0("(?<![:alnum:])", wrong_name, "(?![:alnum:]|\\-)")) #finds index
  df <- gene_summary[var,]
  right_name <- df$approved_symbol
  if (length(var) == 1) {
    return(right_name)
  } else {
    return(wrong_name)
  }
  #fixes 251, leaves 11
}

clean_colnames <- function(dataset) {
  for (name in names(dataset)) {
    if (name %in% gene_summary$approved_symbol == FALSE){
      names(dataset)[names(dataset) == name] <- fix_names(name)
    } else {
      name
    }
  }
}
