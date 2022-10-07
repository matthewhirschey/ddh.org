# Installs all requirements for generating data and running the app

cran_pkgs <- c(
  "biomaRt",
  "Cairo",
  "caret",
  "colorspace",
  "corrr",
  "cowplot",
  "dbscan",
  "dplyr",
  "doParallel",
  "DT",
  "enrichR",
  "furrr",
  "future",
  "ggtext",
  "glue",
  "grDevices",
  "gt",
  "here",
  "hexbin",
  "httr",
  "htmltools",
  "igraph",
  "janitor",
  "jsonlite",
  "knitr",
  "lexicon",
  "lubridate",
  "magick",
  "magrittr",
  "moderndive",
  "r3dmol",
  "optparse",
  "pander",
  "patchwork",
  "paws",
  "plotly",
  "promises",
  "purrr",
  "quarto",
  "rayvertex",
  "RCurl",
  "readxl",
  "remotes",
  "rentrez",
  "reticulate",
  "RISmed",
  "rlang",
  "rmarkdown",
  "rtweet",
  "rvest",
  "scales",
  "shinycssloaders",
  "shinyjs",
  "shinyWidgets",
  "showtext",
  "stopwords",
  "stringr",
  "tibble",
  "tidygraph",
  "tidyr",
  "tidyselect",
  "tidytext",
  "uwot",
  "viridis",
  "visNetwork",
  "vroom",
  "webchem",
  "widyr",
  "XML"
  )

bioc_pkgs <- c(
  "AnnotationDbi",
  "biomaRt",
  "clusterProfiler",
  "impute",
  "org.Hs.eg.db"
  )

github_pkgs <- c(
  "jespermaag/gganatogram",
  "matthewhirschey/ddh",
  "matthewhirschey/tidybiology",
  "tylermorganwall/raymolecule"
  )

# Install CRAN packages
installifnot <- function(pckgName){
    if (!(require(pckgName, character.only = TRUE))) {
        install.packages(pckgName, dep = TRUE)
        library(pckgName, character.only = TRUE)
    }
}

lapply(cran_pkgs, installifnot)

# Install Bioconductor packages
installBiocifnot <- function(pckgName){
    if (!(require(pckgName, character.only = TRUE))) {
        BiocManager::install(pckgName)
        library(pckgName, character.only = TRUE)
    }
}

lapply(bioc_pkgs, installBiocifnot)

# Install GitHub packages
installGitHubifnot <- function(pckgName, pckgShortName){
    if (!(require(pckgShortName, character.only = TRUE))) {
        remotes::install_github(pckgName)
        library(pckgShortName, character.only = TRUE)
    }
}

github_pkgs_name <- sub('.+/(.+)', '\\1', github_pkgs)
mapply(installGitHubifnot, github_pkgs, github_pkgs_name)

