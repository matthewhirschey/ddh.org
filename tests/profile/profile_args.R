gene_data_args <- list(data <- function() {
  list(
    type="gene",
    subtype="gene",
    id="ROCK2",
    gene_symbols=c("ROCK2")
  )
})

pathway_data_args <- list(data <- function() {
  list(
    type="gene",
    subtype="pathway",
    id="1902965",
    gene_symbols=c("ROCK2")
  )
})

gene_list_data_args <- list(data <- function() {
  list(
    type="gene",
    subtype="gene_list",
    id="BRCA1,ROCK2",
    gene_symbols=c("BRCA1", "ROCK2")
  )
})

cell_data_args <- list(data <- function () {
  list(
    type="cell",
    subtype="cell",
    id="JHH6",
    cell_line=c("JHH6")
  )
})

lineage_data_args <- list(data=function () {
  list(
    type="cell",
    subtype="lineage",
    id="Liver",
    cell_line=c("JHH6")
  )
})

lineage_subtype_data_args <- list(data=function() {
  list(
    type="cell",
    subtype="lineage_subtype",
    id="Hepatocellular Carcinoma",
    cell_line=c("JHH6", "SNU398")
  )
})

cell_list_data_args <- list(data=function() {
  list(
    type="cell",
    subtype="cell_list",
    id="JHH6",
    cell_line=c("JHH6")
  )
})

compound_data_args <- list(data=function() {
  list(
    type="compound",
    subtype="compound",
    id="zaltoprofen",
    compound=c("zaltoprofen")
  )
})

moa_data_args <- list(data=function() {
  list(
    type="compound",
    subtype="moa",
    id="cyclooxygenase inhibitor",
    compound=c("floctafenine")
  )
})

compound_list_data_args <- list(data=function() {
  list(
    type="compound",
    subtype="compound_list",
    id="floctafenine, tropesin",
    compound=c("floctafenine", "tropesin")
  )
})
