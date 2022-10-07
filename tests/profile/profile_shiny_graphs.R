source(here::here("code", "fun_graphs.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)

testShinyModule("geneNetworkGraphServer", args = gene_data_args, {
  c(output$text_graph, output$graph)
})

testShinyModule("compoundNetworkGraphServer", args = compound_data_args, {
  c(output$text_compound_graph, output$compound_graph)
})
