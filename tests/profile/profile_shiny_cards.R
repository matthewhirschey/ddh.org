source(here::here("code", "fun_plots.R"), local = TRUE)
source(here::here("code", "fun_graphs.R"), local = TRUE)
source(here::here("code", "fun_tables.R"), local = TRUE)
source(here::here("code", "fun_structures.R"), local = TRUE)
source(here::here("code", "shiny_cards.R"), local = TRUE)

testShinyModule("cellDependenciesPlotDashServer", args = gene_data_args, {
   c(
      output$cell_depdash
   )
})

testShinyModule("cellDependenciesTableDashServer", args = gene_data_args, {
   c(
      output$deptabledash
   )
})

testShinyModule("cellDependenciesGraphDashServer", args = gene_data_args, {
   c(
      output$depgraphdash
   )
})

testShinyModule("ideogramPlotDashServer", args = gene_data_args, {
   c(output$ideogramdash)
})
testShinyModule("tissueAnatogramPlotDashServer", args = gene_data_args, {
   c(output$tissueanatogramdash)
})
testShinyModule("drugGenesCorTableDashServer", args = compound_data_args, {
   c(output$druggenescortabledash)
})

testShinyModule("geneDrugsCorTableDashServer", args = gene_data_args, {
   c(output$genedrugscortabledash)
})

testShinyModule("cellLineDependenciesPlotDashServer", args = cell_data_args, {
   c(output$cell_depdash)
})

testShinyModule("compoundStructureDashServer", args = compound_data_args, {
   c(output$compounddash)
})

testShinyModule("compoundDependenciesPlotDashServer", args = compound_data_args, {
   c(output$compound_depdash)
})

testShinyModule("compoundDependenciesTableDashServer", args = compound_data_args, {
   c(output$depcompoundtabledash)
})

testShinyModule("compoundDependenciesGraphDashServer", args = compound_data_args, {
   c(output$compoundgraphdash)
})
