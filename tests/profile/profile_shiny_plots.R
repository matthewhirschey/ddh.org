source(here::here("code", "fun_helper.R"))
source(here::here("code", "fun_text.R"))
source(here::here("code", "fun_tables.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "shiny_plots.R"))

testShinyModule("cellDependenciesPlotServer", args = gene_data_args, {
  c(
      output$text_cell_dep_plot,
      output$essential_num,
      output$cell_deps
  )
})

testShinyModule("cellBinsPlotServer", args = gene_data_args, {
  c(
    output$cell_bins
  )
})

testShinyModule("cellDepsLinPlotServer", args = gene_data_args, {
  c(
      output$cell_deps_lin
  )
})

testShinyModule("cellDepsSubLinPlotServer", args = gene_data_args, {
  c(
      output$cell_deps_sublin
  )
})

testShinyModule("geneCorrelationPlotServer", args = gene_data_args, {
  c(
      output$text_genecorr_plot,
      output$gene_correlations
  )
})

testShinyModule("compoundDependenciesPlotServer", args = compound_data_args, {
  c(
      output$text_compound_dep_plot,
      output$essential_num_compound,
      output$cell_deps
  )
})

testShinyModule("compoundBinsPlotServer", args = compound_data_args, {
  c(
      output$compound_bins
  )
})

testShinyModule("compoundLinPlotServer", args = compound_data_args, {
  c(
      output$compound_lin
  )
})

testShinyModule("compoundSubLinPlotServer", args = compound_data_args, {
  c(
      output$compound_sublin
  )
})

testShinyModule("compoundCorrelationPlotServer", args = compound_data_args, {
  c(
      output$text_compoundcorr_plot,
      output$compound_correlations
  )
})

testShinyModule("cellLineDependenciesPlotServer", args = cell_data_args, {
  c(
      output$text_cell_dep_cell_plot,
      output$essential_num_cell,
      output$cell_deps_cell
  )
})

testShinyModule("cellLineBinsPlotServer", args = cell_data_args, {
  c(
      output$cell_line_bins
  )
})

testShinyModule("cellCorrelationPlotServer", args = cell_data_args, {
  c(
      output$text_cell_corr_plot,
      output$cell_correlations
  )
})

testShinyModule("cellAnatogramPlotServer", args = gene_data_args, {
  c(
      output$text_subcell_exp_plot,
      output$cellanatogram
  )
})

testShinyModule("cellAnatogramFacetPlotServer", args = gene_data_args, {
  c(
      output$cellanatogramfacet
  )
})

testShinyModule("maleAnatogramPlotServer", args = gene_data_args, {
  c(
      output$maleanatogram
  )
})

testShinyModule("femaleAnatogramPlotServer", args = gene_data_args, {
  c(
      output$femaleanatogram
  )
})

testShinyModule("tissuePlotServer", args = gene_data_args, {
  c(
      output$text_tissue_exp_plot,
      output$tissueplot
  )
})

testShinyModule("cellGeneExpressionPlotServer", args = gene_data_args, {
  c(
      output$text_cell_gene_plot,
      output$cell_gene_plot
  )
})

testShinyModule("cellProteinExpressionPlotServer", args = gene_data_args, {
  c(
      output$text_cell_protein_plot,
      output$cell_protein_plot
  )
})

testShinyModule("cellGeneProteinPlotServer", args = gene_data_args, {
  c(
      output$text_cell_gene_protein_plot,
      output$cell_gene_protein_plot
  )
})

testShinyModule("cellLineGeneExpressionPlotServer", args = cell_data_args, {
  c(
      output$text_cellLine_gene_plot,
      output$cellLine_gene_plot
  )
})

testShinyModule("cellLineProteinExpressionPlotServer", args = cell_data_args, {
  c(
      output$text_cellLine_protein_plot,
      output$cellLine_protein_plot
  )
})

testShinyModule("cellLineGeneProteinPlotServer", args = cell_data_args, {
  c(
      output$text_cellLine_gene_protein_plot,
      output$cellLine_gene_protein_plot
  )
})

testShinyModule("cellexpdepPlotServer", args = gene_data_args, {
  c(
      output$text_cellexpdep_plot,
      output$cellexpdep
  )
})

testShinyModule("ideogramPlotServer", args = gene_data_args, {
  c(
      output$text_ideogram_plot,
      output$ideogram
  )
})

testShinyModule("pubmedGenePlotServer", args = gene_data_args, {
  c(
      output$title_pubmed_plot,
      output$text_pubmed_plot,
      output$pubmed_plot
  )
})

testShinyModule("pubmedCompoundPlotServer", args = compound_data_args, {
  c(
      output$title_pubmed_compound_plot,
      output$text_pubmed_compound_plot,
      output$pubmed_compound_plot
  )
})
