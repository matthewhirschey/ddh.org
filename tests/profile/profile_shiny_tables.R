source(here::here("code", "fun_helper.R"), local = TRUE)
source(here::here("code", "fun_tables.R"), local = TRUE)
source(here::here("code", "shiny_tables.R"), local = TRUE)


testShinyModule("pathwayListServer", args = pathway_data_args, {
  c(
      output$text_pathway_list,
      output$pathway_list
  )
})

testShinyModule("pathwayGeneListServer", args = pathway_data_args, {
  c(
      output$text_pathway_gene_list,
      output$pathway_gene_list
  )
})

testShinyModule("cellDependenciesTableServer", args = gene_data_args, {
  session$setInputs(dep_filter_click = FALSE)
  c(
      output$text_cell_dep_table,
      output$target_achilles
  )
})

testShinyModule("compoundDependenciesTableServer", args = compound_data_args, {
  session$setInputs(compound_dep_filter_click = FALSE)
  c(
      output$text_compound_dep_table,
      output$compound_dep_table
  )
})

testShinyModule("cellLineDependenciesTableServer", args = cell_data_args, {
  session$setInputs(dep_cell_line_filter_click = FALSE)
  c(
      output$text_cell_line_dep_table,
      output$cell_line_achilles
  )
})

testShinyModule("similarGenesTableServer", args = gene_data_args, {
  c(
      output$text_dep_top,
      output$dep_top,
      output$text_pos_enrich,
      output$pos_enrich
  )
})

testShinyModule("similarCellsTableServer", args = cell_data_args, {
  c(
      output$text_cells_dep_top,
      output$cells_dep_top
  )
})

testShinyModule("similarCompoundsTableServer", args = compound_data_args, {
  c(
      output$text_compound_dep_top,
      output$compound_dep_top
  )
})

testShinyModule("dissimilarGenesTableServer", args = gene_data_args, {
  c(
      output$text_dep_bottom,
      output$dep_bottom,
      output$text_neg_enrich,
      output$neg_enrich
  )
})

# TODO add back in once dissimilarCompoundsTableServer is fixed
#testShinyModule("dissimilarCompoundsTableServer", args = compound_data_args, {
#  c(
#      output$text_compound_dep_bottom,
#      output$compound_dep_bottom
#  )
#})

testShinyModule("cellAnatogramTableServer", args = gene_data_args, {
  c(
      output$cellanatogram_table
  )
})

testShinyModule("tissueAnatogramTableServer", args = gene_data_args, {
  session$setInputs(tissue_filter_click = FALSE)
  c(
      output$tissueanatogram_table
  )
})

testShinyModule("cellGeneExpressionTableServer", args = gene_data_args, {
  c(
      output$text_cell_gene_table,
      output$cell_gene_table
  )
})

testShinyModule("cellProteinExpressionTableServer", args = gene_data_args, {
  c(
      output$text_cell_protein_table,
      output$cell_protein_table
  )
})

testShinyModule("cellLineGeneExpressionTableServer", args = cell_data_args, {
  c(
      output$text_cellLine_gene_table,
      output$cellLine_gene_table
  )
})

testShinyModule("cellLineProteinExpressionTableServer", args = cell_data_args, {
  c(
      output$text_cellLine_protein_table,
      output$cellLine_protein_table
  )
})

testShinyModule("pubmedGeneTableServer", args = gene_data_args, {
  c(
      output$pubmed_table
  )
})

testShinyModule("pubmedCompoundTableServer", args = compound_data_args, {
  c(
      output$pubmed_compound_table
  )
})

testShinyModule("geneDrugsTableServer", args = gene_data_args, {
  c(
      output$title_gene_drugs_table,
      output$text_gene_drugs_table,
      output$gene_drugs_table
  )
})

testShinyModule("drugGenesTableServer", args = compound_data_args, {
  c(
      output$title_drug_genes_table,
      output$text_drug_genes_table,
      output$drug_genes_table
  )
})

testShinyModule("geneDrugsCorTableServer", args = gene_data_args, {
  c(
      output$title_gene_drugs_cor_table,
      output$text_gene_drugs_cor_table,
      output$gene_drugs_cor_table
  )
})

testShinyModule("drugGenesCorTableServer", args = compound_data_args, {
  c(
      output$title_drug_genes_cor_table,
      output$text_drug_genes_cor_table,
      output$drug_genes_cor_table
  )
})