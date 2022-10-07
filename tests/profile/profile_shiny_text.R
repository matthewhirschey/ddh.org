source(here::here("code", "fun_helper.R"))
source(here::here("code", "fun_text.R"))
source(here::here("code", "shiny_text.R"))

testShinyModule("geneTitleServer", args = gene_data_args, {
  c(output$gene_summary_title)
})

testShinyModule("pathwayTitleServer", args = pathway_data_args, {
  c(output$pathway_title)
})

testShinyModule("geneListTitleServer", args = gene_list_data_args, {
  c(output$custom_gene_list)
})

testShinyModule("cellTitleServer", args = cell_data_args, {
  c(output$cell_summary_title)
})

testShinyModule("lineageTitleServer", args = lineage_data_args, {
  c(output$lineage_title)
})

testShinyModule("lineageSubtypeTitleServer", args = lineage_subtype_data_args, {
  c(output$lineage_subtype_title)
})

testShinyModule("cellListTitleServer", args = cell_list_data_args, {
  c(output$custom_cell_list)
})

testShinyModule("tissueTitleServer", args = gene_data_args, {
  c(output$tissue_title)
})

testShinyModule("compoundTitleServer", args = compound_data_args, {
  c(output$compound_title)
})

testShinyModule("moaTitleServer", args = moa_data_args, {
  c(output$moa_title)
})

testShinyModule("compoundListTitleServer", args = compound_list_data_args, {
  c(output$custom_compound_list)
})

testShinyModule("geneTextServer", args = gene_data_args, {
  c(
    output$gene_summary_title,
    output$gene_summary_approved_symbol,
    output$gene_summary_approved_name,
    output$gene_summary_aka,
    output$gene_length,
    output$gene_summary_entrez_summary,
    output$ncbi_link
  )
})

testShinyModule("proteinTextServer", args = gene_data_args, {
  c(
    output$protein_summary_title,
    output$ec_link,
    output$protein_summary_uniprot_summary,
    output$uniprot_link,
    output$protein_summary_mass,
    output$protein_summary_seq
  )
})

testShinyModule("pathwayTextServer", args = pathway_data_args, {
  c(output$pathway_title,
    output$pathway_gene_symbols,
    output$pathway_def)
})

testShinyModule("geneListTextServer", args = gene_list_data_args, {
  c(output$custom_gene_list)
})

testShinyModule("cellSummaryTextServer", args = cell_data_args, {
  c(
    output$cell_summary_title,
    output$cell_summary_lineage,
    output$cell_summary_lineage_subtype,
    output$cell_description,
    output$cell_age,
    output$cell_sex
  )
})

testShinyModule("lineageSummaryTextServer", args = lineage_data_args, {
  c(output$lineage_summary_title,
    output$lineage_summary_cell_lines)
})

testShinyModule("lineageSubtypeSummaryTextServer", args = lineage_subtype_data_args, {
  c(
    output$lineage_subtype_summary_title,
    output$lineage_subtype_summary_cell_lines
  )
})

testShinyModule("cellListSummaryTextServer", args = cell_list_data_args, {
  c(output$title)
})

testShinyModule("compoundSummaryTextServer", args = compound_data_args, {
  c(
    output$compound_summary_title,
    output$compound_summary_description,
    output$compound_summary_moa,
    output$compound_summary_phase,
    output$compound_summary_diseasearea,
    output$compound_summary_indication,
    # output$compound_summary_targets, TODO put back in once "missing value where TRUE/FALSE needed" is fixed
    output$compound_summary_cid,
    output$cid_link
  )
})

testShinyModule("moaSummaryTextServer", args = moa_data_args, {
  c(output$moa_summary_title, output$moa_summary_compounds)
})

testShinyModule("compoundListSummaryTextServer", args = compound_list_data_args, {
  c(output$title)
})
