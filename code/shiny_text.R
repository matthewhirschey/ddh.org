# HELPER -----
#DUMMY
nameText <- function (id) {
  ns <- NS(id)
  "test"
}

nameTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
    }
  )
}

##TEMPLATE
# nameText <- function (id) {
#     ns <- NS(id)
# }
# 
# nameTextServer <- function (id, data) {
#   moduleServer(
#     id,
#     function(input, output, session) {
#     }
#   )
# }

#text for pages
# GLOBAL -----
## Titles -----
summaryTitle <- function (id) {
  ns <- NS(id)
  list(
    h3(textOutput(outputId = ns("summary_title")))
  )
}

summaryTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$summary_title <- renderText({paste0(str_c(data()$content, collapse = ", "))})
    }
  )
}

summaryListTitle <- function (id) {
  ns <- NS(id)
  list(
    h3(textOutput(outputId = ns("custom_list_title")))
  )
}

summaryListTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$custom_list_title <- renderText({paste0("Custom ", 
                                                     ifelse(data()$type == "gene", "Gene", "Cell Line"),
                                                     " List: ", 
                                                     str_c(data()$content, collapse = ", "))})
    }
  )
}

## Summary -----
summaryListText <- function (id) {
  ns <- NS(id)
  htmlOutput(outputId = ns("query_list_summary"))
}

summaryListTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$query_list_summary <- renderText({
        summary_list(input = data()) %>% 
          lit_linkr(summary_table = gene_summary) # see fun_helper.R
      })
    }
  )
}

# GENE -----
## Titles -----
geneTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title")))
  )
}

geneTitleServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))})
    }
  )
}
pathwayTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("pathway_title")))
  )
}

pathwayTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_title <- renderText({paste0("Pathway: ", summary_pathway(summary_table = pathways, input = data(), var = "pathway"), " (GO:", summary_pathway(summary_table = pathways, input = data(), var = "go"), ")")})
    }
  )
}

## Summary -----
geneText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    tags$dl(
      tags$dt("Gene Summary"), tags$dd(htmlOutput(outputId = ns("gene_summary_entrez_summary"))),
      tags$dt("Entrez ID"), tags$dd(htmlOutput(outputId = ns("ncbi_link"))), 
      tags$dt("Coding Sequence Length"), tags$dd(textOutput(outputId = ns("gene_length"))), 
      tags$dt("aka"), tags$dd(textOutput(outputId = ns("gene_summary_aka")))
    )
  )
}

geneTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))})
      output$gene_summary_approved_symbol <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"))
      output$gene_summary_approved_name <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))
      output$gene_summary_aka <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "aka"))
      output$gene_length <- renderText(paste0(summary_gene(summary_table = gene_location, input = data(), var = "cds_length"), " bp"))
      output$gene_summary_entrez_summary <- renderText({
        shiny::validate(
          need(!is.na(summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary")), "No data found for this gene."))
        summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary") %>% 
          lit_linkr(summary_table = gene_summary)})
      output$ncbi_link <- renderText(paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', 
                                            summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"),
                                            '" target="_blank">', 
                                            summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"),
                                            '</a>'))
    }
  )
}

#pathways
pathwayText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("pathway_title"))),
    tags$dl(
      tags$dt("Genes"), tags$dd(htmlOutput(outputId = ns("pathway_gene_symbols"))),
      tags$dt("Pathway Description"), tags$dd(textOutput(outputId = ns("pathway_def")))
    )
  )
}

pathwayTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_title <- renderText({paste0("Pathway: ", summary_pathway(summary_table = pathways, input = data(), var = "pathway"), " (GO:", summary_pathway(summary_table = pathways, input = data(), var = "go"), ")")})
      output$pathway_gene_symbols <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "data") %>% internal_link()})
      output$pathway_def <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "def")})
    }
  )
}

#protein page
#this is the protein tab for a gene search
proteinText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_summary_title"))),
    tags$dl(
      tags$dt("Protein Summary"), tags$dd(htmlOutput(outputId = ns("protein_summary_uniprot_summary"))), 
      tags$dt("Uniprot ID"), tags$dd(htmlOutput(outputId = ns("uniprot_link"))),
      tags$dt("Enzyme Commission", a(href = "https://en.wikipedia.org/wiki/Enzyme_Commission_number", img(src="link out_25.png", width="10", height="10"),  target="_blank")), tags$dd(htmlOutput(outputId = ns("ec_link"))),
      tags$dt("Protein Mass"), tags$dd(textOutput(outputId = ns("protein_summary_mass"))),
    ),
  )
}

proteinTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_summary_title <- renderText({paste0(summary_protein(summary_table = proteins, input = data(), var = "gene_name"), ": ", summary_protein(summary_table = proteins, input = data(), var = "protein_name"))})
      output$ec_link <- renderText(
        if (is.na(summary_protein(summary_table = proteins, input = data(), var = "ec"))) {paste0('<a href="https://enzyme.expasy.org" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        } else {paste0('<a href="https://enzyme.expasy.org/EC/', summary_protein(summary_table = proteins, input = data(), var = "ec"), '">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        })
      output$protein_summary_uniprot_summary <- renderText({
        shiny::validate(
          need(!is.na(summary_protein(summary_table = proteins, input = data(), var = "function_cc")), "No data found for this protein"))
        summary_protein(summary_table = proteins, input = data(), var = "function_cc") %>% 
          lit_linkr(summary_table = gene_summary)})
      output$uniprot_link <- renderText(paste0('<a href="https://www.uniprot.org/uniprot/', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"), '" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"),'</a>'))
      output$protein_summary_mass <- renderText(paste0(summary_protein(summary_table = proteins, input = data(), var = "mass"), " kDa"))
    }
  )
}

proteinSeq <- function (id) {
  ns <- NS(id)
  tagList(
    tags$dl(
      tags$dt("Protein Sequence"), tags$dd(verbatimTextOutput(outputId = ns("protein_summary_seq"))),
      tags$dt("Protein Length"), tags$dd(textOutput(outputId = ns("protein_length")))
    )
  )
}

proteinSeqServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_summary_seq <- renderText(summary_protein(summary_table = proteins, input = data(), var = "sequence"))
      output$protein_length <- renderText(glue::glue('{stringr::str_length(summary_protein(summary_table = proteins, input = data(), var = "sequence"))} AA'))
          }
  )
}

# CELL -----
## Titles -----
cellTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_summary_title")))
  )
}

cellTitleServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summary_title <- renderText({paste0(summary_cell(input = data(), var = "cell_line"), ": ", summary_cell(input = data(), var = "lineage_subtype"))})
    }
  )
}
lineageTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_title")))
  )
}

lineageTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_title <- renderText({paste0("Lineage: ", summary_cell(input = data(), var = "lineage"))})
    }
  )
}
lineageSubtypeTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_subtype_title")))
  )
}

lineageSubtypeTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_subtype_title <- renderText({paste0("Lineage: ", summary_cell(input = data(), var = "lineage_subtype"))})
    }
  )
}

tissueTitle <- function (id) {
  ns <- NS(id)
  list(
    h4(textOutput(outputId = ns("tissue_title")))
  )
}

tissueTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tissue_title <- renderText({paste0("Human anatogram of ", str_c(data()$content, collapse = ", "))})
    }
  )
}
###
## Summary -----

cellSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Lineage"), tags$dd(htmlOutput(outputId = ns("cell_summary_lineage"))),
      tags$dt("Lineage subtype"), tags$dd(htmlOutput(outputId = ns("cell_summary_lineage_subtype"))), 
      tags$dt("Description"), tags$dd(htmlOutput(outputId = ns("cell_description"))),
      tags$dt("Age"), tags$dd(textOutput(outputId = ns("cell_age"))),
      tags$dt("Sex"), tags$dd(textOutput(outputId = ns("cell_sex")))
    )
  )
}

cellSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summary_title <- renderText({
        data()$query
      })
      output$cell_summary_lineage <- renderText({
        summary_cell(input = data(), var = "lineage") %>% 
          map_chr(cell_linkr, type = "lineage")
      })
      output$cell_summary_lineage_subtype <- renderText({
        summary_cell(input = data(), var = "lineage_subtype") %>% 
          map_chr(cell_linkr, type = "lineage_subtype")
      })
      output$cell_description <- renderText({
        summary_cellosaurus(input = data(), var = "CC") %>% 
          lit_linkr(summary_table = gene_summary)
      })
      output$cell_age <- renderText({
        summary_cellosaurus(input = data(), var = "AG")
      })
      output$cell_sex <- renderText({
        summary_cellosaurus(input = data(), var = "SX")
      })
    })
}


lineageSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Cell lines"), tags$dd(textOutput(outputId = ns("lineage_summary_cell_lines"))),
    )
  )
}

lineageSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_summary_title <- renderText({
        data()$query
      })
      output$lineage_summary_cell_lines <- renderText({
        paste0(data()$content, collapse = ", ")
      })
    })
}

lineageSubtypeSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_subtype_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Cells"), tags$dd(textOutput(outputId = ns("lineage_subtype_summary_cell_lines"))),
    )
  )
}

lineageSubtypeSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_subtype_summary_title <- renderText({
        data()$query
      })
      output$lineage_subtype_summary_cell_lines <- renderText({
        paste0(data()$content, collapse = ", ")
      })
    })
}

cellListSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("title"))),
  )
}

cellListSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title <- renderText({
        data()$query
      })
    })
}



# COMPOUND -----
## Titles -----
compoundTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("compound_title")))
  )
}

compoundTitleServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_title <- renderText({
        paste0(summary_compound(summary_table = prism_meta, input = data(), var = "name") %>% stringr::str_to_title(), 
               ": ", 
               summary_compound(summary_table = prism_meta, input = data(), var = "moa") %>% stringr::str_to_title())
      })
    }
  )
}

metaboliteTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("metabolite_title")))
  )
}

metaboliteTitleServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$metabolite_title <- renderText({
        paste0(summary_metabolite(summary_table = hmdb_names, input = data(), var = "name") %>% stringr::str_to_title(),
               ": ",
               summary_metabolite(summary_table = hmdb_names, input = data(), var = "class") %>% stringr::str_to_title())
      })
    }
  )
}

moaTitle <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("moa_title")))
  )
}

moaTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$moa_title <- renderText({
        paste0("MOA: ", 
               summary_compound(summary_table = prism_meta, input = data(), var = "moa") %>% stringr::str_to_title())
      })
    }
  )
}
compoundListTitle <- function (id) {
  ns <- NS(id)
  list(
    h3(textOutput(outputId = ns("custom_compound_list")))
  )
}

compoundListTitleServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$custom_compound_list <- renderText({paste0("Custom Compound List: ", 
                                                        summary_compound_list(summary_table = prism_meta, input = data())  %>% stringr::str_to_title())
      })
    }
  )
}
###
## Summary -----
compoundSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("compound_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Description"), tags$dd(textOutput(outputId = ns("compound_summary_description"))),
      tags$dt("MOA"), tags$dd(htmlOutput(outputId = ns("compound_summary_moa"))),
      tags$dt("Phase"), tags$dd(textOutput(outputId = ns("compound_summary_phase"))),
      tags$dt("Disease Area"), tags$dd(textOutput(outputId = ns("compound_summary_diseasearea"))),
      tags$dt("Indication"), tags$dd(textOutput(outputId = ns("compound_summary_indication"))),
      tags$dt("Targets"), tags$dd(htmlOutput(outputId = ns("compound_summary_targets"))),
      tags$dt("CID"), tags$dd(htmlOutput(outputId = ns("cid_link")))
    )
  )
}

compoundSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_summary_title <- renderText({
        data()$query %>% 
          stringr::str_to_title()
      })
      output$compound_summary_description <- renderText({
        summary_compound(prism_meta, input = data(), var = "description") %>% 
          stringr::str_to_sentence()
      })
      output$compound_summary_moa <- renderText({
        summary_compound(prism_meta, input = data(), var = "moa") %>% 
          map_chr(moa_linkr)
      })
      output$compound_summary_phase <- renderText({
        summary_compound(prism_meta, input = data(), var = "phase") %>% 
          stringr::str_to_title()
      })
      output$compound_summary_diseasearea <- renderText({
        summary_compound(prism_meta, input = data(), var = "disease_area") %>% 
          stringr::str_to_title()
      })
      output$compound_summary_indication <- renderText({
        summary_compound(prism_meta, input = data(), var = "indication") %>% 
          stringr::str_to_title()
      })
      output$compound_summary_targets <- renderText({
        summary_compound(prism_meta, input = data(), var = "target") %>% 
          stringr::str_to_upper() %>% 
          map_chr(internal_link)
      })
      output$compound_summary_cid <- renderText({
        summary_compound(prism_meta, input = data(), var = "cid")
      })
      output$cid_link <- renderText(paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/', 
                                           summary_compound(summary_table = prism_meta, input = data(), var = "cid"),
                                           '" target="_blank">', 
                                           summary_compound(summary_table = prism_meta, input = data(), var = "cid"),
                                           '</a>')
      )
    })
}

## MOA QUERY -----
moaSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("moa_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("compounds"), tags$dd(textOutput(outputId = ns("moa_summary_compounds"))),
    )
  )
}

moaSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$moa_summary_title <- renderText({
        data()$query
      })
      output$moa_summary_compounds <- renderText({
        paste0(data()$content, collapse = ", ")
      })
    })
}

## CUSTOM COMPOUND LIST QUERY -----
compoundListSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("title"))),
  )
}

compoundListSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title <- renderText({
        data()$query
      })
    })
}

# GENE TAB -----
## Gene cards -----
geneGoTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("go_table_title")))
  )
}

geneGoTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$go_table_title <- renderText({glue::glue('GO pathways for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

## Protein cards -----
#this is the protein card titles
proteinSizeText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_size_title")))
  )
}

proteinSizeTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_size_title <- renderText({glue::glue('Size Information for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

proteinSeqText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_seq_title"))),
    tags$dl(
      tags$dt("Protein Sequence"), tags$dd(verbatimTextOutput(outputId = ns("protein_summary_seq"))),
    )
  )
}

proteinSeqTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_seq_title <- renderText({glue::glue('Sequence Information for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

proteinSignatureText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_signature_title")))
  )
}

proteinSignatureTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_signature_title <- renderText({glue::glue('Signature Information for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

proteinStructureText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_structure_title")))
  )
}

proteinStructureTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_structure_title <- renderText({glue::glue('Structural Information for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

## Expression cards -----
#this makes card titles
cellGeneExpressionTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_gene_table_text")))
  )
}

cellGeneExpressionTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_gene_table_text <- renderText({glue::glue('Gene expression table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

cellProteinExpressionPlotText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_protein_text")))
  )
}

cellProteinExpressionPlotTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_protein_text <- renderText({glue::glue('Protein expression plot for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

tissuePlotText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("tissueplot_text")))
  )
}

tissuePlotTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tissueplot_text <- renderText({glue::glue('Human tissue expression plot for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

tissueTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("tissuetable_text")))
  )
}

tissueTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tissuetable_text <- renderText({glue::glue('Human tissue expression table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

## Dependency cards -----
cellDepsLinPlotText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("dep_lin_title")))
  )
}

cellDepsLinPlotTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$dep_lin_title <- renderText({glue::glue('Dependency lineage plot for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

cellDepsSubLinPlotText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("dep_sublin_title")))
  )
}

cellDepsSubLinPlotTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$dep_sublin_title <- renderText({glue::glue('Dependency sublineage plot for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

cellDependenciesTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("dep_table_title")))
  )
}



cellDependenciesTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$dep_table_title <- renderText({glue::glue('Dependency table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

#DEPENDENCY CARD COESSENTIAL TAB TITLES
similarGenesTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("sim_text_title")))
  )
}

similarGenesTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$sim_text_title <- renderText({glue::glue('Similar gene table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

similarPathwaysTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("sim_pathways_text_title")))
  )
}

similarPathwaysTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$sim_pathways_text_title <- renderText({glue::glue('Similar pathway table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

dissimilarGenesTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("dsim_text_title")))
  )
}

dissimilarGenesTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$dsim_text_title <- renderText({glue::glue('Inverse gene table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

dissimilarPathwaysTableText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("dsim_pathways_text_title")))
  )
}

dissimilarPathwaysTableTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$dsim_pathways_text_title <- renderText({glue::glue('Inverse pathway table for {summary_protein(summary_table = proteins, input = data(), var = "gene_name")}')})
    }
  )
}

# CELL TAB -----

# COMPOUND TAB -----




