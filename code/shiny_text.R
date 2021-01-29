#Summary pages
#Gene summary
geneSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene Summary"), tags$dd(textOutput(outputId = ns("gene_summary_entrez_summary")))
    ),
  )
}

geneSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))})
      output$gene_summary_entrez_summary <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary"))
    })
}

#pathways
pathwaySummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("pathway_summary_title"))),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(textOutput(outputId = ns("pathway_summary_gene_symbols"))),
      tags$dt("Pathway Description"),
      tags$dd(textOutput(outputId = ns("pathway_summary_def")))
    ),
  )
}

pathwaySummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_summary_title <- renderText({paste0("Pathway: ", summary_pathway(summary_table = pathways, input = data(), var = "pathway"), " (GO:", summary_pathway(summary_table = pathways, input = data(), var = "go"), ")")})
      output$pathway_summary_gene_symbols <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "data")})
      output$pathway_summary_def <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "def")})
    }
  )
}

#gene_list
geneListSummaryText <- function (id) {
  ns <- NS(id)
  list(
    h3("Custom Gene List"),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(textOutput(outputId = ns("custom_gene_list")))
    )
  )
}

geneListSummaryTextServer <- function(id, data) { #what is data here?
  moduleServer(
    id,
    function(input, output, session) {
      output$custom_gene_list <- renderText({summary_gene_list(summary_table = gene_summary, input = data())})
    }
  )
}

#gene page
#this is the gene tab for a gene search
geneText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    tags$dl(
      tags$dt("Gene Summary"), tags$dd(textOutput(outputId = ns("gene_summary_entrez_summary"))),
      tags$dt("Entrez ID"), tags$dd(htmlOutput(outputId = ns("ncbi_link"))), 
      tags$dt("aka"), tags$dd(htmlOutput(outputId = ns("gene_summary_aka")))
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
      output$gene_summary_entrez_summary <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary"))
      output$ncbi_link <- renderText(paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"), '" target="_blank">', summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"),'</a>'))
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
      tags$dt("Protein Summary"), tags$dd(textOutput(outputId = ns("protein_summary_uniprot_summary"))), 
      tags$dt("Uniprot ID"), tags$dd(htmlOutput(outputId = ns("uniprot_link"))),
      tags$dt("Enzyme Commission", a(href = "https://en.wikipedia.org/wiki/Enzyme_Commission_number", img(src="link out_25.png", width="10", height="10"),  target="_blank")), tags$dd(htmlOutput(outputId = ns("ec_link"))),
      tags$dt("Protein Mass"), tags$dd(textOutput(outputId = ns("protein_summary_mass"))),
      tags$dt("Protein Sequence"), tags$dd(verbatimTextOutput(outputId = ns("protein_summary_seq"))),
    ),
  )
}

proteinTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_summary_title <- renderText({paste0(summary_protein(summary_table = proteins, input = data(), var = "gene_name"), ": ", summary_protein(summary_table = proteins, input = data(), var = "protein_name"))})
      #output$protein_summary_approved_symbol <- renderText(summary_protein(summary_table = proteins, input = data(), var = "gene_name"))
      #output$protein_summary_approved_name <- renderText(summary_protein(summary_table = proteins, input = data(), var = "protein_name"))
      output$ec_link <- renderText(
        if (is.na(summary_protein(summary_table = proteins, input = data(), var = "ec"))) {paste0('<a href="https://enzyme.expasy.org" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        } else {paste0('<a href="https://enzyme.expasy.org/EC/', summary_protein(summary_table = proteins, input = data(), var = "ec"), '">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        })
      output$protein_summary_uniprot_summary <- renderText(summary_protein(summary_table = proteins, input = data(), var = "function_cc"))
      output$uniprot_link <- renderText(paste0('<a href="https://www.uniprot.org/uniprot/', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"), '" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"),'</a>'))
      output$protein_summary_mass <- renderText(paste0(summary_protein(summary_table = proteins, input = data(), var = "mass"), " kDa"))
      output$protein_summary_seq <- renderText(summary_protein(summary_table = proteins, input = data(), var = "sequence"))
    }
  )
}

#Cell Lines
cellSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Lineage"), tags$dd(textOutput(outputId = ns("cell_summary_lineage"))),
      tags$dt("Lineage subtype"), tags$dd(textOutput(outputId = ns("cell_summary_lineage_subtype")))
    )
  )
}

cellSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summary_title <- renderText({
        data()$id
      })
      output$cell_summary_lineage <- renderText({
        summary_cell(expression_names, input = data(), var = "lineage")
      })
      output$cell_summary_lineage_subtype <- renderText({
        summary_cell(expression_names, input = data(), var = "lineage_subtype")
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
        data()$id
      })
      output$lineage_summary_cell_lines <- renderText({
        paste0(data()$cell_line, collapse = ", ")
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
        data()$id
      })
      output$lineage_subtype_summary_cell_lines <- renderText({
        paste0(data()$cell_line, collapse = ", ")
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
        data()$id
      })
    })
}

# compound
compoundSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("compound_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("MOA"), tags$dd(textOutput(outputId = ns("compound_summary_moa"))),
      tags$dt("CID"), tags$dd(textOutput(outputId = ns("compound_summary_cid")))
    )
  )
}

compoundSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_summary_title <- renderText({
        data()$id
      })
      output$compound_summary_moa <- renderText({
        summary_compound(prism_names, input = data(), var = "moa")
      })
      output$compound_summary_cid <- renderText({
        summary_compound(prism_names, input = data(), var = "cid")
      })
    })
}

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
        data()$id
      })
      output$moa_summary_compounds <- renderText({
        paste0(data()$compound, collapse = ", ")
      })
    })
}


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
        data()$id
      })
    })
}


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
