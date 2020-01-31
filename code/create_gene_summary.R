library(optparse)
library(tidyverse)
library(vroom)
library(janitor)
library(rentrez)
library(feather)

gene_names_url <- "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=md_eg_id&col=gd_locus_type&col=md_mim_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
gene_summary_output_filename <- "gene_summary.Rds"
fetch_batch_size <- 100   # batch size of 500 was too high. 200 worked, using 100 as reasonable default

# returns a data frame based on gene_names_url and summaries from entrez "gene" database.
build_gene_summary <- function(gene_names_url, entrez_key) {
  if (!is.null(entrez_key)) {
    set_entrez_key(entrez_key)
  }

  hugo <- vroom(gene_names_url, delim = "\t", col_names = TRUE) %>%
    clean_names()

  ids <- hugo %>%
    drop_na(ncbi_gene_id_supplied_by_ncbi) %>%
    pull(ncbi_gene_id_supplied_by_ncbi)

  # split the list of NCBI gene IDs into batches of size fetch_batch_size
  # so that entrez_summary can be called with batches instead of single IDs
  ids_batches <-split(ids, ceiling(seq_along(ids)/fetch_batch_size))

  entrez <- tibble(
    ncbi_gene_id_supplied_by_ncbi = numeric(),
    entrez_summary = character()
  )
  fetched_cnt <- 0
  # Fetch each batch
  for (id_batch in ids_batches) {
    summary_batch <- entrez_summary(db="gene", id=id_batch)
    # summary_batch is a list of summaries, so loop and handle each individually
    for (summary in summary_batch) {
      fetched_cnt <- fetched_cnt + 1
      tmp <- tibble(ncbi_gene_id_supplied_by_ncbi=as.numeric(summary$uid), entrez_summary=summary$summary)
      entrez <- entrez %>%
        bind_rows(tmp)
      }
    if (fetched_cnt %% 1000 == 0) {
      message("Fetched ", fetched_cnt)
    }
  }

  # Join hugo and the entrez fetched results
  gene_summary <- hugo %>%
    left_join(entrez) %>%
    rename(ncbi_gene_id = ncbi_gene_id_supplied_by_ncbi)

  gene_summary$entrez_summary <- replace_na(gene_summary$entrez_summary, "No NCBI summary")
  gene_summary$locus_type <- str_to_sentence(gene_summary$locus_type)
  gene_summary$locus_type <- str_replace(gene_summary$locus_type, "Rna", "RNA")
  gene_summary$approved_name <- str_to_sentence(gene_summary$approved_name)
  gene_summary %>%
    unite("aka", previous_symbols:alias_symbols, sep = ", ", na.rm = TRUE)
}

create_gene_summary <- function(gene_names_url, entrez_key, gene_summary_output_path) {
  gene_summary = build_gene_summary(gene_names_url, entrez_key)
  saveRDS(gene_summary, file = gene_summary_output_path)
}

# Command line argument parser that will let a user optionally specify:
# - entrez key to speed up fetching data. To create an entrez key see "Using API Keys" at https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# - alternate destination directory
option_list = list(
  make_option(c("--entrezkey"), type="character", default=NULL,
              help="NCBI entrez key [default= %default]", metavar="character"),
  make_option(c("--destdir"), type="character", default="data",
              help="Destination directory [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# fetch data and write gene summary feather file
create_gene_summary(gene_names_url, opt$entrezkey,
                    here::here(opt$destdir, gene_summary_output_filename))
