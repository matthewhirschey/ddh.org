library(optparse)
library(tidyverse)
library(vroom)
library(janitor)
library(rentrez)
library(feather)

gene_names_url <- "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=md_eg_id&col=gd_locus_type&col=md_mim_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
pubtator_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
gene2pubmed_url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
gene_summary_output_filename <- "gene_summary.Rds"
fetch_batch_size <- 100   # batch size of 500 was too high. 200 worked, using 100 as reasonable default

#FOR TESTING: import gene_summary, and select original columns
#gene_summary <- readRDS(here::here("data", "gene_summary.Rds")) %>% 
#  select ("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary")

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

update_gene_summary <- function(gene_summary_input, gene2pubmed_url, pubtator_url) {
  gene_summary_input <- as_tibble(gene_summary_input)
  
  #gene2pubmedurl defined in current_release.R
  gene2pubmed_raw <- read_tsv(gene2pubmed_url, col_names = TRUE) %>% 
    clean_names()
  #add pubmed_count to gene_summary
  gene_summary_output <- gene2pubmed_raw %>% 
    filter(number_tax_id == 9606) %>%  #only the rows corresponding to humans (#tax_id = 9606) 
    group_by(gene_id) %>% 
    count(sort = TRUE) %>% 
    right_join(gene_summary_input, by = c("gene_id" = "ncbi_gene_id"))  %>% 
    rename("ncbi_gene_id" = "gene_id", 
           "pubmed_count" = "n") %>% 
    mutate(pubmed_count = replace_na(pubmed_count, 0)) %>% 
    ungroup()
  
  #add concept count
  gene2pubtator <- read_tsv(pubtator_url, col_names = c("pmid", "type", "concept_id", "mentions", "resource")) %>% 
    select(pmid, concept_id, mentions) %>% 
    filter(concept_id %in% gene_summary_input$ncbi_gene_id == TRUE)
  
  gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)
  
  gene_summary_output <- gene2pubtator %>% 
    group_by(concept_id) %>% 
    count(sort = TRUE) %>% 
    right_join(gene_summary_output, by = c("concept_id" = "ncbi_gene_id"))  %>% 
    rename("ncbi_gene_id" = "concept_id", 
           "concept_count" = "n") %>% 
    mutate(concept_count = replace_na(concept_count, 0)) %>% 
    ungroup()
  
  gene_summary_output <- gene_summary_output %>%
    mutate(pubmed_count_rank = percent_rank(pubmed_count),
           concept_count_rank = percent_rank(concept_count)) %>% 
    select("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary", "pubmed_count", "pubmed_count_rank", "concept_count", "concept_count_rank")
  return(gene_summary_output)
}

create_gene_summary <- function(gene_names_url, entrez_key, gene_summary_output_path) {
  gene_summary <- build_gene_summary(gene_names_url, entrez_key) 
  gene_summary <- update_gene_summary(gene_summary, gene2pubmed_url, pubtator_url)
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
