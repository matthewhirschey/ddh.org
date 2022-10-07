# read ENTREZ_KEY from environment variable
entrez_key = Sys.getenv("ENTREZ_KEY")

gene_names_url <- "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=md_eg_id&col=gd_locus_type&col=md_mim_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
source(here::here("code", "current_release.R"))
gene_summary_output_filename <- "gene_summary.Rds"
fetch_batch_size <- 100   # batch size of 500 was too high. 200 worked, using 100 as reasonable default

#FOR TESTING: import gene_summary, and select original columns
#gene_summary <- readRDS(here::here("data", "gene_summary.Rds")) %>% 
#  select ("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary")

retry_entrez_summary <- function(db, id, retries=entrez_retries, retry_sleep_seconds=entrez_retry_sleep_seconds) {
  result <- NULL
  for (i in 0:retries) {
    tryCatch({
      result <- entrez_summary(db=db, id=id)
    }, error = function(err_msg) {
      message(err_msg)
    })
    if (!is.null(result)) {
      break
    }
  }
  if (is.null(result)) {
    stop("Unable to fetch entrez_summary for ", paste0(id, collapse=" "))
  }
  result
}

# returns a data frame based on gene_names_url and summaries from entrez "gene" database.
build_gene_summary <- function(gene_names_url, entrez_key, gene_symbol = NULL) {
  if (!is.null(entrez_key)) {
    set_entrez_key(entrez_key)
  }
  
  if (is.null(gene_symbol)) {
    hugo <- vroom::vroom(gene_names_url, delim = "\t", col_names = TRUE) %>%
      janitor::clean_names()
    
    ids <- hugo %>%
      tidyr::drop_na(ncbi_gene_id_supplied_by_ncbi) %>%
      dplyr::pull(ncbi_gene_id_supplied_by_ncbi)
    
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
      summary_batch <- retry_entrez_summary(db="gene", id=id_batch)
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
    
  } else {
    #this allows gene_symbol to be specified for methods generation in create_methods.R
    hugo <- vroom::vroom(gene_names_url, delim = "\t", col_names = TRUE) %>%
      janitor::clean_names() %>% 
      dplyr::filter(approved_symbol %in% gene_symbol)
    
    gene_symbol_id <- hugo$ncbi_gene_id_supplied_by_ncbi
    
    summary <- retry_entrez_summary(db = "gene", id = gene_symbol_id, retries = 1, retry_sleep_seconds = 1)
    entrez <- tibble::tibble(ncbi_gene_id_supplied_by_ncbi=as.numeric(summary$uid), entrez_summary=summary$summary)
  }
  
  # Join hugo and the entrez fetched results
  gene_summary <- hugo %>%
    left_join(entrez) %>%
    rename(ncbi_gene_id = ncbi_gene_id_supplied_by_ncbi)
  
  gene_summary$entrez_summary <- tidyr::replace_na(gene_summary$entrez_summary, "No NCBI summary")
  gene_summary$locus_type <- stringr::str_to_sentence(gene_summary$locus_type)
  gene_summary$locus_type <- stringr::str_replace(gene_summary$locus_type, "Rna", "RNA")
  gene_summary$approved_name <- stringr::str_to_sentence(gene_summary$approved_name)
  gene_summary %>%
    tidyr::unite("aka", previous_symbols:alias_symbols, sep = ", ", na.rm = TRUE)
}

update_gene_summary <- function(gene_summary_input, 
                                pubmed_count_url = gene2pubmed_url, 
                                concept_count_url = pubtator_url) {
  gene_summary_input <- as_tibble(gene_summary_input)
  
  #gene2pubmedurl defined in current_release.R
  gene2pubmed_raw <- read_tsv(pubmed_count_url, col_names = TRUE) %>% 
    janitor::clean_names()
  #add pubmed_count to gene_summary
  gene_summary_output <- gene2pubmed_raw %>% 
    dplyr::filter(number_tax_id == 9606) %>%  #only the rows corresponding to humans (#tax_id = 9606) 
    dplyr::group_by(gene_id) %>% 
    dplyr::count(sort = TRUE) %>% 
    dplyr::right_join(gene_summary_input, by = c("gene_id" = "ncbi_gene_id"))  %>% 
    dplyr::rename("ncbi_gene_id" = "gene_id", 
           "pubmed_count" = "n") %>% 
    dplyr::mutate(pubmed_count = replace_na(pubmed_count, 0)) %>% 
    dplyr::ungroup()
  
  #add concept count
  gene2pubtator <- read_tsv(concept_count_url, col_names = c("pmid", "type", "concept_id", "mentions", "resource")) %>% 
    dplyr::select(pmid, concept_id, mentions) %>% 
    dplyr::filter(concept_id %in% gene_summary_input$ncbi_gene_id == TRUE)
  
  gene2pubtator$concept_id <- as.numeric(gene2pubtator$concept_id)
  
  gene_summary_output <- gene2pubtator %>% 
    dplyr::group_by(concept_id) %>% 
    dplyr::count(sort = TRUE) %>% 
    dplyr::right_join(gene_summary_output, by = c("concept_id" = "ncbi_gene_id"))  %>% 
    dplyr::rename("ncbi_gene_id" = "concept_id", 
           "concept_count" = "n") %>% 
    dplyr::mutate(concept_count = replace_na(concept_count, 0)) %>% 
    dplyr::ungroup()
  
  gene_summary_output <- gene_summary_output %>%
    dplyr::mutate(pubmed_count_rank = percent_rank(pubmed_count),
           concept_count_rank = percent_rank(concept_count)) %>% 
    dplyr::select("approved_symbol", "approved_name", "aka", "ncbi_gene_id", "hgnc_id", "chromosome", "ref_seq_i_ds", "locus_type", "omim_id_supplied_by_omim", "uni_prot_id_supplied_by_uni_prot", "entrez_summary", "pubmed_count", "pubmed_count_rank", "concept_count", "concept_count_rank")
  return(gene_summary_output)
}

generate_gene_summary <- function(gene_names_url, entrez_key, gene_summary_output_path, gene_symbol = NULL, regenerate = TRUE) {
  if (regenerate == FALSE) {
    #get Rds from zenodo for saved gene_summary information
    tmp_gene_location <- tempfile()
    download.file(gene_summary_url, tmp_gene_location)
    gene_summary <- readRDS(tmp_gene_location)
    saveRDS(gene_summary, file = gene_summary_output_path)
  } else {
    # Specify an entrez key to speed up fetching data. To create an entrez key see "Using API Keys" at https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
    gene_summary <- build_gene_summary(gene_names_url, entrez_key, gene_symbol) 
    gene_summary <- update_gene_summary(gene_summary, gene2pubmed_url, pubtator_url)
    saveRDS(gene_summary, file = gene_summary_output_path)
  }
}

#testing for full gene set (top) or single gene for methods (bottom); these are now called in generate_data or create_methods
#generate_gene_summary(gene_names_url, entrez_key, here::here("data", paste0(release, "_", gene_summary_output_filename)), gene_symbol = NULL)
#generate_gene_summary(gene_names_url, entrez_key, here::here("data", paste0(release, "_", gene_summary_output_filename)), gene_symbol = "TP53")

