#LOAD IMAGES.ZIP-----
#load current data/images/zip file in first step of step7 sh script

#READ CURRENT RELEASE information-----
source(here::here("code", "current_release.R"))

#LOAD LIBRARIES-----
source(here::here("code", "install_libraries.R"))

#set plan
num_workers = strtoi(Sys.getenv("SLURM_CPUS_PER_TASK", "4")) # Use Slurm setting or fall back to 4
if(supportsMulticore()) {
  plan(multicore, workers=num_workers)  
} else {
  plan(multisession, workers=num_workers)
}

#SET test VARS-----
test <- FALSE #set to TRUE for test dataset only, and to bypass futures
purge <- FALSE #set to TRUE when updating data on quarterly release; set to FALSE everyother time (when cards don't need regeneration)

#LOAD DATA-----
source(here::here("code", "app_params.R")) #needed for app_data to load correctly
source(here::here("code", "app_data.R")) #not best practice, but easy way to get all the data; REPLACE

#LOAD STRUCTURES INTO CACHE
#zip is loaded into cache at the end of step6
structure_tar_path <- here::here(app_cache_dir, "protein_pdbs.tar.gz")
structure_path <- here::here(app_cache_dir, "protein_pdbs")
untar(tarfile = structure_tar_path, 
      exdir = structure_path)

#LOAD BARCODES INTO CACHE
#zip is loaded into cache at the end of step6
barcode_zip_path <- here::here(app_cache_dir, "barcodes.zip")
barcode_path <- here::here(app_cache_dir)
unzip(zipfile = barcode_zip_path, 
      exdir = barcode_path)

#LOAD CELL IMAGES INTO CACHE
#zip is loaded into cache at the end of step6
cell_images_zip_path <- here::here(app_cache_dir, "cell_images.zip")
cell_images_path <- here::here(app_cache_dir)
unzip(zipfile = cell_images_zip_path, 
      exdir = cell_images_path)

#LOAD FUNCTIONS-----
source(here::here("code", "fun_colors.R"), local = TRUE)
source(here::here("code", "fun_helper.R"), local = TRUE)
source(here::here("code", "fun_plots.R"))

#REMOVE DATA-----
if(purge == TRUE){
all_files <- list.files(here::here("data", "images"), full.names = T, recursive = TRUE)
#list files that don't need to be updated
files_to_keep <- c("barcode",
                   "ideogram", 
                   "size",
                   "sequence",
                   "structure",
                   "molecule_structure")
files_to_remove <- list(all_files[str_detect(all_files, paste(files_to_keep,collapse = "|"), negate = TRUE)]) #paste converts to regex
walk(.x = files_to_remove, .f = unlink)
}

#SAVE GENES-----
#single function
save_gene <- function(gene_symbol,
                      fun_name, 
                      type = "gene", 
                      plot_type = "card") {
  #make master dirs
  if(!dir.exists(here::here("data", "images", type))) {dir.create(here::here("data", "images", type), recursive = TRUE)}
  
  #make gene-specific dir
  save_path <- here::here("data", "images", type, gene_symbol)
  if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
  
  #file prefix
  file_name <- glue::glue('{save_path}/{gene_symbol}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg')

  #if file exists, skip
  if(file.exists(file_name)) {
    message(glue::glue('{stringr::str_remove(fun_name, "make_")} {plot_type} for {gene_symbol} already exists'))
  } else {
    #build the input for the raw fun() & call function
    if(plot_type == "card"){
      image <- rlang::exec(fun_name, 
                           input = list(type = type, query = gene_symbol),
                           card = TRUE)
      image_width = 1000 
      image_height = 1500
      dpi = "retina"
    } else {
      #get file size
      plot_size <- plot_size_finder(function_name = fun_name)
      
      #if-else here b/c function won't like getting a card var if not there
      image <- rlang::exec(fun_name, 
                           input = list(type = type, query = gene_symbol))
      image_width = plot_size$plot_width
      image_height = plot_size$plot_height
      dpi = "print"
    }
    
    ggsave(filename = glue::glue('{gene_symbol}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg'),
           plot = image, 
           width = image_width, 
           height = image_height, 
           units = "px",
           path = save_path, 
           device = "jpeg",
           dpi = dpi, 
           bg = "white")
    
    message(glue::glue('saved the {stringr::str_remove(fun_name, "make_")} {plot_type} for {gene_symbol}'))
  }
}

#test single case
# save_gene(gene_symbol = "ROCK1", fun_name = "make_ideogram")
# save_gene(gene_symbol = "ROCK1", fun_name = "make_ideogram", plot_type = "plot")
# save_gene(gene_symbol = "ROCK1", fun_name = "make_celldeps")
# save_gene(gene_symbol = "ROCK1", fun_name = "make_sequence")
# save_gene(gene_symbol = "empty", fun_name = "make_ideogram")

#paramaterize function
save_gene_images <- function(gene_symbol, 
                             fun_name, 
                             test, 
                             image_type = "card") {
  #https://vctrs.r-lib.org/reference/vec-rep.html
  genes <- vctrs::vec_rep(gene_symbol, length(fun_name)) #rep along the vec
  functions <- vctrs::vec_rep_each(fun_name, length(gene_symbol)) #rep along the vec grouping by 'each'
  
  #fix futures
  #https://cran.r-project.org/web/packages/future/vignettes/future-4-issues.html
  
  if(test == TRUE){
    purrr::walk2(.x = genes, .y = functions, ~ save_gene(format_path_part(.x), .y, plot_type = image_type))
  } else {
    furrr::future_walk2(.x = genes, .y = functions, ~ save_gene(format_path_part(.x), .y, plot_type = image_type))
  }
}

#test
# save_gene_card(gene_symbol = "ROCK2", fun_name = "make_ideogram")
# save_gene_card(gene_symbol = "C10orf113", fun_name = "make_ideogram")
# save_gene_cards(gene_symbol = gene_vec, fun_name = fun_name_vec)

#define functions to make
fun_name_genes <- c("make_ideogram", "make_colorful_anatogram", "make_celldeps")

##Define genes to get ----
if(test == TRUE) {
  gene_vec <- c("ROCK1", "ROCK2") #sample(gene_summary$approved_symbol, 10) #test
} else {
  #dont grab full....#gene_vec <- gene_summary$approved_symbol
  #grab top n most popular genes
  n_genes <- 1000
  gene_vec <- 
    read_tsv(gene2pubmed_url, col_names = TRUE) %>% #gene2pubmed_url read from current_release.R
    janitor::clean_names() %>%
    dplyr::filter(number_tax_id == 9606) %>%  #only the rows corresponding to humans (#tax_id = 9606) 
    dplyr::group_by(gene_id) %>% 
    dplyr::count(sort = TRUE) %>% 
    dplyr::left_join(gene_summary, by = c("gene_id" = "ncbi_gene_id")) %>% 
    dplyr::select(gene_id, approved_symbol, number_pubs = n, approved_name) %>% 
    dplyr::filter(!is.na(approved_symbol))  %>% 
    dplyr::ungroup() %>% 
    dplyr::slice_max(number_pubs, n = n_genes) %>% 
    dplyr::pull(approved_symbol)
  
  #filter top 1000 for presence in achilles #!MIR21
  gene_vec <- gene_vec[gene_vec %in% achilles_cor_nest$fav_gene]
  
  #save
  saveRDS(gene_vec, file = here::here("data", "gene_vec.Rds"))
  
  #add test genes
  # if(file.exists(here::here("tests", "data", "all_genes.Rds"))){ #file is generated AFTER this step...
  #   all_genes <- readRDS(here::here("tests", "data", "all_genes.Rds"))
  #   gene_vec <- c(gene_vec, all_genes)
  # }
}

##Save Gene Plots --------
fun_name_gene_plots <- c("make_ideogram",
                         "make_female_anatogram",
                         "make_male_anatogram",
                         "make_celldeps", 
                         # "make_radial", #toggle-dependent 
                         "make_proteinsize", 
                         "make_pubmed", 
                         "make_cellanatogram", 
                         "make_tissue",
                         "make_cellbins", 
                         "make_lineage", 
                         "make_sublineage", 
                         "make_correlation", 
                         "make_cellexpression", 
                         "make_expdep") 
                         

save_gene_images(gene_symbol = gene_vec, 
                 fun_name = fun_name_gene_plots, 
                 test, 
                 image_type = "plot")

##Special Save Protein Structure --------
source(here::here("code", "generate_structures.R"))

##Special Save Gene Barcode --------
source(here::here("code", "generate_barcodes.R"))

##Save Gene Cards ----
#each needs card var in function name
fun_name_gene_cards <- c("make_ideogram", 
                         "make_female_anatogram", #only female
                         #"make_male_anatogram",
                         "make_celldeps", 
                         "make_sequence", #only card
                         "make_radial", 
                         "make_proteinsize", 
                         "make_structure", #only card, plot already exists as jpeg
                         "make_pubmed", 
                         "make_cellanatogram", 
                         "make_tissue",
                         "make_cellbins", 
                         "make_lineage", 
                         "make_sublineage", 
                         "make_cellexpression", 
                         "make_cellgeneprotein", 
                         "make_expdep")

save_gene_images(gene_symbol = gene_vec, 
                 fun_name = fun_name_gene_cards, 
                 test, 
                 image_type = "card")

#SAVE CELL LINES-----
#single function
#single function
save_cell <- function(cell_line,
                      fun_name, 
                      type = "cell", 
                      plot_type = "card") {
  
  #make master dirs
  if(!dir.exists(here::here("data", "images", type))) {dir.create(here::here("data", "images", type), recursive = TRUE)}
  
  #make cell-specific dir
  save_path <- here::here("data", "images", type, cell_line)
  if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
  
  #file prefix
  file_name <- glue::glue('{save_path}/{cell_line}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg')
  
  #if file exists, skip
  if(file.exists(file_name)) {
    message(glue::glue('{stringr::str_remove(fun_name, "make_")} {plot_type} for {cell_line} already exists'))
  } else {
    #build the input for the raw fun() & call function
    if(plot_type == "card"){
      image <- rlang::exec(fun_name, 
                           input = list(type = type, cell_line = cell_line),
                           card = TRUE)
      image_width = 1000 
      image_height = 1500
    } else {
      #if-else here b/c function won't like getting a card var if not there
      image <- rlang::exec(fun_name, 
                           input = list(type = type, cell_line = cell_line))
      image_width = NA
      image_height = NA
    }
    
    ggsave(filename = glue::glue('{cell_line}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg'),
           plot = image, 
           width = image_width, 
           height = image_height, 
           units = "px",
           path = save_path, 
           device = "jpeg",
           dpi = "print")
    
    message(glue::glue('saved the {stringr::str_remove(fun_name, "make_")} {plot_type} for {cell_line}'))
  }
}

#test single case
# save_cell(cell_line = "HEPG2", fun_name = "make_celldeps", plot_type = "plot")
# save_cell(cell_line = "HEPG2", fun_name = "make_celldeps", plot_type = "card")

#paramaterize function
save_cell_images <- function(cell_line, 
                             fun_name, 
                             test, 
                             image_type = "card") {
  #https://vctrs.r-lib.org/reference/vec-rep.html
  cells <- vctrs::vec_rep(cell_line, length(fun_name)) #rep along the vec
  functions <- vctrs::vec_rep_each(fun_name, length(cell_line)) #rep along the vec grouping by 'each'
  
  #fix futures
  #https://cran.r-project.org/web/packages/future/vignettes/future-4-issues.html
  
  if(test == TRUE){
    purrr::walk2(.x = cells, .y = functions, ~ save_cell(format_path_part(.x), .y, plot_type = image_type))
  } else {
    furrr::future_walk2(.x = cells, .y = functions, ~ save_cell(format_path_part(.x), .y, plot_type = image_type))
  }
}

#define cells to get
if(test == TRUE) {
  cell_vec <- c("HEPG2", "HUH7") #test
} else {
  cell_vec <- expression_meta$cell_line #full
}

##Save Cell Plots --------
fun_name_cell_plots <- c("make_celldeps")

save_cell_images(cell_line = cell_vec, 
                 fun_name = fun_name_cell_plots, 
                 test, 
                 image_type = "plot")

##Special Save Cell Images--------
source(here::here("code", "generate_cell_images.R"))

##Save Cell Cards ----
#each needs card var in function name
fun_name_cell_cards <- c("make_celldeps")

save_cell_images(cell_line = cell_vec, 
                 fun_name = fun_name_cell_cards, 
                 test, 
                 image_type = "card")

#SAVE COMPOUNDS-----
#single function to save images using compound input
save_compound <- function(compound,
                          fun_name, 
                          type = "compound", 
                          plot_type = "card") {
  
  #make master dirs
  if(!dir.exists(here::here("data", "images", type))) {dir.create(here::here("data", "images", type), recursive = TRUE)}
  
  #make gene-specific dir
  save_path <- here::here("data", "images", type, compound)
  if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
  
  #file prefix
  file_name <- glue::glue('{save_path}/{compound}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg')
  
  #if file exists, skip
  if(file.exists(file_name)) {
    message(glue::glue('{stringr::str_remove(fun_name, "make_")} {plot_type} for {compound} already exists'))
  } else {
    #build the input for the raw fun() & call function
    if(plot_type == "card"){
      image <- rlang::exec(fun_name, 
                           input = list(type = type, compound = compound),
                           card = TRUE)
      image_width = 1000 
      image_height = 1500
    } else {
      #if-else here b/c function won't like getting a card var if not there
      image <- rlang::exec(fun_name, 
                           input = list(type = type, compound = compound))
      image_width = NA
      image_height = NA
    }
    
    #special case for saving non-ggplot image
    if(fun_name == "make_molecule_structure"){
      magick <- image_read(image)
      #CARD
      magick::image_write(image = magick, 
                          path = glue::glue('{save_path}/{compound}_{stringr::str_remove(fun_name, "make_")}_{plot_type.jpeg'),
                          format = "jpeg")
    } else {
      ggsave(filename = glue::glue('{compound}_{stringr::str_remove(fun_name, "make_")}_{plot_type}.jpeg'),
             plot = image, 
             width = image_width, 
             height = image_height, 
             units = "px",
             path = save_path, 
             device = "jpeg",
             dpi = "print")
    }
    message(glue::glue('saved the {stringr::str_remove(fun_name, "make_")} {plot_type} for {compound}'))
  }
}

#test single case
# save_compound(compound = "aspirin", fun_name = "make_celldeps")
# save_compound(compound = "aspirin", fun_name = "make_molecule_structure")

save_compound_images <- function(compound, 
                                 fun_name, 
                                 test, 
                                 image_type = "card") {
  #https://vctrs.r-lib.org/reference/vec-rep.html
  compounds <- vctrs::vec_rep(compound, length(fun_name)) #rep along the vec
  functions <- vctrs::vec_rep_each(fun_name, length(compound)) #rep along the vec grouping by 'each'
  
  #fix futures
  #https://cran.r-project.org/web/packages/future/vignettes/future-4-issues.html
  
  if(test == TRUE){
    purrr::walk2(.x = compounds, .y = functions, ~ save_compound(format_path_part(.x), .y, plot_type = image_type))
  } else {
    furrr::future_walk2(.x = compounds, .y = functions, ~ save_compound(format_path_part(.x), .y, plot_type = image_type))
  }
}

##Define compounds to get ----
if(test == TRUE) {
  compound_vec <- c("aspirin", "ibuprofen") #test
} else {
  compound_vec <- prism_names$name #full
}

##Save Compound Cards ----
#each needs card var in function name
fun_name_compound_cards <- c("make_celldeps", 
                             #"make_molecule_structure", 
                             "make_pubmed")


save_compound_images(compound = compound_vec, 
                     fun_name = fun_name_compound_cards, 
                     test, 
                     image_type = "card")

##Save Compound Plots ----
fun_name_compound_plots <- c("make_celldeps", 
                             #"make_molecule_structure", 
                             "make_pubmed")


save_compound_images(compound = compound_vec, 
                     fun_name = fun_name_compound_plots, 
                     test, 
                     image_type = "plot")



# MESSAGE -----
message("finished generating all cards")
