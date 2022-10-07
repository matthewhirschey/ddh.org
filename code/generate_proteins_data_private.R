#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#generate public protein file
source(here::here("code", "generate_proteins_data.R"))

#load data already generated, req'd below
expression_names <- readRDS(file=here::here("data", paste0(release, "_expression_names.Rds")))
expression_long <- readRDS(file=here::here("data", paste0(release, "_expression_long.Rds")))
expression_meta <- readRDS(file=here::here("data", paste0(release, "_expression_meta.Rds")))

#get protein datasets for CCLE from Gygi
tmp_proteins_ccle_info <- tempfile()
download.file(proteins_ccle_info_url, tmp_proteins_ccle_info)
proteins_ccle_info <- readxl::read_xlsx(tmp_proteins_ccle_info,  sheet = 2, col_names = TRUE)

tmp_proteins_ccle <- tempfile()
download.file(proteins_ccle_url, tmp_proteins_ccle)
proteins_ccle_raw <- readxl::read_xlsx(tmp_proteins_ccle, sheet = 2, col_names = TRUE)

#cleaning
proteins_ccle <-
  proteins_ccle_raw %>% 
  dplyr::select(-starts_with("Ten"), #remove 10plex code
                -c("Protein_Id", "Description", "Group_ID", "Uniprot", "Uniprot_Acc"), #remove extra metadata
                -contains("Column")) %>% #remove empty columns that were parsed
  `colnames<-`(str_remove_all(names(.), "\\_TenPx\\d+")) %>% 
  #dplyr::select(1:10) %>% #for testing
  tidyr::pivot_longer(cols = -Gene_Symbol, names_to = "ccle_code", values_to = "protein_expression") %>% 
  dplyr::filter(protein_expression != "NA", 
                !is.na(Gene_Symbol)) %>% 
  dplyr::mutate(ccle_code = stringr::str_remove_all(ccle_code, "\\_TenPx\\d+"), 
                protein_expression = as.numeric(protein_expression)) %>%
  dplyr::left_join(expression_meta, by = c("ccle_code" = "ccle_name")) %>% 
  dplyr::select(gene = Gene_Symbol, cell_line, protein_expression) %>% 
  dplyr::arrange(dplyr::desc(protein_expression), dplyr::desc(cell_line)) %>% #set-up distinct
  dplyr::distinct(cell_line, gene, .keep_all = TRUE) #remove a few cell line & gene_symbol duplicates
  
#combine with expression long
proteins_ccle_long <-  #capture here to make df for plotting
  proteins_ccle %>% 
  dplyr::left_join(expression_names) %>% 
  dplyr::select(X1, gene, protein_expression)

expression_long <- #update expression long with ccle protein data
  expression_long %>% 
  dplyr::full_join(proteins_ccle_long, by = c("X1", "gene"))

#save files
saveRDS(expression_long, file = here::here("data", paste0(release, "_expression_long.Rds"))) #updates expression long by adding to it

