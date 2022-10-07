#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#FIND IMAGES -----
current_file_path <- here::here("cache", "protein_pdbs") #set in sh script
all_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)

#MAKE master dirs -----
target_file_path <- here::here("data", "images", "gene")
if(!dir.exists(here::here(target_file_path))) {dir.create(here::here(target_file_path), recursive = TRUE)}

#MOVE IMAGE LOOP -----
for (i in all_files) {
  #file_name to gene_name
  gene_name <- stringr::str_extract(i, "[[:graph:]]+(?=\\.)") #go from file name to gene_name
  image_name <- glue::glue('{gene_name}_structure_plot.jpg')
  pdb_name <- glue::glue('{gene_name}.pdb')
  
  #if file exists, skip
  if(file.exists(glue::glue('{target_file_path}/{gene_name}/{image_name}')) && file.exists(glue::glue('{target_file_path}/{gene_name}/{pdb_name}'))) {
    message(glue::glue('structure files for {gene_name} already exist'))
    unlink(glue::glue('{current_file_path}/{i}')) #deletes cached image, so it can counts remaining unmatched
    unlink(glue::glue('{current_file_path}/{pdb_name}'))
  } else {
    #make gene-specific dir
    save_path <- here::here(target_file_path, gene_name)
    if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
    #move image
    file.rename(from = here::here(current_file_path, i), 
                to = here::here(save_path, image_name))
    #move pdb file
    file.rename(from = here::here(current_file_path, pdb_name), 
                to = here::here(save_path, pdb_name))
    #update progress
    message(glue::glue('moved {gene_name} files'))
  }
}
#finish message
remaining_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)
message("all structure files moved")
message(glue::glue('these files remain: {glue::glue_collapse(remaining_files, ", ")}'))
