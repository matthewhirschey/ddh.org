#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

#FIND IMAGES -----
current_file_path <- here::here("cache", "barcodes") #manually set to cache
all_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)

#MAKE master dirs -----
target_file_path <- here::here("data", "images", "gene")
if(!dir.exists(here::here(target_file_path))) {dir.create(here::here(target_file_path), recursive = TRUE)}

#MOVE IMAGE LOOP -----
for (i in all_files) {
  #file_name to gene_name
  gene_name <- stringr::str_extract(i, "[[:graph:]]+(?=_barcode)") #go from file name to gene_name
  save_name <- i
  
  #if file exists, skip
  if(file.exists(glue::glue('{target_file_path}/{gene_name}/{save_name}'))) {
    message(glue::glue('barcode for {gene_name} already exists'))
    unlink(glue::glue('{current_file_path}/{i}')) #deletes cached image, so it can counts remaining unmatched
  } else {
    #make gene-specific dir
    save_path <- here::here(target_file_path, gene_name)
    if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
    
    #move image
    file.rename(from = here::here(current_file_path, i), 
                to = here::here(save_path, save_name))
    
    #update progress
    message(glue::glue('moved {i}'))
  }
}
#finish message
remaining_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)
message("all barcodes moved")
message(glue::glue('{length(remaining_files)} files remain'))
