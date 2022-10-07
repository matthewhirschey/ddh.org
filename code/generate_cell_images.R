# To remake images
# 1. uncomment all commented code
# 2. navigate to cache and unzip cell_images_raw
# 3. load cell_image_cropper() which takes the raw images and crops them for plots
# 4. get all_cell_images and purrr::walk() to crop them
# 5. reset default card == TRUE in cell_image_cropper
# 6. run purrr::walk again to make cards
# 7. recomment all this code

# original code to scrape ATCC for images; not image URLs -----
# scrape <- FALSE
# if(scrape == TRUE) {
# 
#   get_cell_url <- function(atcc_id) {
#     message(glue::glue("getting {atcc_id}"))
#     if(is.na(atcc_id)) {return(NA)}
#     atcc_url <- glue::glue('https://www.atcc.org/products/all/{atcc_id}.aspx#characteristics')
# 
#     cell_urls <-
#       read_html(atcc_url) %>%
#       html_elements("img") %>%
#       html_attr("src") %>%
#       purrr::keep(str_detect(., pattern = "micrograph"))
# 
#     cell_url <- cell_urls[1]
#     if(is.na(cell_url)){return(NA)}
# 
#     base_url <- "https://www.atcc.org"
#     image_url <- glue::glue('{base_url}{cell_url}')
# 
#     return(image_url)
#   }
# 
#   ##for testing
#   # get_cell_url(atcc_id = "CCL-2")
#   # get_cell_url(atcc_id = "HB-8065")
#   # get_cell_url("CRL-8024")
#   # get_cell_url(atcc_id = "HB-8064")
#   # get_cell_url(atcc_id = "HTB-108")
# 
#   get_cell_url_possibly <- possibly(get_cell_url, otherwise = NA)
# 
#   get_cell_image <- function(atcc_id){
# 
#     if(is.na(atcc_id)){
#       return(glue::glue("No ATCC id for {atcc_id}"))
#     } else {
#       url <- get_cell_url_possibly(atcc_id)
#       if(is.na(url)){return(glue::glue("No image retrieved for {atcc_id}"))}
#       file_name <- glue::glue('{atcc_id}.jpg')
#       download.file(url = url,
#                     destfile = here::here("data", "cells", file_name), #shouldve named it for cell line name
#                     quiet = TRUE)
#       print(glue::glue("finished {atcc_id}"))
#       Sys.sleep(.05)
#       return(file_name)
#     }
#   }
# 
#   #alt approach to download files
#   #uses some of the functions in the ifelse statement
#   #the unmatched `get_these` do not have images at atcc
#   all_ids <- cellosaurus %>% filter(!is.na(ATCC)) %>% pull(ATCC)
#   downloads <- list.files(path = here::here("data", "cells")) %>%
#     stringr::str_remove(., pattern = "\\.jpg")
#   get_these <- all_ids[!(all_ids %in% downloads)]
# 
# 
#   for (i in get_these[1:length(get_these)]) {
#     file_name <- glue::glue("{i}.jpg")
#       url <- get_cell_url_possibly(atcc_id = i)
#       if(!is.na(url)){
#       download.file(url = url,
#                     destfile = here::here("data", "cells", file_name),
#                     quiet = TRUE)
#       } else {
#         print(glue::glue("cannot get {file_name}"))
#       }
#   }
# 
#   cell_rename <- function(file_name){
#     atcc_id <- stringr::str_remove(file_name, pattern = "\\.jpg")
#     cell_name <-
#       cellosaurus %>%
#       dplyr::filter(ATCC %in% atcc_id) %>%
#       pull(name)
#     if(length(cell_name) == 0){return(glue::glue("didn't rename {file_name}"))}
#     new_file_name <- glue::glue("{cell_name}.jpg")
# 
#     file.rename(from = here::here("data", "cells", file_name),
#                 to = here::here("cache", "cell_images", new_file_name)) #also move it here
#     message(glue::glue("renamed {new_file_name}"))
#   }
# 
#   downloads <- list.files(path = here::here("data", "cells"))
#   purrr::walk(.x = downloads, .f = cell_rename)
# 
#   ##fix images
#   #library(magick)
# 
# cell_image_cropper <- function(file_name,
#                                card = FALSE){ #change to true for purrr::walk() below
#   path <- here::here("cache", "cell_images_raw")
#   cell_name <- stringr::str_remove(file_name, pattern = "\\.jpg")
#   cell <- image_read(glue::glue("{path}/{file_name}")) #
#   #print(cell)
# 
#   zoom_factor <- "x600"
#   cell_image <-
#     image_trim(cell) %>%
#     image_scale(glue::glue("{zoom_factor}")) %>%
#     image_rotate(90)
# 
#   cell_info <- image_info(cell_image)
# 
#   image_size <- "400x400"#400x400 square
#   top <-
#     cell_image %>%
#     image_crop(glue::glue("{image_size}+120+10")) #120 from left, 10 from top
# 
#   bottom_anchor <- cell_info$height - 400
#   bottom <-
#     cell_image %>%
#     image_crop(glue::glue("{image_size}+120+{bottom_anchor}"))
# 
#   trim_scale <- function(img){
#     tmp <-
#       img %>%
#       image_trim() %>% #add extra to get rid of white space
#       image_scale("400") %>%
#       image_crop("400x400") %>%
#       image_trim() %>%
#       image_scale("x400") %>%
#       image_crop("400x400") %>%
#       image_border(color = "white", geometry = "3x3")
#     return(tmp)
#   }
#   top <- trim_scale(top)
#   bottom <- trim_scale(bottom)
# 
#   cell_stack <- c(top, bottom)
#   cell_final <- image_append(cell_stack, stack = TRUE)
#   #print(cell_final)
# 
#   new_path <- here::here("cache", "cell_images")
# 
#   if(card == TRUE){
#     image_write(top,
#                 path = glue::glue("{new_path}/{cell_name}_cell_image_card.jpeg")) #jpeg! not jpg
#     message(glue::glue("finished {cell_name} card"))
#   } else {
#     image_write(cell_final,
#                 path = glue::glue("{new_path}/{cell_name}_cell_image_plot.jpeg")) #jpeg! not jpg
#     message(glue::glue("finished {cell_name} plot"))
#   }
# }
# all_cell_images <- list.files(path = here::here("cache", "cell_images_raw"))
# purrr::walk(.x = all_cell_images, .f = cell_image_cropper)
# #stringr::str_which(all_cell_images, "T84")
# 
# }

#FIND IMAGES -----
current_file_path <- here::here("cache", "cell_images") #manually set to cache
all_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)

#MAKE master dirs -----
target_file_path <- here::here("data", "images", "cell")
if(!dir.exists(here::here(target_file_path))) {dir.create(here::here(target_file_path), recursive = TRUE)}

#MOVE IMAGE LOOP -----
for (i in all_files) {
  #file_name to gene_name
  cell_name <- stringr::str_extract(i, "[[:graph:]]+(?=_cell_image)") #go from file name to cell_name for either plot or card

  #if file exists, skip
  if(file.exists(glue::glue('{target_file_path}/{cell_name}/{i}'))) {
    message(glue::glue('image for {cell_name} already exists'))
    unlink(glue::glue('{current_file_path}/{i}')) #deletes cached image, so it can counts remaining unmatched
  } else {
    #make gene-specific dir
    save_path <- here::here(target_file_path, cell_name)
    if(!dir.exists(save_path)) {dir.create(save_path, recursive = FALSE)}
    
    #move image
    file.rename(from = here::here(current_file_path, i), 
                to = here::here(save_path, i))
    
    #update progress
    message(glue::glue('moved {i}'))
  }
}
#finish message
remaining_files <- list.files(current_file_path, full.names = FALSE, recursive = TRUE)
message("all cell images moved")
message(glue::glue('{length(remaining_files)} files remain'))
