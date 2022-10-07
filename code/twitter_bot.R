#twitter_bot R script
source(here::here("code", "install_libraries.R"))
#library(rtweet)

#LOAD DATA --------------------------------------------------------------------
source(here::here("code", "current_release.R")) #read current release information
app_data_dir <- "data" #hard code
source(here::here("code", "app_data.R"))

#LOAD FUNS --------------------------------------------------------------------
sd_threshold <- 3 #need for fun_plots.R
mean_virtual_achilles <- 0 #need for fun_plots.R
source(here::here("code", "token.R"), local = TRUE)
source(here::here("code", "fun_plots.R"), local = TRUE)
source(here::here("code", "fun_graphs.R"), local = TRUE)
source(here::here("code", "fun_colors.R"), local = TRUE) #colors
source(here::here("code", "fun_ggplot_theme.R"), local = TRUE) #themes

#FUNS DF --------------------------------------------------------------------
function_df <- 
  tibble(fun_name = c("make_celldeps", #always include this one, 
                      "make_ideogram", 
                      "make_colorful_female", 
                      "make_sequence", #only card
                      "make_radial", 
                      "make_proteinsize", 
                      #"make_structure", #only card, plot already exists as jpeg
                      "make_pubmed", 
                      "make_cellanatogram", 
                      "make_tissue",
                      "make_cellbins", 
                      "make_lineage", 
                      "make_sublineage", 
                      "make_cellexpression", 
                      "make_cellgeneprotein"), 
         fun_title = c("Is my gene essential in any cells?", #always include this one
                       "Where is my gene located?", 
                       "Where is my gene expressed?", 
                       "What can the protein sequence tell me?", #only card
                       "Does my protein have a signature?", 
                       "How big is my protein?", 
                       #"What could my protein look like?", #only card, plot already exists as jpeg
                       "What's known about my gene?", 
                       "Where is my protein located in a cell?", 
                       "Where is my protein expressed?",
                       "How many cells is my gene essential in?", 
                       "Which cell lineages care about my gene?", 
                       "Which cell sublineages care about my gene?", 
                       "Is my gene highly expressed in cancer cells?", 
                       "How does gene vs. protein expression compare?"))

#TWITTER FUNS --------------------------------------------------------------------
make_status <- function(summary_df = gene_summary, 
                        input = list()){
  twitter_summary <- 
    summary_df %>% 
    dplyr::filter(approved_symbol %in% input$query)
  
  adverb_list <- c("an interesting", "an unexpected", "an unusual", "a surprising", "a strong", "a mind-bending")
  adverb <- sample(adverb_list, 1)
  
  text <- glue::glue('{twitter_summary$approved_symbol}:{twitter_summary$approved_name} has {adverb} association with other genes. ')
  text_url <- glue::glue('Check it out at https://www.datadrivenhypothesis.org/?show=gene&query={input$query}')
  
  tweet_body <- stringr::str_c(text, text_url, collapse ="\n")
  return(tweet_body)
}

get_functions <- function(df = function_df, 
                          n_functions = 3){
  #possible functions
  random_functions <- sample(df$fun_name[2:length(df$fun_name)], n_functions, replace = FALSE) #subset make sure to skip 1
  selected_functions <- c(random_functions, "make_celldeps") #because it's always included
  return(selected_functions)
}


#save parameters
twitter_save <- function(tmp_file, plot_id) {
  ggsave(tmp_file, plot_id, width = 16, height = 9, units = "cm", dpi = 150, device = "png")
}

#generate content
make_media <- function(input_query = query,
                       df = function_df,
                       function_name){
  image <- rlang::exec(function_name, 
                       input = input_query, 
                       card = TRUE)
  plot_title <- 
    df %>% 
    dplyr::filter(fun_name %in% function_name) %>% 
    dplyr::pull(fun_title)
  
  image <-
    image +
    labs(title = plot_title)
  
  twitter_plot <- tempfile(fileext = ".png")
  twitter_save(twitter_plot, plot = image)
  return(twitter_plot)
}

#test single function
#make_media(input_query = query, function_name = "make_ideogram")
#make_media(input_query = query, function_name = "make_tissue")
#make_media(input_query = query, function_name = "make_structure")

#GET QUERY GENE --------------------------------------------------------------------
query <- list(type = "gene",
              query = sample(surprise_genes, 1))


selected_functions <- get_functions()
tweet_vec <- 
  purrr::map_chr(.x = selected_functions, ~ make_media(input_query = query, function_name = .x))

#POST TWEET --------------------------------------------------------------------
post_tweet(status = make_status(input = query), 
           token = token,
           media = c(tweet_vec[1], tweet_vec[2], tweet_vec[3], tweet_vec[4])
)
