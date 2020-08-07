library(tidyverse)
library(tidygraph)
library(networkD3)
library(ggraph)
library(viridis)

source(here::here("code", "fun_tables.R")) #for make_table funs

setup_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10) {
  #make empty tibble
  dep_network <- tibble()
  #either find top/bottom correlated genes if given single gene, or take list to fill gene_list
  if(length(gene_symbol) == 1){
    #find top and bottom correlations for fav_gene
    dep_top <- make_top_table(toptable_data, gene_symbol) %>%
      slice(1:threshold)
    
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>%
      slice(1:threshold) #limit for visualization?
    
    #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
    gene_list <- dep_top %>%
      bind_rows(dep_bottom) %>%
      dplyr::pull("Gene")
  } else {
    gene_list <- toptable_data %>% #this code ensures that the list of genes from a pathway are in the data
      dplyr::filter(fav_gene %in% gene_symbol) %>%
      dplyr::pull(fav_gene)
  }
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in gene_list){
    message("Getting correlations from ", i)
    dep_top_related <- toptable_data %>%
      dplyr::filter(fav_gene == i) %>%
      tidyr::unnest(data) %>%
      dplyr::select(-fav_gene) %>%
      dplyr::arrange(desc(r2)) %>%
      dplyr::slice(1:threshold) %>% 
      dplyr::mutate(x = i, origin = "pos") %>% 
      dplyr::rename(y = gene) %>%
      dplyr::select(x, y, r2, origin)
    
    dep_bottom_related <- bottomtable_data %>%
      dplyr::filter(fav_gene == i) %>%
      tidyr::unnest(data) %>%
      dplyr::select(-fav_gene) %>%
      dplyr::arrange(r2) %>%
      dplyr::slice(1:threshold) %>% 
      dplyr::mutate(x = i, origin = "neg") %>% 
      dplyr::rename(y = gene) %>%
      dplyr::select(x, y, r2, origin)
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  return(dep_network)
}
  
make_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold)
  
  if(length(gene_symbol) == 1){
    dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
  } else {
    dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  if(length(gene_symbol) == 1){
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_top$Gene ~ "Positive",
                               name %in% dep_bottom$Gene ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
        arrange(group)
    } else {
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                 name %in% dep_network_top ~ "Positive",
                                 name %in% dep_network_bottom ~ "Negative",
                                 TRUE ~ NA),
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
        arrange(group) 
    }
  
  links <- graph_network %>%
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>%
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
    as.data.frame
  
  # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1
  
  #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
  if(nrow(links_filtered) == 0) {links_filtered <- tibble("from" = 0, "to" = 0, "r2" = 1, "origin" = "pos")}
  
  #use color meter to get hexdec color values
  node_color <- 'd3.scaleOrdinal(["#EDA555","#AD677D", "#0C2332", "#544097"])'
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene")
    nodes_filtered <- bind_rows(nodes_filtered, dummy)
    node_color <- 'd3.scaleOrdinal(["#AD677D", "#0C2332", "#544097", "#EDA555"])' #change color order for consistency
  }

  forceNetwork(Links = links_filtered, 
               Nodes = nodes_filtered, 
               Source = "from", 
               Target ="to", 
               NodeID = "name", 
               Group = "group", 
               zoom = TRUE, 
               bounded = TRUE, 
               opacity = 0.8,
               opacityNoHover = 100, 
               Nodesize = "degree", 
               colourScale = node_color, 
               legend = TRUE)
}

make_graph_report <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold)
  
  #make some objs for below
  if(length(gene_symbol) == 1){
    dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
  } else {
    dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  if(length(gene_symbol) == 1){
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_top$Gene ~ "Positive",
                               name %in% dep_bottom$Gene ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
      arrange(group)
  } else {
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_network_top ~ "Positive",
                               name %in% dep_network_bottom ~ "Negative",
                               TRUE ~ NA),
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
      arrange(group) 
  }
  
  links <- graph_network %>%
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>%
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
    as.data.frame
  
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id)
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id)
  
  colors <- c("#EDA555","#AD677D", "#0C2332", "#544097")
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene")
    nodes_filtered <- bind_rows(nodes_filtered, dummy)
    #colors <- c("#AD677D", "#0C2332", "#544097", "#EDA555") #no need to reset colors; do that in 'breaks' below
    }
  
  graph_network_ggraph <- tidygraph::tbl_graph(nodes = nodes_filtered, edges = links_filtered)
  
  graph_network_ggraph %>%
    ggraph::ggraph(layout = "auto") +
    geom_edge_fan() + #edge_width = aes(abs(r2)), alpha = 0.3
    geom_node_point(aes(size = degree, color = group), alpha = 0.8) +
    geom_node_label(aes(filter = group != "Connected", label = name), repel = TRUE) +
    scale_colour_manual(values = colors, breaks = c("Query Gene", "Positive", "Negative", "Connected")) +
    theme_graph(base_family = 'Helvetica') +
    guides(size = "none", color = guide_legend(""))
}

#figure legend
graph_title <- "Network Graph."
graph_legend <- "Each point represents a single gene taken from the top associated genes with the query gene. Genes with only one connection were removed."
graph_legend_list <- "Each point represents one of the queried genes, and then the top and bottom associated genes with it. Genes with only one connection were removed."
