library(tidyverse)
library(tidygraph)
library(networkD3)
library(ggraph)
library(viridis)

make_graph <- function(gene_symbol, threshold = 10, deg = 2) {
  dep_network <- tibble()
  
  #find top and bottom correlations for fav_gene
  dep_top <- make_top_table(gene_symbol) %>%
    slice(1:threshold)
  
  dep_bottom <- make_bottom_table(gene_symbol) %>%
    slice(1:threshold) #limit for visualization?
  
  #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
  related_genes <- dep_top %>%
    bind_rows(dep_bottom) %>%
    dplyr::pull("Gene")
  
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in related_genes){
    message("Getting correlations from ", i, " related to ", gene_symbol)
    dep_top_related <- make_top_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    dep_bottom_related <- make_bottom_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  nodes <-  as_tibble(graph_network) %>%
    rowid_to_column("id") %>%
    mutate(degree = igraph::degree(graph_network),
           group = case_when(name %in% gene_symbol == TRUE ~ "query",
                             name %in% dep_top$Gene ~ "pos",
                             name %in% dep_bottom$Gene ~ "neg",
                             TRUE ~ "connected")) %>%
    arrange(desc(degree))
  
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
  
  #Pos #74D055, Neg #3A568C, Q #FDE825, C #450D53
  node_color <- 'd3.scaleOrdinal(["#74D055", "#3A568C", "#FDE825", "#450D53"])'
  
  forceNetwork(Links = links_filtered, Nodes = nodes_filtered, Source = "from", Target ="to", NodeID = "name", Group = "group", zoom = TRUE, bounded = TRUE, opacityNoHover = 100, Nodesize = "degree", colourScale = node_color)
}

make_graph_report <- function(gene_symbol, threshold = 10, deg = 2) {
  dep_network <- tibble()
  
  #find top and bottom correlations for fav_gene
  dep_top <- make_top_table(gene_symbol) %>%
    slice(1:threshold)
  
  dep_bottom <- make_bottom_table(gene_symbol) %>%
    slice(1:threshold) #limit for visualization?
  
  #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
  related_genes <- dep_top %>%
    bind_rows(dep_bottom) %>%
    dplyr::pull("Gene")
  
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in related_genes){
    message("Getting correlations from ", i, " related to ", gene_symbol)
    dep_top_related <- make_top_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    dep_bottom_related <- make_bottom_table(i) %>% 
      slice(1:threshold) %>% 
      mutate(x = i, origin = "pos") %>% 
      rename(y = Gene, r2 = `R^2`) %>%
      select(x, y, r2, origin)
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  nodes <-  as_tibble(graph_network) %>%
    rowid_to_column("id") %>%
    mutate(degree = igraph::degree(graph_network),
           group = case_when(name %in% gene_symbol == TRUE ~ "query",
                             name %in% dep_top$Gene ~ "pos",
                             name %in% dep_bottom$Gene ~ "neg",
                             TRUE ~ "connected")) %>%
    arrange(desc(degree))
  
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
  
  graph_network_ggraph <- tidygraph::tbl_graph(nodes = nodes_filtered, edges = links_filtered)
  
  graph_network_ggraph %>%
    ggraph::ggraph(layout = "auto") +
    geom_edge_fan(aes(edge_width = abs(r2)), alpha = 0.3) +
    geom_node_point(aes(size = degree, color = group), alpha = 0.7) +
    geom_node_label(aes(label = name), repel = TRUE) +
    scale_colour_viridis(discrete = TRUE, name = "Group", labels = c("Query", "Positive", "Negative", "Connected")) +
    theme_graph(base_family = 'Helvetica')
}
