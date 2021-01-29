library(tidyverse)
library(tidygraph)
library(ggraph)
library(viridis)
library(visNetwork)

source(here::here("code", "fun_tables.R")) #for make_table funs

setup_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, corrType = "Positive and Negative") {
  #make empty tibble
  dep_network <- tibble()
  
  # make the correct graph including only correlations of the designated type
  if(corrType == "Positive and Negative"){
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
        dplyr::ungroup(.) %>% 
        dplyr::arrange(desc(r2)) %>%
        dplyr::slice(1:threshold) %>% 
        dplyr::mutate(x = i, origin = "pos") %>% 
        dplyr::rename(y = gene) %>%
        dplyr::select(x, y, r2, origin)
      
      dep_bottom_related <- bottomtable_data %>%
        dplyr::filter(fav_gene == i) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup(.) %>% 
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
  else if(corrType == "Positive"){
    #either find top correlated genes if given single gene, or take list to fill gene_list
    if(length(gene_symbol) == 1){
      #find top correlations for fav_gene
      dep_top <- make_top_table(toptable_data, gene_symbol) %>%
        slice(1:threshold)
      
      #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
      gene_list <- dep_top %>%
        dplyr::pull("Gene")
    } else {
      gene_list <- toptable_data %>% #this code ensures that the list of genes from a pathway are in the data
        dplyr::filter(fav_gene %in% gene_symbol) %>%
        dplyr::pull(fav_gene)
    }
    #this loop will take each gene, and get their top correlations, and build a df containing the top n number of genes for each gene
    for (i in gene_list){
      message("Getting correlations from ", i)
      dep_top_related <- toptable_data %>%
        dplyr::filter(fav_gene == i) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup(.) %>% 
        dplyr::arrange(desc(r2)) %>%
        dplyr::slice(1:threshold) %>% 
        dplyr::mutate(x = i, origin = "pos") %>% 
        dplyr::rename(y = gene) %>%
        dplyr::select(x, y, r2, origin)
      
      #Bind to the final df for graphing
      dep_network <- dep_network %>%
        bind_rows(dep_top_related)
    }
    return(dep_network)
  } 
  else if(corrType == "Negative"){  
    #either find bottom correlated genes if given single gene, or take list to fill gene_list
    if(length(gene_symbol) == 1){
      #find bottom correlations for fav_gene
      dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>%
        slice(1:threshold) #limit for visualization?
      
      #this takes the genes from the  bottom, and pulls them to feed them into a for loop
      gene_list <- dep_bottom %>%
        dplyr::pull("Gene")
    } else {
      gene_list <- bottomtable_data %>% #this code ensures that the list of genes from a pathway are in the data
        dplyr::filter(fav_gene %in% gene_symbol) %>%
        dplyr::pull(fav_gene)
    }
    #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
    for (i in gene_list){
      message("Getting correlations from ", i)
      
      dep_bottom_related <- bottomtable_data %>%
        dplyr::filter(fav_gene == i) %>%
        tidyr::unnest(data) %>%
        dplyr::ungroup(.) %>% 
        dplyr::arrange(r2) %>%
        dplyr::slice(1:threshold) %>% 
        dplyr::mutate(x = i, origin = "neg") %>% 
        dplyr::rename(y = gene) %>%
        dplyr::select(x, y, r2, origin)
      
      #Bind to the final df for graphing
      dep_network <- dep_network %>%
        bind_rows(dep_bottom_related)
    }
    return(dep_network)}
}
#tests
#setup_graph(gene_symbol = "SDHA")
#setup_graph(gene_symbol = c("SDHA", "SDHB"))
#setup_graph(gene_symbol = c("GSS", "SST"))


#' Create network graph visualization using visNetwork
#' 
#' This function takes in dependency correlations and a gene query list to then output a dependency network graph 
#' visualization containing the top/bottom threshold for each of the top/bottom threshold of the gene query list 
#' using visNetwork.
#'
#' @param toptable_data A tibble of genes and their associated top correlated genes
#' @param bottomtable_data A tibble of genes and their associated bottom correlated genes
#' @param gene_symbol A character or character vector of gene_symbols used to create network graph
#' @param threshold A numerical representing the number of genes to pull from top and bottom tables
#' @param deg A numerical representing the minimum number of connections for a gene to be connected to the network
#' @param corrType A string that describes what type of correlations to include, options are: "Positive and Negative", "Positive", or "Negative"
#' @param displayHeight Default to "90vh". The height of the network in pixels("500px"), as a percentage("100%"), or as a percentage of the viewport("70vh", where 70 represents 70% of the viewport)
#' @param displayWidth Default to "100%". The width of the network in pixels("500px"), as a percentage("100%"), or as a percentage of the viewport("70vh", where 70 represents 70% of the viewport)
#' @param tooltipLink Boolean to denote whether or not to include a link in the tooltip for a gene. Default to false.
#'
#' @return NULL - Outputs a complete network graph
#' @export
#'
#' @examples
make_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2, corrType = "Positive and Negative", displayHeight = '90vh', displayWidth = '100%', tooltipLink = FALSE) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold, corrType)
  
  if(corrType == "Positive and Negative"){
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
                                 name %in% dep_top$Gene == TRUE ~ "Positive",
                                 name %in% dep_bottom$Gene == TRUE ~ "Negative",
                                 TRUE ~ "Connected"), 
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
        arrange(group)
    } else {
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                        name %in% dep_network_top == TRUE ~ "Positive",
                                        name %in% dep_network_bottom == TRUE ~ "Negative"),
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
        arrange(group) 
    }
  }  else if(corrType == "Positive"){
    if(length(gene_symbol) == 1){
      dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    } else {
      dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    }
    
    #make graph
    graph_network <- tidygraph::as_tbl_graph(dep_network)
    if(length(gene_symbol) == 1){
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                 name %in% dep_top$Gene == TRUE ~ "Positive",
                                 TRUE ~ "Connected"), 
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Positive", "Connected")))  %>%
        arrange(group)
    } else {
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                        name %in% dep_network_top == TRUE ~ "Positive"),
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Positive")))  %>% #you don't end up with "connected" in a multi-gene list
        arrange(group) 
    }
  }  else if (corrType == "Negative"){
    if(length(gene_symbol) == 1){
      dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
    } else {
      dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
    }
    
    #make graph
    graph_network <- tidygraph::as_tbl_graph(dep_network)
    if(length(gene_symbol) == 1){
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                 name %in% dep_bottom$Gene == TRUE ~ "Negative",
                                 TRUE ~ "Connected"), 
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Negative", "Connected")))  %>%
        arrange(group)
    } else {
      nodes <-  as_tibble(graph_network) %>%
        rowid_to_column("id") %>%
        mutate(degree = igraph::degree(graph_network),
               group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                        name %in% dep_network_bottom == TRUE ~ "Negative"),
               group = as_factor(group), 
               group = fct_relevel(group, c("Query Gene", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
        arrange(group) 
    }
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
  if(nrow(links_filtered) == 0) {
    if(corrType == "Negative"){
      links_filtered <- tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "neg")  
    }    else{links_filtered <- tibble("from" = -1, "to" = -1, "r2" = 1, "origin" = "pos")}
  }

  # shift id values properly to work with visNetwork
  nodes_filtered <- nodes_filtered %>%
    dplyr::mutate(id=0:(dim(nodes_filtered)[1]-1))
  
  # recreate the network using filtered edges to use degree function
  graph_network_filtered <- igraph::graph_from_data_frame(links_filtered, directed = F) %>% 
    igraph::simplify()
  degVec <- igraph::degree(graph_network_filtered)
  
  # make node size a function of their degree in the current network
  if(names(degVec)[1] == "-1"){ # if there are no links remaining then assign degree of each node to 0
    nodes_filtered <- nodes_filtered %>%  
      mutate(value = 0)
    }else if(length(gene_symbol > 1) & length(degVec) != dim(nodes_filtered)[1]){ # handle cases where some query genes are disconnected
      genesWithConnections <- degVec %>% 
        names() %>% 
        as.numeric()
      nodes_filtered <- nodes_filtered %>% 
        add_column(value = 0)
      for(gene in 1:dim(nodes_filtered)[1]){
        if(nodes_filtered[gene,"id"] %in% genesWithConnections){
          nodes_filtered[gene,"value"] <- degVec[nodes_filtered[gene,"id"] %>% toString()]
        }
      }
    }else{
    nodes_filtered <- nodes_filtered[degVec %>%  names() %>% as.numeric() +1,] %>%  
      mutate(value = degVec)
    }
  
  # make sure query gene is at the start so legend shows proper order
  nodes_filtered <- nodes_filtered %>% 
    arrange(id)
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  disconnected <- F
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene", "value" = 0)
    nodes_filtered <- bind_rows(dummy, nodes_filtered)
    disconnected <- T
  }
  
  # get the approved_name of each gene from the gene_summary table - will be added to nodes tibble for tooltip
  nameTable <- tibble(name=character())
  for(gene in nodes_filtered$name){
    newVal <- gene_summary %>% 
      dplyr::filter(approved_symbol==gene) %>% 
      dplyr::pull(approved_name)
    if(length(newVal)==0){
      nameTable <- add_row(nameTable, name = "No gene summary information available")# handles cases where the gene is not in the gene summary table
    } else{
      nameTable <- add_row(nameTable, name=newVal)
    }
  }
  
  # add title information (tooltip that appears on hover)  
  if(!tooltipLink){ # Do not form a url when just making the standalone graph while testing or making reports
    nodes_filtered <- nodes_filtered %>%
      dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name, '</p>'),
                    label = nodes_filtered$name )
  }else{
    nodes_filtered <- nodes_filtered %>%
      dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name ,'<br><a target="_blank" href="?show=gene&query_type=gene&symbol=',nodes_filtered$name,'">Gene Link</a></p>'),
                    label = nodes_filtered$name )
  }
  
  
  # colors used within the network
  queryGeneColor <-"rgba(237, 165, 85, 0.8)"
  positiveColor <-"rgba(173, 103, 125, 0.8)"
  negativeColor <- "rgba(12, 35, 50, 0.8)"
  connectedColor <- "rgba(84, 64, 151, 0.8)"
  borderColor <- "rgba(204, 204, 204, 0.8)" #(gray80), formerly white 255, 255, 255 
  edgeColor <- "rgba(84, 84, 84, 1)"

  # Physics parameters defining the visNetwork
  iter <- 150 # number of iterations to perform of stablization before display
  gravity <- 0.5
  damping <- 0.11
  timestep <- 0.25 # reducing the timestep reduces the jitteriness of the graph and can help stabilize it
  if(disconnected){ # set up function to be called at the end of stabilization, changes the zoom to handle when query gene is disconnected and flies out of the network
    stabilizationZoomFn <- "function() {this.moveTo({scale:0.35})}"
  } else{
    stabilizationZoomFn <- "function() {}"
  }
  
  # build the network visualization
  if(corrType == "Positive and Negative"){
    visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>% 
      visOptions(highlightNearest = list(enabled = T)) %>% 
      visGroups(groupname = "Query Gene", color = list(background = queryGeneColor, border =borderColor, highlight = queryGeneColor, hover = queryGeneColor ), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
      visEdges(color = edgeColor, smooth = F) %>% 
      visNodes(scaling = list(min = 10, max =20)) %>% 
      visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>% 
      visEvents(stabilizationIterationsDone = stabilizationZoomFn)
    # visConfigure(enabled=TRUE) # use to test out new features to add from visNetwork
  }  else if(corrType == "Positive"){
    visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>% 
      visOptions(highlightNearest = list(enabled = T)) %>% 
      visGroups(groupname = "Query Gene", color = list(background = queryGeneColor, border =borderColor, highlight = queryGeneColor, hover = queryGeneColor ), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Positive", color = list(background = positiveColor, border = borderColor, highlight = positiveColor, hover = positiveColor), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
      visEdges(color = edgeColor, smooth = F) %>% 
      visNodes(scaling = list(min = 10, max =20)) %>% 
      visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>% 
      visEvents(stabilizationIterationsDone = stabilizationZoomFn) 
  }  else if(corrType == "Negative"){
    visNetwork(nodes = nodes_filtered, edges = links_filtered, width = displayWidth, height = displayHeight) %>% 
      visOptions(highlightNearest = list(enabled = T)) %>% 
      visGroups(groupname = "Query Gene", color = list(background = queryGeneColor, border =borderColor, highlight = queryGeneColor, hover = queryGeneColor ), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Negative", color = list(background = negativeColor, border = borderColor, highlight = negativeColor, hover = negativeColor), shape='dot', borderWidth = 2) %>%
      visGroups(groupname = "Connected", color = list(background = connectedColor, border = borderColor, highlight = connectedColor, hover = connectedColor), shape='dot', borderWidth = 2) %>%
      visEdges(color = edgeColor, smooth = F) %>% 
      visNodes(scaling = list(min = 10, max = 20)) %>% 
      visPhysics(barnesHut = list(damping = damping, centralGravity = gravity), timestep = timestep, stabilization = list(iterations = iter)) %>% 
      visEvents(stabilizationIterationsDone = stabilizationZoomFn) 
  }
}  

# Test Cases
# make_graph(gene_symbol = "SDHA", threshold = 10)
# make_graph(gene_symbol = "GLI2") # disconnected query gene

# make_graph(gene_symbol = "SDHA", corrType = "Positive")
# make_graph(gene_symbol = "SDHA", corrType = "Negative")
# make_graph(gene_symbol = "CS",threshold = 18, deg = 2)
# make_graph(gene_symbol = c("GSS", "SST"))

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
                               name %in% dep_top$Gene == TRUE ~ "Positive",
                               name %in% dep_bottom$Gene == TRUE ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
      arrange(group)
  } else {
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                      name %in% dep_network_top == TRUE ~ "Positive",
                                      name %in% dep_network_bottom == TRUE ~ "Negative"),
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
  
  #readjust
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id)
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id)
  
  #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
  if(nrow(links_filtered) == 0) {links_filtered <- tibble("from" = 1, "to" = 1, "r2" = 1, "origin" = "pos")}
  
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
    geom_edge_fan(edge_colour = "grey45") + #edge_width = aes(abs(r2)), alpha = 0.3
    geom_node_point(aes(size = degree, color = group), alpha = 0.8) +
    geom_node_label(aes(filter = group != "Connected", label = name), repel = TRUE, family = "Roboto Slab", size = 3) +
    scale_colour_manual(values = colors, breaks = c("Query Gene", "Positive", "Negative", "Connected"), name = NULL) +
    scale_size(range = c(2, 8)) +
    theme_graph(base_family = 'Nunito Sans', base_size = 14) +
    guides(size = "none", color = guide_legend(override.aes = list(size = 3)))
}

#figure legend
graph_title <- "Network Graph."
graph_legend <- "Each point represents a single gene taken from the top associated genes with the query gene. Genes with only one connection were removed."
graph_legend_list <- "Each point represents one of the queried genes, and then the top and bottom associated genes with it. Genes with only one connection were removed."

#make_graph_report(gene_symbol = "TP53")
#make_graph_report(gene_symbol = c("SDHA", "SDHB"))
#make_graph_report(gene_symbol = c("GSS", "SST"))