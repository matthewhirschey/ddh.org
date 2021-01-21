library(tidyverse)

make_summary <- function(data_values) { #do I need to carry over summary table vars?
  if(data_values$type == "gene"){
    summary_table <- tibble(
      identifier = summary_gene(summary_table = gene_summary, data_values, var = "approved_symbol"),
      name = summary_gene(summary_table = gene_summary, data_values, var = "approved_name"),
      summary = summary_gene(summary_table = gene_summary, data_values, var = "entrez_summary"))
  } else if (data_values$type == "pathway") {
    summary_table <- tibble(
      identifier = paste0("GO:", summary_pathway(summary_table = pathways, data_values, var = "go")),
      name = summary_pathway(summary_table = pathways, data_values, var = "pathway"),
      summary = summary_pathway(summary_table = pathways, data_values, var = "def"))
  } else { #gene_list
    summary_table <- tibble(
      identifier = c("custom gene list"),
      name = summary_gene_list(summary_table = gene_summary, data_values),
      summary = c("User defined gene input list"))
  }
  return(summary_table)
}

#tests
#make_summary(data_values = list(id="GSS", type="gene", gene_symbols=c("GSS")))
#make_summary(data_values = list(id="0060148", type="pathway", gene_symbols=c("SDHA", "SDHB")))
#make_summary(data_values = list(id="SDHA,SDHB", gene_symbols=c("SDHA", "SDHB"), type="gene_list"))

#render in temp dir replaces usual render function
render_rmarkdown_in_tempdir <- function(data_values, rmd_path, output_file, envir = parent.frame()) {
  # The rmd_path variable must be an absolute path.

  # make sure the base report directory exists
  report_base_dir = here::here("report")
  if (!file.exists(report_base_dir)) {
    dir.create(report_base_dir)
  }
  # determine the filename of the Rmd file we will use for rendering
  rmd_filename <- basename(rmd_path)
  # create a temporary directory and make it our working directory
  temp_dir <- tempfile(pattern="tmpdir", tmpdir=report_base_dir)
  dir.create(temp_dir)
  owd <- setwd(temp_dir)
  on.exit(setwd(owd))
  on.exit(unlink(temp_dir, recursive = TRUE))
  # copy the Rmd file into our temporary(current) directory
  file.copy(rmd_path, rmd_filename, overwrite = TRUE)

  #good file names
  good_file_name <- data_values$id
  if (data_values$type == "gene_list") {
    good_file_name <- paste0("custom_", paste0(data_values$gene_symbols, collapse="_"))
  }
  #zip
  output_html_filename <- paste0(good_file_name, "_report.html")
  zip_filenames <- c(output_html_filename)

  # bring in network from parent environment
  network <- get("network", envir = envir)
  network_filename <- paste0(good_file_name, "_", "graph")

  # Save base network with legend added for the final zip
  visSave(network %>% visLegend(position = "right", width = .25, zoom = F), file = paste0(network_filename, "Interactive.html"))
  zip_filenames <- append(zip_filenames, paste0(network_filename, "Interactive.html"))

  # assign filename information to envir so the network can be found within the report_gene.Rmd file
  assign("networkPath", paste0(network_filename, "Interactive.html"), envir = envir)

  rmarkdown::render(rmd_filename, output_file = output_html_filename, envir = envir)
  # get the names of all the items included for rendering
  for (name in names(envir)) {
    env_item = envir[[name]]
    # if the env_item is a plot
    if ("plot_env" %in% names(env_item)) {
      # save the plot to a png
      plot_filename <- paste0(good_file_name, "_", name, ".png")
      # custom heights for lineage plots
      if (name == "lineage") {
        ggsave(plot_filename, width = 12, height = 10.5, env_item, dpi = 300, type = "cairo")
      }
      if (name == "sublineage") {
        ggsave(plot_filename, env_item, width = 12, height = 22, dpi = 300, type = "cairo")
      }
      # dynamic height depending on # of genes for cellbins plot
      if (name == "cellbins") {
        ggsave(plot_filename, env_item, width = 12, dpi = 300, type = "cairo",
               height = length(unique(data_values$gene_symbols)) + 3, limitsize = FALSE)
      }
      # smaller plots for anatogram and network
      if (name == "cellanatogram") {
        ggsave(plot_filename, env_item, width = 8, height = 7, dpi = 300, type = "cairo")
      }
      if (name == "graph") {
        ggsave(plot_filename, env_item, width = 12, height = 10.5, dpi = 300, type = "cairo")
      }
      # landscape aspect ratio for all other plots
      if (!name %in% c("lineage", "sublineage", "cellbins", "cellanatogram")) {
        ggsave(plot_filename, env_item, width = 12, height = 7.5, dpi = 300, type = "cairo")
      }
      # include the plot png in the zip download
      zip_filenames <- append(zip_filenames, plot_filename)
    }
  }
  zip(zipfile = output_file, files = zip_filenames)
}

#specific instructions to render reports based on query type and report template
render_gene_report <- function(data_values, output_file) {
  if(data_values$type == "gene" | data_values$type == "pathway" | data_values$type == "gene_list") {
    gene_symbol <- data_values$gene_symbols
  } else {
    stop("declare your type!")
  }
  num <- length(achilles$X1)
  summary <- make_summary(data_values)
  cellanatogram <- make_cellanatogram(cellanatogram_data = subcell, gene_symbol)
  cellanatogram_table <- make_cellanatogram_table(cellanatogram_data = subcell, gene_symbol)
  celldeps <- make_celldeps(celldeps_data = achilles, expression_data = expression_names, gene_symbol, mean = mean_virtual_achilles)
  cellbins <- make_cellbins(cellbins_data = achilles, expression_data = expression_names, gene_symbol)
  genecorrelation <- make_correlation(table_data = achilles_cor_nest, gene_symbol, mean = mean_virtual_achilles, upper_limit = achilles_upper, lower_limit = achilles_lower)
  lineage <- make_lineage(celldeps_data = achilles, expression_data = expression_names, gene_symbol)
  sublineage <- make_sublineage(celldeps_data = achilles, expression_data = expression_names, gene_symbol)
  expression_plot <- make_cellexpression(expression_data = expression, expression_join = expression_names, gene_symbol, mean = mean_virtual_expression, upper_limit = expression_upper, lower_limit = expression_lower)
  expression_table <- make_expression_table(expression_data = expression, expression_join = expression_names, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles_data = achilles, expression_data = expression_names, gene_symbol) %>% head(10)
  target_achilles_top <- make_achilles_table(achilles_data = achilles, expression_data = expression_names, gene_symbol) %>% tail(10)
  dep_top <- make_top_table(toptable_data = master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_top(enrichmenttop_data = master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(bottomtable_data = master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_bottom(enrichmentbottom_data = master_negative, gene_symbol)
  network <- make_graph(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2, corrType = "Positive and Negative", displayHeight = '80vh', displayWidth = '100%')
  render_rmarkdown_in_tempdir(data_values, here::here("code", "report_gene.Rmd"), output_file)
}

# render_gene_report(data_values = list(id="GSS", type="gene", gene_symbols=c("GSS")), output_file=here::here("/report/dataLarge.zip"))
# #render_gene_report(input = "0060148", type = "pathway", output_file = "0060148.zip")
# #render_gene_report(input = c("GSS", "SST"), type = "gene_list")

#logic to matching query type to rendered content
render_report_to_file <- function(data_values, file) {
  if (data_values$type == "gene" | data_values$type == "pathway" | data_values$type == "gene_list") {
    render_gene_report(data_values, output_file = file)
  } else {
    stop("no report for you!")
  }
}

#render_report_to_file(input = "GSS", type = "gene", file = "gss_trial.pdf")
