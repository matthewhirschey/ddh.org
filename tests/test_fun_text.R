library(here)
library(tidyverse)

# create dummy gene_summary tibble
approved_symbol <- c("RDX","ROCK2")
approved_name = c("Radixin", "Rho associated kinase 2")
gene_summary <- tibble(approved_symbol=approved_symbol, approved_name=approved_name)

# create dummy pathways tibble
go <- c("123456")
pathway <- c("My Protein")
gene <- c("RDX", "ROCK2")
pathways <- tibble(go=go, pathway=pathway, gene=gene) %>%
  group_by(go, pathway) %>%
  nest(data = gene)

source(here::here("code/fun_text.R"), local=TRUE)

test_that("summary_gene pulls approved_symbol from gene_summary by default", {
  data_params <- list(type="gene", id="RDX", gene_symbols=c("RDX"))
  result <- summary_gene(input=data_params)
  expect_equal(result, "RDX")
})

test_that("summary_gene pulls approved_name from gene_summary when specified", {
  data_params <- list(type="gene", id="RDX", gene_symbols=c("RDX"))
  result <- summary_gene(input=data_params, var="approved_name")
  expect_equal(result, "Radixin")
})

test_that("summary_pathway pulls pathway from pathways by default", {
  data_params <- list(type="pathway", id="123456", gene_symbols=c("RDX"))
  result <- summary_pathway(input=data_params)
  expect_equal(result, "My Protein")
})

test_that("summary_pathway pulls data seperated by commas from pathways when specified", {
  data_params <- list(type="pathway", id="123456", gene_symbols=c("RDX"))
  result <- summary_pathway(input=data_params, var="data")
  expect_equal(result, "RDX, ROCK2")
})

test_that("summary_gene_list filters out invalid symbols from gene_summary", {
  data_params <- list(type="gene_list", id="RDX,ROCK2,OTHER", gene_symbols=c("RDX","ROCK2","OTHER"))
  result <- summary_gene_list(input=data_params)
  expect_equal(result, "RDX, ROCK2")
})

