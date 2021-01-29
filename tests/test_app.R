library(here)
source(here::here("tests", "util.R"))

test_that("app_data_dir is set to 'data'", {
  source_code_ignoring_errors("app.R")
  expect_equal(app_data_dir, "data")
})

test_that("private is set to TRUE", {
  source_code_ignoring_errors("app.R")
  expect_equal(private, TRUE)
})
