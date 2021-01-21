library(here)

test_that("data_dir is set to 'data'", {
  source(here::here("code/current_release.R"))
  expect_equal(app_data_dir, "data")
})
