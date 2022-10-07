#GENERATE STATIC HTML PAGES BEFORE LOADING APP

static_page_path <- here::here("code", "www")

#download methods.zip from S3
methods_file <- "methods.zip"
methods_path_full <- glue::glue("{static_page_path}/{methods_file}")

s3 <- paws::s3()
s3$download_file(Bucket = Sys.getenv("AWS_DATA_BUCKET_ID"),
                 Key = methods_file, 
                 Filename = as.character(methods_path_full))
unzip(zipfile = methods_path_full, 
      overwrite = TRUE, 
      exdir = static_page_path)

#render public reports page
rmarkdown::render(input = here::here("code", "public.Rmd"), 
                  output_format = "html_document", 
                  output_dir = static_page_path, 
                  clean = TRUE)

#run an additional render to create the welcome page
rmarkdown::render(input = here::here("code", "welcome.Rmd"), 
                  output_format = "html_document", 
                  output_dir = static_page_path, 
                  clean = TRUE)

#run an additional render to create the secret report page
rmarkdown::render(input = here::here("code", "report.Rmd"), 
                  output_format = "html_document", 
                  output_dir = static_page_path, 
                  clean = TRUE)
