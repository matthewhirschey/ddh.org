FROM rocker/tidyverse:3.6.1
RUN install2.r --repos http://archive.linux.duke.edu/cran/ --deps TRUE here janitor corrr beepr enrichR moderndive pander vroom rentrez feather optparse tidytext widyr doParallel tidyselect dbplyr dplyr haven jsonlite modelr tidyr webchem
RUN Rscript -e 'devtools::install_github("jespermaag/gganatogram")'
