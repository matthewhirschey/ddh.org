FROM rocker/tidyverse:3.6.1
RUN install2.r --deps TRUE here janitor corrr beepr enrichR moderndive pander vroom rentrez feather optparse tidytext widyr doParallel
