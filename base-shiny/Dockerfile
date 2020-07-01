FROM rocker/shiny-verse:3.6.3
LABEL maintainer="Dan Leehr <dan.leehr@duke.edu>"

# Install texlive-xetex for PDF generation
RUN apt-get update && apt-get install -y texlive-xetex

RUN install2.r --repos http://archive.linux.duke.edu/cran/ --deps TRUE \
 ggplot2 \
 tidyverse \
 here \
 corrr  \
 lubridate \
 plotly \
 networkD3 \
 feather \
 tidygraph \
 cowplot \
 digest \
 patchwork \
 viridis \
 pander \
 DT \
 ggraph \
 future \
 shinyWidgets \
 promises

RUN Rscript -e 'devtools::install_github("jespermaag/gganatogram")'
Run Rscript -e 'install.packages(c("dbplyr", "dplyr", "haven", "jsonlite", "modelr", "tidyr"))'

COPY ./shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN chown -R shiny /var/lib/shiny-server/

# OpenShift gives a random uid for the user and some programs try to find a username from the /etc/passwd.
# Let user to fix it, but obviously this shouldn't be run outside OpenShift
RUN chmod ug+rw /etc/passwd
COPY ./fix-username.sh /fix-username.sh
COPY ./shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod a+rx /usr/bin/shiny-server.sh

# Make sure the directory for individual app logs exists and is usable
RUN chmod -R a+rwX /var/log/shiny-server
RUN chmod -R a+rwX /var/lib/shiny-server

