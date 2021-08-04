FROM rocker/shiny-verse:3.6.3
LABEL maintainer="Dan Leehr <dan.leehr@duke.edu>"

# Install texlive-xetex for PDF generation
RUN apt-get update && apt-get install -y texlive-xetex

RUN install2.r --repos http://archive.linux.duke.edu/cran/ --deps TRUE \
 ggplot2 \
 tidyselect \
 here \
 corrr  \
 lubridate \
 plotly \
 visNetwork \
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
 promises \
 shiny \
 showtext \
 dbplyr \
 dplyr \
 haven \
 jsonlite \
 modelr \
 tidyr \
 ggdist \
 htmltools \
 gt

# install rtweet for twitter bot
RUN install2.r --repos http://archive.linux.duke.edu/cran/ --deps TRUE rtweet

RUN Rscript -e 'devtools::install_github("jespermaag/gganatogram")'

# Install and register "Roboto Slab" and "Nunito Sans" fonts
RUN mkdir /usr/share/fonts/googlefonts \
    && wget --output-document=/tmp/RobotoSlab.zip https://fonts.google.com/download?family=Roboto%20Slab \
    && unzip -d /usr/share/fonts/googlefonts /tmp/RobotoSlab.zip \
    && wget --output-document=/tmp/NunitoSans.zip https://fonts.google.com/download?family=Nunito%20Sans \
    && unzip -d /usr/share/fonts/googlefonts /tmp/NunitoSans.zip \
    && fc-cache -fv

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
