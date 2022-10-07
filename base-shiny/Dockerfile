FROM rocker/shiny-verse:4.2.0
LABEL maintainer="John Bradley <john.bradley@duke.edu>"

# Install requirements for R libraries
RUN apt-get update && apt-get install -y libglpk-dev libmagick++-dev libavfilter-dev libpoppler-cpp-dev python3-dev python3.8-venv libjpeg-dev tcl8.6-dev

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

# Load methods and pre-render pages
RUN R -e 'source("render_pages.R")'

CMD /usr/bin/shiny-server.sh
