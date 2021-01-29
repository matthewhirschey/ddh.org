FROM rocker/shiny-verse:3.6.3

RUN install2.r --deps TRUE here testthat here corrr enrichR optparse
ADD entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

