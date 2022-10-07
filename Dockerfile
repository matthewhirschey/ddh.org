FROM rocker/shiny-verse:4.2.0

# Install requirements for R libraries
RUN apt-get update && apt-get install -y libglpk-dev libmagick++-dev libavfilter-dev libpoppler-cpp-dev python3-dev python3.8-venv

# Specify directory so reticulate will install in a shared location
# ENV WORKON_HOME=/opt/virtualenvs

# Install and verify libraries
ADD code/install_libraries.R /code/install_libraries.R
RUN Rscript /code/install_libraries.R

# Install python modules for reticulate
# ADD code/install_python_libraries.R /code/install_python_libraries.R
# RUN Rscript /code/install_python_libraries.R

# Fix reticulate segfault - https://github.com/rstudio/reticulate/issues/1133#issuecomment-1021783041
# RUN /opt/virtualenvs/gls_regression/bin/pip3 install --no-binary="numpy" numpy --ignore-installed
