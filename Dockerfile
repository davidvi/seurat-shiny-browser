FROM r-base:4.2.2


# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    libglpk40 \
    libglpk-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('preprocessCore', 'Seurat'))"

# Install CRAN packages
RUN R -e "install.packages(c('tidyverse', 'ggplot2', 'DT', 'shinyWidgets', 'shinyjs', 'readxl'))"

RUN R -e "installed_packages <- installed.packages(); print(installed_packages[, 'Package'])"

# Create directories
RUN mkdir -p /home/shiny-app/data

# Create directories for app structure
RUN mkdir -p /home/shiny-app/tabs

# Copy application files
COPY inst/shiny-app/global.R /home/shiny-app/global.R
COPY inst/shiny-app/ui.R /home/shiny-app/ui.R
COPY inst/shiny-app/server.R /home/shiny-app/server.R
COPY inst/shiny-app/app.R /home/shiny-app/app.R
COPY inst/shiny-app/tabs/*.R /home/shiny-app/tabs/
COPY inst/shiny-app/server/*.R /home/shiny-app/server/

# Copy package files for running as package
COPY R/run_app.R /home/shiny-app/run_app.R
COPY R/utils.R /home/shiny-app/utils.R

# Expose port
EXPOSE 3030

# Create directory for server files
RUN mkdir -p /home/shiny-app/server

# Set working directory
WORKDIR /home/shiny-app

# Create a startup script
RUN echo '#!/usr/bin/env Rscript\n\
args <- commandArgs(trailingOnly = TRUE)\n\
data_folder <- args[1]\n\
port <- as.integer(args[2])\n\
if(is.na(port)) port <- 3030\n\
source("run_seurat_browser.R")\n\
run_seurat_browser(data_folder = data_folder, port = port)\n\
' > /home/shiny-app/docker_start.R

# Make it executable
RUN chmod +x /home/shiny-app/docker_start.R

# Set entry point that allows for command-line arguments
ENTRYPOINT ["Rscript", "docker_start.R"]

# Default command - will be overridden if arguments are passed to docker run
CMD ["/home/shiny-app/data", "3030"]
