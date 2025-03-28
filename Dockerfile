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
COPY global.R /home/shiny-app/global.R
COPY ui.R /home/shiny-app/ui.R
COPY server.R /home/shiny-app/server.R
COPY app.R /home/shiny-app/app.R
COPY tabs/*.R /home/shiny-app/tabs/
COPY server/*.R /home/shiny-app/server/

# Expose port
EXPOSE 3030

# Create directory for server files
RUN mkdir -p /home/shiny-app/server

# Set working directory
WORKDIR /home/shiny-app

# Set entry point that allows for command-line arguments
ENTRYPOINT ["Rscript", "app.R"]

# Default command - will be overridden if arguments are passed to docker run
CMD ["/home/shiny-app/data"]
