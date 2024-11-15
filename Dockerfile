FROM r-base:latest

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

# Install BiocManager and R packages
RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install(c('preprocessCore', 'Seurat'))" \
    && R -e "install.packages(c('shinyWidgets', 'shinyjs', 'readxl', 'limma', 'ggplot2', 'tidyverse', 'DT'))"

# Create directories
RUN mkdir -p /home/shiny-app/data

# Copy server script
COPY server.R /home/shiny-app/server.R

# Expose port
EXPOSE 3030

# Run the server script
CMD ["Rscript", "/home/shiny-app/server.R"]
