# Development script for running the app locally
# This script is mainly for development and testing
# 
# Usage: 
#   Rscript dev_run.R [data_folder_path] [port]
#
# Examples:
#   Rscript dev_run.R                            # Use default data folder and port
#   Rscript dev_run.R /path/to/data              # Specify data folder
#   Rscript dev_run.R /path/to/data 8080         # Specify data folder and port

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set the data folder path from command line or use default
if (length(args) >= 1 && !is.na(args[1]) && args[1] != "") {
  data_folder <- args[1]
  message("Using data folder from command line: ", data_folder)
} else {
  data_folder <- "./data"
  message("Using default data folder: ", data_folder)
}

# Set the port from command line or use default
if (length(args) >= 2 && !is.na(args[2]) && as.integer(args[2]) > 0) {
  port <- as.integer(args[2])
  message("Using port from command line: ", port)
} else {
  # Default port - change this if 3030 is already in use on your system
  port <- 3030
  message("Using default port: ", port)
}

# Set options for the app
options(future.globals.maxSize = 10000 * 1024^2) # Set to ~10GB
options(seurat_browser_data_folder = data_folder)

# Load necessary libraries
library(shiny)
library(Seurat)
library(SeuratObject)
library(shinyWidgets)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(future)
library(readxl)

# Try to load SingleR and celldex if they are installed
tryCatch({
  suppressWarnings(library(SingleR))
  message("SingleR loaded successfully")
}, error = function(e) {
  message("SingleR package not found. Auto Name tab functionality will be limited.")
  message("To install SingleR: BiocManager::install('SingleR')")
})

tryCatch({
  suppressWarnings(library(celldex))
  message("celldex loaded successfully")
}, error = function(e) {
  message("celldex package not found. Auto Name tab functionality will be limited.")
  message("To install celldex: BiocManager::install('celldex')")
})

# Load BiocParallel if available
tryCatch({
  suppressWarnings(library(BiocParallel))
  message("BiocParallel loaded successfully")
}, error = function(e) {
  message("BiocParallel package not found. Consider installing it for better performance.")
})

# Define direct path to the Shiny app directory
app_dir <- file.path(getwd(), "inst", "shiny-app")
message("Running app directly from: ", app_dir)

# Run the app in development mode with the specified port
message("Starting Seurat Shiny Browser on port: ", port)
message("You can access the app at: http://localhost:", port)
shiny::runApp(app_dir, port = port, launch.browser = FALSE, host = "0.0.0.0")