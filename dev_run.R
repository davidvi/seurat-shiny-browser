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

# Source the run_app function
source("R/run_app.R")

# Set options for the app
options(future.globals.maxSize = 10000 * 1024^2) # Set to ~10GB
options(seurat_browser_data_folder = data_folder)

# Load necessary libraries
library(shiny)
library(Seurat)
library(shinyWidgets)
library(shinyjs)
library(ggplot2)
library(tidyverse)
library(DT)
library(future)
library(readxl)

# Run the app in development mode with the specified port
message("Starting Seurat Shiny Browser on port: ", port)
message("You can access the app at: http://localhost:", port)
run_app(data_folder = data_folder, port = port, launch.browser = TRUE)