# Development script for running the app locally
# This script is mainly for development and testing

# Set the data folder path - change this to your data location
data_folder <- "./data"

# Set the port - change this if 3030 is already in use on your system
# You can use any valid port number between 1024-65535 that's not already in use
port <- 3030

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