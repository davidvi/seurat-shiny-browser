suppressPackageStartupMessages({
  library(Seurat)
  library(shiny)
  library(shinyWidgets)
  library(shinyjs)
  library(readxl)
  library(preprocessCore)
  library(ggplot2)
  library(tidyverse)
  library(DT)
  library(future)
})

# load data
# base_folder <- "/home/shiny-app/data/"
base_folder <- "/Users/vanijzen/Developer/seurat-shiny-browser/data/"

# Get folders in data directory
all_folders <- list.dirs(path = base_folder, full.names = FALSE, recursive = FALSE)

# Set the initial folder and get its files
current_folder <- all_folders[1]

# Function to get RDS files from a specific folder
get_rds_files <- function(folder) {
  folder_path <- paste0(base_folder, folder, "/")
  
  # Check if directory exists
  if (!dir.exists(folder_path)) {
    message("Directory does not exist: ", folder_path)
    return(character(0))
  }
  
  # Get the RDS files
  files <- list.files(path = folder_path, pattern = "\\.rds$", full.names = FALSE)
  message("Folder ", folder, " contains ", length(files), " RDS files")
  return(files)
}

# Get initial files from the first folder
all_rds_files <- get_rds_files(current_folder)

sample <- NULL
sample_name <- NULL
cluster_names <- NULL

load_sample <- function(sample_file, folder) {
  # Check that file name is provided
  if (is.null(sample_file) || sample_file == "" || length(sample_file) == 0) {
    message("ERROR: No sample file provided to load_sample()")
    return(NULL)
  }
  
  # Construct the file path
  file_path <- paste0(base_folder, folder, "/", sample_file)
  
  message("DEBUG: Attempting to load file from: ", file_path)
  
  # Check that file exists
  if (!file.exists(file_path)) {
    message("ERROR: File does not exist: ", file_path)
    return(NULL)
  }
  
  # Load the sample with error handling
  tryCatch({
    sample <- readRDS(file = file_path)
    message("DEBUG: Successfully loaded RDS file")
    
    # Set the default assay to RNA if it exists
    tryCatch({
      DefaultAssay(sample) <- "RNA"
      message("DEBUG: Set default assay to RNA")
    }, error = function(e) {
      message("DEBUG: Could not set default assay to RNA: ", e$message)
    })
    
    return(sample)
  }, error = function(e) {
    message("ERROR: Failed to load RDS file: ", e$message)
    return(NULL)
  })
}


# Initialize with empty/default values in case there are no files
sample <- NULL
sample_name <- "No sample loaded"
cluster_names <- c()

# Only try to load a sample if RDS files exist
if (length(all_rds_files) > 0) {
  tryCatch({
    sample <- load_sample(all_rds_files[1])
    sample_name <- all_rds_files[1]
    cluster_names <- levels(Idents(sample))
  }, error = function(e) {
    message("Error loading initial sample: ", e$message)
  })
}

# Get available reductions
available_reductions <- c("umap") # Default

if (!is.null(sample)) {
  tryCatch({
    reductions_list <- Reductions(sample)
    if (length(reductions_list) > 0) {
      available_reductions <- reductions_list
    }
  }, error = function(e) {
    message("Error getting reductions: ", e$message)
  })
}
