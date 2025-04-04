#' Get Seurat Objects from a Folder
#'
#' @description
#' Lists RDS files in a given folder that may contain Seurat objects.
#'
#' @param folder The path to the folder to search in.
#' @param base_folder The base data directory path.
#'
#' @return A character vector of file names.
#' @keywords internal
get_rds_files <- function(folder, base_folder) {
  folder_path <- file.path(base_folder, folder, "/")
  
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

#' Load a Seurat Object
#'
#' @description
#' Loads a Seurat object from an RDS file.
#'
#' @param sample_file The name of the RDS file to load.
#' @param folder The folder containing the file.
#' @param base_folder The base data directory path.
#'
#' @return A Seurat object or NULL if loading fails.
#' @keywords internal
load_sample <- function(sample_file, folder, base_folder) {
  # Check that file name is provided
  if (is.null(sample_file) || sample_file == "" || length(sample_file) == 0) {
    message("ERROR: No sample file provided to load_sample()")
    return(NULL)
  }
  
  # Construct the file path
  file_path <- file.path(base_folder, folder, sample_file)
  
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
      Seurat::DefaultAssay(sample) <- "RNA"
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