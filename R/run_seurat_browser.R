#' Run Seurat Shiny Browser App
#'
#' @description
#' This function launches the Seurat Shiny Browser application.
#' 
#' @param data_folder Path to the folder containing Seurat data files. If NULL, a default path will be used.
#' @param port The port for the Shiny app to run on. Default is 3030. Change this if the default port is already in use.
#' @param launch.browser If TRUE, the system's default web browser will be launched automatically.
#' @param ... Additional arguments to pass to shiny::runApp.
#'
#' @return A Shiny app object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run with default settings
#' run_seurat_browser(data_folder = "~/my_seurat_data/")
#' 
#' # Run on a different port
#' run_seurat_browser(data_folder = "~/my_seurat_data/", port = 8080)
#' }
run_seurat_browser <- function(data_folder = NULL, port = 3030, launch.browser = TRUE, ...) {
  # Increase max global size
  options(future.globals.maxSize = 10000 * 1024^2) # Set to ~10GB
  
  # Set the data folder path
  if (is.null(data_folder)) {
    # Default path
    if (dir.exists("/home/shiny-app/data/")) {
      data_folder <- "/home/shiny-app/data/"
    } else {
      # Try to use current working directory
      data_folder <- getwd()
    }
  }
  
  # Store the data folder path for global.R to use
  options(seurat_browser_data_folder = data_folder)
  message("Setting data folder to: ", data_folder)
  
  # Set Shiny options
  options(shiny.autoreload = TRUE)
  options(shiny.port = port)
  
  # Get the path to the app directory
  app_dir <- system.file("shiny-app", package = "seuratShinyBrowser")
  
  # Check if the app directory exists
  if (app_dir == "") {
    # If running in development mode, check where the app files are located
    if (dir.exists(file.path(getwd(), "inst", "shiny-app"))) {
      # For development with package structure
      app_dir <- file.path(getwd(), "inst", "shiny-app")
      message("Running app from inst/shiny-app directory")
    } else {
      # For development with files at top level
      app_dir <- file.path(getwd())
      message("Running app from top-level directory")
    }
  }
  
  # Launch the Shiny app
  shiny::runApp(app_dir, 
                launch.browser = launch.browser, 
                port = port, 
                host = "0.0.0.0", 
                ...)
}

#' @export
#' @rdname run_seurat_browser
run_app <- function(data_folder = NULL, port = 3030, launch.browser = TRUE, ...) {
  .Deprecated("run_seurat_browser")
  run_seurat_browser(data_folder = data_folder, port = port, launch.browser = launch.browser, ...)
}