# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set the data folder path
# Default is "/home/shiny-app/data/"
# Can be overridden with the first command-line argument
if (length(args) > 0) {
  data_folder <- args[1]
  message("Setting data folder to: ", data_folder)
} else {
  data_folder <- "/home/shiny-app/data/"
  message("Using default data folder: ", data_folder)
}

# Store the data folder path for global.R to use
options(seurat_browser_data_folder = data_folder)

# Load global settings and functions
source("global.R")

# Load UI definition
source("ui.R")

# Load server function
source("server.R")

### run app ###

options(shiny.autoreload = TRUE)
options(shiny.host = "0.0.0.0")
options(shiny.port = 3030)
shinyApp(ui, server)
