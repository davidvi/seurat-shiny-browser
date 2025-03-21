# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set the data folder path and port
# Default is "/home/shiny-app/data/" and port 3030
# Can be overridden with command-line arguments

# Initialize defaults
data_folder <- "/home/shiny-app/data/"
port <- 3030

# Parse arguments
if (length(args) > 0) {
  data_folder <- args[1]
  message("Setting data folder to: ", data_folder)
  
  # Check if a second argument for port is provided
  if (length(args) > 1) {
    port_arg <- as.integer(args[2])
    # Validate port is a number and in valid range
    if (!is.na(port_arg) && port_arg > 0 && port_arg < 65536) {
      port <- port_arg
      message("Setting port to: ", port)
    } else {
      message("Invalid port specified. Using default port: ", port)
    }
  } else {
    message("Using default port: ", port)
  }
} else {
  message("Using default data folder: ", data_folder)
  message("Using default port: ", port)
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
options(shiny.port = port)
shinyApp(ui, server)
