# Load global settings and functions
source("global.R")

# Load UI definition
source("ui.R")

# Load server function
source("server.R")

# Create and return the Shiny app
shinyApp(ui, server)