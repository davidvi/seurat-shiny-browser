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
