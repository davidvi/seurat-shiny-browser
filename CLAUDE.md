# Seurat Shiny Browser Instructions

## Build & Run Commands
- Build Docker image: `./build.sh`
- Run container: `./run.sh`
- Access app: http://localhost:3030
- Data location: Mount local data to `/home/shiny-app/data` in container

## Code Style Guidelines
- Use tidyverse-style R code with proper indentation (2 spaces)
- Organize code sections with comment headers (e.g., `### rename cluster ###`)
- Use reactive values (`reactiveValues`) for shared state
- Apply proper error handling with tryCatch
- Structure UI with fluidPage and sidebarLayout patterns
- Use progress indicators for long-running operations
- Follow Seurat function naming conventions (CamelCase)
- Include meaningful notifications for user feedback
- Place libraries at top with suppressPackageStartupMessages
- Ensure proper cleanup in observer functions