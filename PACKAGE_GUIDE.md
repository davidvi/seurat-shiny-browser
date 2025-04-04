# Seurat Shiny Browser Package Guide

## Package Structure

The Seurat Shiny Browser has been converted into a proper R package with the following structure:

```
seurat-shiny-browser/
├── DESCRIPTION              # Package metadata and dependencies
├── LICENSE                  # MIT License file
├── NAMESPACE                # Exported functions and imports
├── R/                       # R functions for the package
│   ├── run_app.R            # Main function to launch the app
│   └── utils.R              # Helper functions
├── README.md                # Project overview and documentation
├── R_PACKAGE_INSTALL.md     # Installation instructions
├── PACKAGE_CLEANUP.md       # File structure cleanup instructions
├── inst/                    # Additional package resources
│   └── shiny-app/           # Shiny application files
│       ├── app.R            # Main app entry point
│       ├── global.R         # Global settings and functions
│       ├── server.R         # Server logic
│       ├── ui.R             # UI definition
│       ├── server/          # Server components
│       └── tabs/            # UI tab components
├── man/                     # Documentation files
└── dev_run.R                # Development script
```

> **Note**: You may currently see duplicate files in both top-level directories (`server/`, `tabs/`, etc.) and in the `inst/shiny-app/` directory. This duplication occurred during the package conversion process. See `PACKAGE_CLEANUP.md` for instructions on how to resolve this duplication and maintain a clean package structure.

## Key Components

1. **Package Functions (`R/`)**:
   - `run_app.R`: The main exported function to launch the application
   - `utils.R`: Helper functions for data handling

2. **Shiny App (`inst/shiny-app/`)**:
   - The complete Shiny application code, organized into modular components
   - Core app files: app.R, global.R, server.R, ui.R
   - UI components in tabs/ directory
   - Server logic broken down in server/ directory

3. **Documentation**:
   - `README.md`: Project overview and usage instructions
   - `R_PACKAGE_INSTALL.md`: Detailed installation guide
   - `man/`: Auto-generated function documentation
   - `PACKAGE_GUIDE.md`: This guide explaining the package structure

## Usage Workflow

1. **Installation**:
   ```r
   devtools::install_github("yourusername/seuratShinyBrowser")
   ```

2. **Loading**:
   ```r
   library(seuratShinyBrowser)
   ```

3. **Running**:
   ```r
   run_app(data_folder = "/path/to/your/data")
   ```

## Development Workflow

When developing or modifying the package:

1. Edit the source files in `R/` or `inst/shiny-app/`
2. Run `devtools::document()` to update documentation
3. Use `devtools::load_all()` for testing without reinstalling
4. For quick testing, use `Rscript dev_run.R`
5. Install the updated package with `devtools::install()`

## Adding New Features

When adding new features to the package:

1. For new functions, add them to appropriate R/ files with roxygen comments
2. For UI changes, modify files in inst/shiny-app/tabs/
3. For server logic, modify files in inst/shiny-app/server/
4. Update dependencies in DESCRIPTION if new packages are required
5. Document any new exported functions with roxygen2 comments

## Package Dependencies

The package has the following dependencies:
- Seurat
- shiny
- shinyWidgets
- shinyjs
- ggplot2
- tidyverse
- DT
- future
- readxl
- preprocessCore

These are specified in the DESCRIPTION file and will be installed automatically when using `devtools::install_github()`.

## Deployment Options

The Seurat Shiny Browser package can be deployed in several ways:

1. **Local Installation**: Install the package and run on your local machine
2. **Docker Container**: Use the provided Dockerfile to build a container
3. **Shiny Server**: Deploy on a Shiny Server instance for multi-user access
4. **shinyapps.io**: Deploy to RStudio's shinyapps.io service