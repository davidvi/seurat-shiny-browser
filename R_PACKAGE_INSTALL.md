# Building and Installing the Seurat Shiny Browser Package

This document provides instructions for building and installing the Seurat Shiny Browser package from source.

## Prerequisites

Before you begin, ensure you have the following installed:

1. R (version ≥ 4.0.0)
2. RStudio (recommended, but not required)
3. devtools package
4. roxygen2 package
5. Required dependencies:
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

## Installation from GitHub

The easiest way to install the package is directly from GitHub:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install the package with all required dependencies
# This will automatically install or update all necessary packages with compatible versions
devtools::install_github("davidvi/seurat-shiny-browser")
```

The package depends on the following R packages (specific version requirements are handled automatically):

- **Core**: Seurat (>= 5.0.0), SeuratObject (>= 5.0.0)
- **UI**: shiny, shinyWidgets, shinyjs, DT
- **Data manipulation**: tidyverse, dplyr, tidyr, purrr, etc.
- **Visualization**: ggplot2
- **File handling**: readxl, readr
- **Processing**: preprocessCore, future

If you encounter any issues with dependencies during installation, you can install them manually:

```r
# Core packages
install.packages("Seurat")

# UI packages
install.packages(c("shiny", "shinyWidgets", "shinyjs", "DT"))

# Data manipulation and visualization
install.packages(c("tidyverse", "ggplot2"))

# Other utilities
install.packages(c("future", "readxl"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("preprocessCore")
```

## Building from Source

If you prefer to build the package from source:

1. Clone the repository:

```bash
git clone https://github.com/davidvi/seurat-shiny-browser.git
cd seurat-shiny-browser
```

2. Install required development tools:

```r
install.packages(c("devtools", "roxygen2", "pkgdown", "testthat"))
```

3. Generate documentation:

```r
devtools::document()
```

4. Build and install the package:

```r
devtools::install()
```

## Usage

After installation, you can launch the application with:

```r
library(seuratShinyBrowser)

# Launch with default port (3030)
run_app(data_folder = "/path/to/your/data")

# Or specify a custom port (useful if 3030 is already in use)
run_app(data_folder = "/path/to/your/data", port = 8080)
```

You can then access the application in your web browser at:
- Default port: http://localhost:3030 
- Custom port: http://localhost:8080 (or whatever port you specified)

## Common Issues

- **Missing dependencies**: If you encounter missing dependency errors, install the required packages manually.
- **Permission errors**: Make sure you have write access to the R library locations.
- **Documentation errors**: Run `devtools::document()` before installation.
- **Path issues**: Use absolute paths for the data folder to avoid unexpected behavior.

## Customizing the Package

If you want to modify the package:

1. Make your changes to the R files in the `R/` directory or the Shiny app files in `inst/shiny-app/`
2. Run `devtools::document()` to update documentation
3. Rebuild and reinstall the package with `devtools::install()`

## License

This package is licensed under the MIT License.