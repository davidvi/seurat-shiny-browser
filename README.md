# Seurat Shiny Browser

A comprehensive interactive web application for exploring and analyzing single-cell RNA-seq data using the Seurat framework.

![Seurat Shiny Browser](https://img.shields.io/badge/Seurat-Shiny-blue)
![R](https://img.shields.io/badge/R-4.2.2-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

## Overview

Seurat Shiny Browser is a powerful tool designed for biologists, bioinformaticians, and researchers working with single-cell RNA sequencing data. It provides an intuitive interface for performing common Seurat operations without writing code, including:

- Loading and visualizing Seurat objects or raw 10X Genomics data
- Exploring gene expression across different cell populations
- Performing dimensional reduction (PCA, UMAP) and visualization
- Executing step-by-step normalization and clustering workflows
- Finding marker genes through differential expression analysis
- Integrating and batch-correcting multiple samples
- Organizing and managing datasets across multiple folders

## Installation

### Option 1: Installation as an R Package (New)

You can now install Seurat Shiny Browser as a standard R package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install the package with all required dependencies (v5.0.0+ of Seurat and other packages)
devtools::install_github("davidvi/seurat-shiny-browser")

# Launch the application
library(seuratShinyBrowser)  # Package name may differ from the repository name

# Default port (3030)
run_seurat_browser(data_folder = "/path/to/your/data")

# Or specify a custom port
run_seurat_browser(data_folder = "/path/to/your/data", port = 8080)
```

The package will automatically install all required dependencies with compatible versions, including:
- Seurat (>= 5.0.0) and SeuratObject (>= 5.0.0)
- All necessary Shiny and visualization packages
- Data manipulation libraries (tidyverse ecosystem)

See R_PACKAGE_INSTALL.md for more detailed installation instructions.

### Option 2: Using Docker

The easiest way to run Seurat Shiny Browser is with Docker, which eliminates the need to install R packages directly.

#### Prerequisites

- [Docker](https://www.docker.com/products/docker-desktop/) installed on your system

#### Steps

1. Clone the repository:

```bash
git clone https://github.com/davidvi/seurat-shiny-browser.git
cd seurat-shiny-browser
```

2. Build the Docker image:

```bash
./build.sh
```

3. Run the application:

```bash
# Default settings
./run.sh

# Or with custom container data folder, port, and local data directory
./run.sh /home/shiny-app/data 8080 /path/to/your/local/data
```

4. Access the application in your browser at:

```
http://localhost:3030  # If using default port
http://localhost:8080  # If using custom port
```

The run.sh script accepts three optional parameters:
1. Container data path (default: /home/shiny-app/data)
2. Port number (default: 3030)
3. Local data directory to mount (default: ./data)

### Option 3: Running Directly in R

#### Prerequisites

- R (version ≥ 4.2.0)
- Required R packages:
  - Seurat (≥ 5.0.0)
  - SeuratObject (≥ 5.0.0)
  - shiny
  - shinyWidgets
  - shinyjs
  - ggplot2
  - magrittr
  - dplyr, tidyr, and other tidyverse packages
  - DT
  - future
  - readxl
  - preprocessCore

#### Steps

1. Clone the repository:

```bash
git clone https://github.com/davidvi/seurat-shiny-browser.git
cd seurat-shiny-browser
```

2. Install required packages:

```R
install.packages(c("shiny", "shinyWidgets", "shinyjs", "ggplot2", "magrittr", 
                   "dplyr", "tidyr", "purrr", "stringr", "readr", "forcats",
                   "DT", "future", "readxl", "tibble", "lubridate"))
install.packages("BiocManager")
BiocManager::install("preprocessCore")
install.packages("Seurat")
```

3. Run the application:

```bash
# Method 1: Use R/dev/dev_run.R with command-line arguments (RECOMMENDED)
Rscript R/dev/dev_run.R /path/to/your/data         # Specify data folder
Rscript R/dev/dev_run.R /path/to/your/data 8080    # Specify data folder and port
Rscript R/dev/dev_run.R                            # Use default data folder (./data)

# Method 2: Run directly with the package functions
Rscript -e "devtools::load_all(); run_seurat_browser(data_folder = '/path/to/your/data', port = 3030)"
```

4. Access the application in your browser at:

```
http://localhost:3030
```

## Usage

### Setting Up Your Data

The app works with two types of data:

1. **Processed Seurat Objects (.rds files)**  
   Place your .rds files in subdirectories within the data folder.

2. **Raw 10X Genomics Data**  
   Create a folder structure with 10X output (matrix.mtx, features.tsv, barcodes.tsv).

### App Workflow

#### 1. Load/Save Tab

- **Loading Data**:
  - Select a folder and load a processed Seurat object
  - Import raw 10X Genomics data with customizable parameters
  
- **File Management**:
  - Create new folders for organizing your data
  - Move files between folders
  - Delete samples when needed
  - Save the current sample with a new name

#### 2. Visualize Tab

- **Single Gene Visualization**:
  - Search for specific genes of interest
  - View violin plots and feature plots simultaneously
  
- **Multi-Gene Visualization**:
  - Enter multiple genes for comparative analysis
  - View expression patterns across clusters with dotplots
  
- **Metadata Visualization**:
  - Color cells by any metadata column
  - Inspect dimensional reductions (UMAP, t-SNE, etc.)

#### 3. Edit/Rename Clusters Tab

- **Cluster Management**:
  - Rename single or multiple clusters
  - Back up and restore cluster identities
  - Delete unwanted clusters

#### 4. Differential Expression Tab

- **Marker Gene Discovery**:
  - Find markers for specific clusters
  - Compare two clusters for differential expression
  - Find markers across all clusters at once
  - Download results as CSV files

#### 5. Integrate Samples Tab

- **Sample Integration**:
  - Select multiple samples to merge
  - Choose integration methods (CCA, RPCA, Harmony, Joint PCA)
  - Run clustering and UMAP on integrated data

#### 6. Normalize/Cluster Tab

- **Step-by-Step Analysis**:
  - QC filtering with customizable thresholds
  - Normalization (LogNormalize or SCTransform)
  - Variable feature selection
  - Scaling and PCA
  - Cluster identification with adjustable parameters
  - Dimensional reduction (UMAP)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- The [Seurat](https://satijalab.org/seurat/) team for their excellent single-cell analysis framework
- The R Shiny development team
- Contributors to the open-source packages this application relies on