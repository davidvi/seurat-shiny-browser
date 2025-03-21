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

### Option 1: Using Docker (Recommended)

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
./run.sh
```

4. Access the application in your browser at:

```
http://localhost:3030
```

### Option 2: Running Directly in R

#### Prerequisites

- R (version ≥ 4.0.0)
- Required R packages:
  - Seurat (≥ 4.3.0)
  - shiny
  - shinyWidgets
  - shinyjs
  - ggplot2
  - tidyverse
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
install.packages(c("shiny", "shinyWidgets", "shinyjs", "ggplot2", "tidyverse", "DT", "future", "readxl"))
install.packages("BiocManager")
BiocManager::install(c("preprocessCore", "Seurat"))
```

3. Run the application:

```bash
Rscript app.R /path/to/your/data
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
  - Choose integration methods (CCA, RPCA, Harmony, FastMNN)
  - Run clustering and UMAP on integrated data

#### 6. Normalize/Cluster Tab

- **Step-by-Step Analysis**:
  - QC filtering with customizable thresholds
  - Normalization (LogNormalize or SCTransform)
  - Variable feature selection
  - Scaling and PCA
  - Cluster identification with adjustable parameters
  - Dimensional reduction (UMAP)

### Custom Data Location and Port

You can specify a custom data directory and port either through Docker or direct execution:

- **With Docker**:

```bash
# Default port (3030)
./run.sh /path/to/your/data

# Custom port
./run.sh /path/to/your/data 8080
```

- **Direct with R**:

```bash
# Default port (3030)
Rscript app.R /path/to/your/data

# Custom port
Rscript app.R /path/to/your/data 8080
```

If you change the port, make sure to access the application at the correct URL (e.g., `http://localhost:8080` if you specified port 8080).

## Troubleshooting

### Common Issues

1. **App fails to start**:
   - Ensure the default port (3030) is not in use by another application
   - If you encounter a port conflict, you can specify a different port using the second command line argument
   - Check R and package versions are compatible

2. **No folders or files appear**:
   - Verify your data directory path is correct and accessible
   - Check folder permissions

3. **Analysis fails**:
   - Ensure your Seurat objects are compatible with the Seurat version in use
   - For large datasets, increase the memory allocation to Docker

4. **UI appears but plots don't render**:
   - Check the console logs for specific error messages
   - Verify that the loaded Seurat object contains the expected assays and reductions

## Contributing

Contributions to improve the Seurat Shiny Browser are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- The [Seurat](https://satijalab.org/seurat/) team for their excellent single-cell analysis framework
- The R Shiny development team
- Contributors to the open-source packages this application relies on

---

For additional questions or support, please [open an issue](https://github.com/davidvi/seurat-shiny-browser/issues) on the GitHub repository.