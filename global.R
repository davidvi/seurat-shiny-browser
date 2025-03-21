suppressPackageStartupMessages({
  library(Seurat)
  library(shiny)
  library(shinyWidgets)
  library(shinyjs)
  library(readxl)
  library(preprocessCore)
  library(ggplot2)
  library(tidyverse)
  library(DT)
  library(future)
})

# load data
base_folder <- "/home/shiny-app/data/"
all_rds_files <- list.files(path = base_folder, pattern = ".rds", full.names = FALSE)

sample <- NULL
sample_name <- NULL
cluster_names <- NULL

load_sample <- function(sample_file) {
  sample <- readRDS(file = paste0(base_folder, sample_file))
  tryCatch(
    {
      DefaultAssay(sample) <- "RNA"
    },
    error = function(e) {}
  )
  return(sample)
}

sample <- load_sample(all_rds_files[1])
sample_name <- all_rds_files[1]
cluster_names <- levels(Idents(sample))
