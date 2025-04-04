# This file contains package-level functions that are called when the package is loaded

#' @importFrom DT datatable
#' @importFrom SeuratObject DefaultAssay
#' @importFrom dplyr filter select mutate arrange summarize group_by
#' @importFrom purrr map map_df keep
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect str_replace
#' @importFrom forcats fct_reorder
#' @importFrom lubridate ymd_hms
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom readxl read_excel
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom sp CRS
#' @importFrom shinyWidgets pickerInput
#' @importFrom shinyjs useShinyjs
#' @importFrom ggplot2 ggplot
#' @importFrom future plan
#' @importFrom magrittr %>%
#' @import Seurat
#' @import shiny
NULL

# Dummy function to satisfy CRAN check
.onLoad <- function(libname, pkgname) {
  # This is just to prevent R CMD check from complaining about
  # no visible binding for global variable
  utils::globalVariables(c("."))
  invisible()
}