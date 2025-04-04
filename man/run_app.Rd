\name{run_app}
\alias{run_app}
\title{Run Seurat Shiny Browser App}
\usage{
run_app(data_folder = NULL, port = 3030, launch.browser = TRUE, ...)
}
\arguments{
\item{data_folder}{Path to the folder containing Seurat data files. If NULL, a default path will be used.}

\item{port}{The port for the Shiny app to run on. Default is 3030. Change this if the default port is already in use.}

\item{launch.browser}{If TRUE, the system's default web browser will be launched automatically.}

\item{...}{Additional arguments to pass to shiny::runApp.}
}
\description{
This function launches the Seurat Shiny Browser application.
}
\examples{
\dontrun{
# Run with default settings
run_app(data_folder = "~/my_seurat_data/")

# Run on a different port
run_app(data_folder = "~/my_seurat_data/", port = 8080)
}
}