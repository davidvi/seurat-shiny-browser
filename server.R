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
  library(shinyjs)
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
all_fovs <- NULL
all_fovs <- unique(sample@meta.data$fov)
fov <- all_fovs[1]

# UI

jsCode <- "
$(document).on('keypress', '#gene', function(e) {
  if (e.which == 13) {  // Enter key = keycode 13
    $('#search_button').click();  // trigger search_button click event
  }
});
"

ui <- fluidPage(
  useShinyjs(),
  tags$script(jsCode),
  titlePanel("Single Cell Data Explorer"),
  sidebarLayout(
    sidebarPanel(
      h3("Search for a gene"),
      textInput(
        inputId = "gene",
        label = "Gene Name:",
        value = "CSF1R"
      ),
      actionButton("search_button", "Search gene"),
      h3("Load another sample"),
      selectInput(
        inputId = "sample",
        label = "Select Sample:",
        all_rds_files
      ),
      actionButton("load_button", "Load file"),
      h3("Actions on cell clusters"),
      selectInput(
        inputId = "cluster_name_selector",
        label = "Select cluster:",
        choices = cluster_names
      ),
      actionButton("find_markers", "Find markers for selected cluster"),
      actionButton("find_all_markers", "Find markers for all clusters"),
      textInput(
        inputId = "new_cluster_name",
        label = "Rename cluster to:",
        value = ""
      ),
      actionButton("rename_button", "Rename cluster"),
      selectInput(
        inputId = "second_cluster_name_selector",
        label = "Select second cluster:",
        choices = cluster_names
      ),
      actionButton("compare_button", "Differential analysis two cluster"),
      h3("Actions on CoxMX data FOVs"),
      selectInput(
        inputId = "fov",
        label = "Select FOV:",
        choices = all_fovs
      ),
      actionButton("search_fov", "Load FOV"),
      h3("Save the current sample"),
      textInput(
        inputId = "new_save_name",
        label = "Save file as:",
        value = ""
      ),
      actionButton("save_button", "Save file")
    ),
    mainPanel(
      h4(textOutput("sample_name")),

      # First row of plots
      fluidRow(
        column(6, plotOutput(outputId = "violin_plot_sample1")),
        column(6, plotOutput(outputId = "feature_plot_sample1"))
      ),
      h4(textOutput("fov_name")),

      # Second row of plots
      fluidRow(
        column(6, plotOutput(outputId = "cosmx_cell_counts")),
        column(6, plotOutput(outputId = "cosmx_overlay_transcripts"))
      ),
      dataTableOutput("markers_table"),
      downloadButton("download_data", "Download table data")
    )
  )
)

# server

server <- function(input, output, session) {
  rv <- reactiveValues(
    sample = sample,
    sample_name = sample_name,
    cluster_names = cluster_names,
    gene = "CSF1R",
    all_rds_files = all_rds_files,
    markers = data.frame(),
    all_fovs = all_fovs,
    fov = all_fovs[1]
  )


  ### rename cluster ###

  observeEvent(input$rename_button, {
    if (input$new_cluster_name != "") {
      cells.use <- WhichCells(rv$sample, idents = input$cluster_name_selector)
      rv$sample <- SetIdent(rv$sample, cells = cells.use, value = input$new_cluster_name)
      rv$cluster_names <- levels(Idents(rv$sample))
    } else {
      showNotification("cannot be empty")
    }
  })

  ### search gene ###

  observeEvent(input$search_button, {
    rv$gene <- toupper(input$gene)
    rv$fov <- input$fov
  })

  ### search fov ###

  observeEvent(input$search_fov, {
    rv$gene <- toupper(input$gene)
    rv$fov <- input$fov
  })

  ### find clusters ###

  observeEvent(input$find_markers, {
    withProgress(message = "Finding markers for cluster", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Findings markers")
      plan("multicore", workers = 10)
      rv$markers <- FindMarkers(rv$sample, only.pos = TRUE, logfc.threshold = 0.25, ident.1 = input$cluster_name_selector)
      incProgress(0.9, "Done")
    })
  })

  ### compare clusters ###

  observeEvent(input$compare_button, {
    withProgress(message = "DE analysis for cluster", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Finding markers")
      plan("multicore", workers = 10)
      rv$markers <- FindMarkers(rv$sample, only.pos = FALSE, logfc.threshold = 0.25, ident.1 = input$cluster_name_selector, ident.2 = input$second_cluster_name_selector)
      incProgress(0.9, "Done")
    })
  })

  ### find all markers ###

  observeEvent(input$find_all_markers, {
    withProgress(message = "Finding all markers", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Finding markers")
      plan("multicore", workers = 10)
      rv$markers <- FindAllMarkers(rv$sample, only.pos = TRUE, logfc.threshold = 0.25)
      incProgress(0.9, "Done")
    })
  })

  ### download table ###

  output$download_data <- downloadHandler(
    filename = function() {
      paste0(rv$sample_name, "_markers.csv")
    },
    content = function(file) {
      write.csv(rv$markers, file)
    }
  )

  ### Load sample ###

  observeEvent(input$load_button, {
    withProgress(message = "Loading Sample", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Loading Sample")
      rv$sample <- load_sample(input$sample)
      incProgress(0.5, "Set sample name")
      rv$sample_name <- input$sample
      rv$cluster_names <- levels(Idents(rv$sample))
      rv$markers <- NULL
      rv$all_fovs <- unique(rv$sample@meta.data$fov)
      rv$fov <- rv$all_fovs[1]
      incProgress(0.9, "Done")
    })
  })

  ### save as new ###

  observeEvent(input$save_button, {
    if (endsWith(input$new_save_name, ".rds")) {
      withProgress(message = "Saving Sample", detail = "Please wait...", value = 0, {
        incProgress(0.1, "Saving Sample")
        saveRDS(rv$sample, file = paste0(base_folder, input$new_save_name))
        incProgress(0.8, "Setting new name")
        rv$sample_name <- input$new_save_name
        rv$all_rds_files <- list.files(path = base_folder, pattern = ".rds", full.names = FALSE)
        updateSelectInput(session, "sample", choices = rv$all_rds_files, selected = rv$sample_name)
        incProgress(0.9, "Done")
      })
    } else {
      showNotification("cannot be empty and must end in .rds")
    }
  })

  observe({
    rv$all_rds_files <- list.files(path = base_folder, pattern = ".rds", full.names = FALSE)

    updateSelectInput(session, "sample", choices = rv$all_rds_files, selected = rv$sample_name)

    updateSelectInput(session, "fov", choices = rv$all_fovs)

    updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
    updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)

    output$sample_name <- renderText({
      paste("Loaded Sample: ", rv$sample_name)
    })
    output$violin_plot_sample1 <- renderPlot({
      VlnPlot(rv$sample, features = c(rv$gene), pt.size = 0) + NoLegend() + xlab("")
    })
    output$feature_plot_sample1 <- renderPlot({
      FeaturePlot(rv$sample, features = rv$gene, label = T) + NoLegend()
    })
    output$markers_table <- renderDataTable({
      rv$markers
    })

    ### CosMX plots ###

    output$fov_name <- renderText({
      paste("Loaded FOV: ", rv$fov)
    })

    output$cosmx_cell_counts <- renderPlot({
      # Check if 'fov' is a column in the metadata
      if ("fov" %in% colnames(rv$sample@meta.data)) {
        # Render the plot if the condition is true
        ImageFeaturePlot(subset(x = rv$sample, subset = fov == rv$fov),
          features = rv$gene,
          size = 0.75, cols = c("white", "red")
        )
      } else {
        # Optional: Display a message or an alternative output if the condition is false
        plot.new()
        text(0.5, 0.5, "FOV column not found in metadata", cex = 1.2)
      }
    })
    output$cosmx_overlay_transcripts <- renderPlot({
      # Check if 'fov' is a column in the metadata
      if ("fov" %in% colnames(rv$sample@meta.data)) {
        # Render the plot if the condition is true
        ImageDimPlot(subset(x = rv$sample, subset = fov == rv$fov),
          cols = "polychrome",
          alpha = 0.3,
          molecules = rv$gene,
          mols.size = 0.3,
          crop = TRUE,
          nmols = 20000,
          border.color = "black",
          coord.fixed = FALSE
        )
      } else {
        # Optional: Display a message or an alternative output if the condition is false
        plot.new()
        text(0.5, 0.5, "FOV column not found in metadata", cex = 1.2)
      }
    })
  })
}

### run app ###

options(shiny.autoreload = TRUE)
options(shiny.host = "0.0.0.0")
options(shiny.port = 3030)
shinyApp(ui, server)
