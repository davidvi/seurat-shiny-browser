# server

server <- function(input, output, session) {
  rv <- reactiveValues(
    sample = sample,
    sample_name = sample_name,
    cluster_names = cluster_names,
    gene = "CSF1R",
    all_rds_files = all_rds_files,
    markers = data.frame()
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
  })
}
