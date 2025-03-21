# Differential expression analysis module

### find clusters ###

observeEvent(input$find_markers, {
  withProgress(message = "Finding markers for cluster", detail = "Please wait...", value = 0, {
    incProgress(0.1, "Findings markers")
    plan("multicore", workers = 10)
    rv$markers <- FindMarkers(rv$sample, only.pos = TRUE, logfc.threshold = 0.25, ident.1 = input$cluster_name_selector_analyze)
    incProgress(0.9, "Done")
  })
})

### compare clusters ###

observeEvent(input$compare_button, {
  withProgress(message = "DE analysis for cluster", detail = "Please wait...", value = 0, {
    incProgress(0.1, "Finding markers")
    plan("multicore", workers = 10)
    rv$markers <- FindMarkers(rv$sample, only.pos = FALSE, logfc.threshold = 0.25, ident.1 = input$first_cluster_name_selector, ident.2 = input$second_cluster_name_selector)
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

# Render the markers table
output$markers_table <- renderDataTable({
  rv$markers
})