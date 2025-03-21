# Visualization module

# Visualize tab plots
output$violin_plot_visualize <- renderPlot({
  VlnPlot(rv$sample, features = c(rv$gene), pt.size = 0) + NoLegend() + xlab("")
})
output$feature_plot_visualize <- renderPlot({
  tryCatch({
    FeaturePlot(rv$sample, features = rv$gene, reduction = rv$selected_reduction, label = TRUE) + NoLegend()
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cannot display feature on", rv$selected_reduction, "\nError:", e$message), cex = 1.2)
  })
})
output$dim_plot_visualize <- renderPlot({
  tryCatch({
    if(input$metadata_column_selector == "ident") {
      # For identity class, use DimPlot
      DimPlot(rv$sample, reduction = rv$selected_reduction, label = TRUE)
    } else {
      # Check if the selected metadata is continuous or categorical
      meta_values <- rv$sample@meta.data[[input$metadata_column_selector]]
      
      # If it's numeric and has more than a few unique values, treat as continuous
      if(is.numeric(meta_values) && length(unique(meta_values)) > 10) {
        # For continuous data, use FeaturePlot
        FeaturePlot(rv$sample, 
                  features = input$metadata_column_selector, 
                  reduction = rv$selected_reduction)
      } else {
        # For categorical data, use DimPlot
        DimPlot(rv$sample, 
               reduction = rv$selected_reduction, 
               group.by = input$metadata_column_selector)
      }
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cannot display:", input$metadata_column_selector, "on", rv$selected_reduction, "\nError:", e$message), cex = 1.2)
  })
})

# Edit tab plots
output$violin_plot_edit <- renderPlot({
  VlnPlot(rv$sample, features = c(rv$gene), pt.size = 0) + NoLegend() + xlab("")
})
output$feature_plot_edit <- renderPlot({
  tryCatch({
    FeaturePlot(rv$sample, features = rv$gene, reduction = rv$selected_reduction, label = TRUE) + NoLegend()
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cannot display feature on", rv$selected_reduction, "\nError:", e$message), cex = 1.2)
  })
})

# Differential Expression tab plots
output$violin_plot_analyze <- renderPlot({
  VlnPlot(rv$sample, features = c(rv$gene), pt.size = 0) + NoLegend() + xlab("")
})
output$feature_plot_analyze <- renderPlot({
  tryCatch({
    FeaturePlot(rv$sample, features = rv$gene, reduction = rv$selected_reduction, label = TRUE) + NoLegend()
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Cannot display feature on", rv$selected_reduction, "\nError:", e$message), cex = 1.2)
  })
})