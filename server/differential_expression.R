# Differential expression analysis module

### find clusters ###

observeEvent(input$find_markers, {
  # Validate we have a sample loaded
  if(is.null(rv$sample)) {
    showNotification("No sample loaded. Please load a sample first.", type = "error", duration = 10)
    return()
  }
  
  # Validate we have clusters
  if(length(levels(Idents(rv$sample))) == 0) {
    showNotification("No clusters found in sample. Please run clustering first.", type = "error", duration = 10)
    return()
  }
  
  # Validate a cluster is selected
  if(is.null(input$cluster_name_selector_analyze) || input$cluster_name_selector_analyze == "") {
    showNotification("Please select a cluster for analysis", type = "warning", duration = 10)
    return()
  }
  
  withProgress(message = "Finding markers for cluster", detail = "Please wait...", value = 0, {
    incProgress(0.1, "Finding markers")
    
    tryCatch({
      # Initialize markers data frame
      rv$markers <- data.frame()
      
      # Use future plan for parallelization
      plan("multicore", workers = 4)  # Reduced from 10 to 4 for stability
      
      # Run FindMarkers
      markers_result <- FindMarkers(
        rv$sample, 
        only.pos = TRUE, 
        logfc.threshold = 0.25, 
        ident.1 = input$cluster_name_selector_analyze
      )
      
      # If there are results, add a gene column and update rv$markers
      if(!is.null(markers_result) && nrow(markers_result) > 0) {
        # Add gene column if not present
        if(!("gene" %in% colnames(markers_result))) {
          markers_result$gene <- rownames(markers_result)
          # Reorder columns to put gene first
          markers_result <- markers_result[, c("gene", setdiff(colnames(markers_result), "gene"))]
        }
        
        # Update reactive value
        rv$markers <- markers_result
        
        # Show success message
        incProgress(0.9, "Done")
        showNotification(
          paste("Found", nrow(markers_result), "marker genes for cluster", input$cluster_name_selector_analyze),
          type = "message",
          duration = 10
        )
      } else {
        showNotification(
          paste("No marker genes found for cluster", input$cluster_name_selector_analyze, "with current threshold settings"),
          type = "warning",
          duration = 10
        )
      }
      
    }, error = function(e) {
      showNotification(
        paste("Error finding markers:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
})

### compare clusters ###

observeEvent(input$compare_button, {
  # Validate we have a sample loaded
  if(is.null(rv$sample)) {
    showNotification("No sample loaded. Please load a sample first.", type = "error", duration = 10)
    return()
  }
  
  # Validate clusters are selected
  if(is.null(input$first_cluster_name_selector) || is.null(input$second_cluster_name_selector)) {
    showNotification("Please select two clusters to compare", type = "warning", duration = 10)
    return()
  }
  
  # Validate different clusters are selected
  if(input$first_cluster_name_selector == input$second_cluster_name_selector) {
    showNotification("Please select two different clusters to compare", type = "warning", duration = 10)
    return()
  }
  
  withProgress(message = "Running differential expression analysis", detail = "Please wait...", value = 0, {
    incProgress(0.1, "Finding differentially expressed genes")
    
    tryCatch({
      # Initialize markers data frame
      rv$markers <- data.frame()
      
      # Use future plan for parallelization
      plan("multicore", workers = 4)  # Reduced from 10 to 4 for stability
      
      # Run FindMarkers for comparison
      markers_result <- FindMarkers(
        rv$sample, 
        only.pos = FALSE, 
        logfc.threshold = 0.25, 
        ident.1 = input$first_cluster_name_selector, 
        ident.2 = input$second_cluster_name_selector
      )
      
      # If there are results, add a gene column and update rv$markers
      if(!is.null(markers_result) && nrow(markers_result) > 0) {
        # Add gene column if not present
        if(!("gene" %in% colnames(markers_result))) {
          markers_result$gene <- rownames(markers_result)
          # Reorder columns to put gene first
          markers_result <- markers_result[, c("gene", setdiff(colnames(markers_result), "gene"))]
        }
        
        # Update reactive value
        rv$markers <- markers_result
        
        # Show success message
        incProgress(0.9, "Done")
        showNotification(
          paste("Found", nrow(markers_result), "differentially expressed genes between clusters", 
                input$first_cluster_name_selector, "and", input$second_cluster_name_selector),
          type = "message",
          duration = 10
        )
      } else {
        showNotification(
          paste("No differentially expressed genes found between clusters", 
                input$first_cluster_name_selector, "and", input$second_cluster_name_selector, 
                "with current threshold settings"),
          type = "warning",
          duration = 10
        )
      }
      
    }, error = function(e) {
      showNotification(
        paste("Error in differential expression analysis:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
})

### find all markers ###

observeEvent(input$find_all_markers, {
  # Validate we have a sample loaded
  if(is.null(rv$sample)) {
    showNotification("No sample loaded. Please load a sample first.", type = "error", duration = 10)
    return()
  }
  
  # Validate we have clusters
  if(length(levels(Idents(rv$sample))) <= 1) {
    showNotification("Need at least 2 clusters for comparison. Please run clustering first.", type = "error", duration = 10)
    return()
  }
  
  withProgress(message = "Finding markers for all clusters", detail = "This may take a while...", value = 0, {
    incProgress(0.1, "Finding all markers")
    
    tryCatch({
      # Initialize markers data frame
      rv$markers <- data.frame()
      
      # Use future plan for parallelization
      plan("multicore", workers = 4)  # Reduced from 10 to 4 for stability
      
      # Run FindAllMarkers
      markers_result <- FindAllMarkers(
        rv$sample, 
        only.pos = TRUE, 
        logfc.threshold = 0.25,
        min.pct = 0.1
      )
      
      # If there are results, add a gene column and update rv$markers
      if(!is.null(markers_result) && nrow(markers_result) > 0) {
        # Add gene column if not present (FindAllMarkers should already have this)
        if(!("gene" %in% colnames(markers_result)) && "features" %in% colnames(markers_result)) {
          # Rename features to gene
          colnames(markers_result)[colnames(markers_result) == "features"] <- "gene"
        } else if(!("gene" %in% colnames(markers_result)) && !("features" %in% colnames(markers_result))) {
          # Add gene column from rownames if needed
          markers_result$gene <- rownames(markers_result)
          # Reorder columns to put gene first
          markers_result <- markers_result[, c("gene", setdiff(colnames(markers_result), "gene"))]
        }
        
        # Update reactive value
        rv$markers <- markers_result
        
        # Show success message
        incProgress(0.9, "Done")
        showNotification(
          paste("Found", nrow(markers_result), "marker genes across all clusters"),
          type = "message",
          duration = 10
        )
      } else {
        showNotification(
          "No marker genes found with current threshold settings",
          type = "warning",
          duration = 10
        )
      }
      
    }, error = function(e) {
      showNotification(
        paste("Error finding all markers:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
})

### download table ###

output$download_data <- downloadHandler(
  filename = function() {
    paste0(rv$sample_name, "_markers_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    # Check if markers exist
    if(is.null(rv$markers) || nrow(rv$markers) == 0) {
      # Create an empty file with a message
      write.csv(data.frame(message = "No marker genes found"), file)
    } else {
      # Write the markers to CSV
      write.csv(rv$markers, file, row.names = FALSE)
    }
  }
)

# Render the markers table
output$markers_table <- renderDataTable({
  # Check if markers exist
  if(is.null(rv$markers) || nrow(rv$markers) == 0) {
    return(data.frame(message = "No marker genes found yet. Run one of the analyses above."))
  }
  
  # Return the markers table
  rv$markers
}, options = list(
  pageLength = 10,
  scrollX = TRUE,
  scrollY = "400px"
))