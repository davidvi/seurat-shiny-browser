# Cluster operations module

### rename cluster ###

observeEvent(input$rename_button, {
  # Validate input
  if (input$new_cluster_name == "") {
    showNotification("New cluster name cannot be empty", type = "error")
    return()
  }
  
  # Get selected clusters
  selected_clusters <- input$cluster_name_selector
  
  # Validate that at least one cluster is selected
  if (length(selected_clusters) == 0) {
    showNotification("Please select at least one cluster to rename", type = "warning")
    return()
  }
  
  # Begin renaming process
  withProgress(message = "Renaming clusters", detail = "Please wait...", value = 0, {
    
    # If only one cluster is selected, simple rename
    if (length(selected_clusters) == 1) {
      incProgress(0.3, detail = "Renaming cluster...")
      
      # Get cells from the selected cluster
      cells.use <- WhichCells(rv$sample, idents = selected_clusters)
      
      # Set new identity
      rv$sample <- SetIdent(rv$sample, cells = cells.use, value = input$new_cluster_name)
      
      showNotification(paste("Cluster", selected_clusters, "renamed to", input$new_cluster_name), type = "message")
    } 
    # Multiple clusters selected
    else {
      # Check if using numeric suffixes
      add_suffix <- input$add_number_suffix
      
      for (i in seq_along(selected_clusters)) {
        cluster <- selected_clusters[i]
        incProgress(i / length(selected_clusters), detail = paste("Renaming cluster", cluster, "..."))
        
        # Get cells from this cluster
        cells.use <- WhichCells(rv$sample, idents = cluster)
        
        # Create new name
        if (add_suffix) {
          new_name <- paste0(input$new_cluster_name, "_", i)
        } else {
          new_name <- input$new_cluster_name
        }
        
        # Set new identity
        rv$sample <- SetIdent(rv$sample, cells = cells.use, value = new_name)
      }
      
      if (add_suffix) {
        showNotification(
          paste("Renamed", length(selected_clusters), "clusters to", input$new_cluster_name, "with numeric suffixes"),
          type = "message"
        )
      } else {
        showNotification(
          paste("Renamed", length(selected_clusters), "clusters to", input$new_cluster_name),
          type = "message"
        )
      }
    }
    
    # Update cluster names
    rv$cluster_names <- levels(Idents(rv$sample))
    
    # Update all cluster selectors
    updatePickerInput(session, "cluster_name_selector", choices = rv$cluster_names)
    updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
    updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
    updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
    updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
  })
})

### backup cluster names ###

observeEvent(input$backup_cluster_button, {
  withProgress(message = "Backing up cluster identities", detail = "Please wait...", value = 0, {
    # Store current identities in a metadata column
    rv$sample$idents.backup <- Idents(rv$sample)
    
    # Update metadata column selector if present
    if("metadata_column_selector" %in% names(input)) {
      updateSelectInput(session, "metadata_column_selector", 
                      choices = c("ident", colnames(rv$sample@meta.data)),
                      selected = input$metadata_column_selector)
    }
    
    showNotification("Cluster identities backed up to 'idents.backup' metadata column", type = "message")
  })
})

### restore cluster names ###

observeEvent(input$restore_cluster_button, {
  # Check if the backup column exists
  if ("idents.backup" %in% colnames(rv$sample@meta.data)) {
    withProgress(message = "Restoring cluster identities", detail = "Please wait...", value = 0, {
      # Restore identities from the backup column
      rv$sample <- SetIdent(rv$sample, value = "idents.backup")
      
      # Update cluster names
      rv$cluster_names <- levels(Idents(rv$sample))
      
      # Update cluster selectors
      updatePickerInput(session, "cluster_name_selector", choices = rv$cluster_names)
      updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
      updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
      updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
      updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
      
      showNotification("Cluster identities restored from 'idents.backup'", type = "message")
    })
  } else {
    # No backup found
    showNotification("No backup found. Please create a backup first.", type = "error")
  }
})

### delete cluster ###

observeEvent(input$delete_cluster_button, {
  # Check if deletion is confirmed
  if (!input$confirm_delete_cluster) {
    showNotification("Please confirm deletion by checking the confirmation box", type = "warning")
    return()
  }
  
  # Get the clusters to delete
  clusters_to_delete <- input$delete_cluster_selector
  
  # Check if any clusters were selected
  if (length(clusters_to_delete) == 0) {
    showNotification("Please select at least one cluster to delete", type = "warning")
    return()
  }
  
  # Verify all selected clusters exist
  all_clusters <- levels(Idents(rv$sample))
  invalid_clusters <- clusters_to_delete[!clusters_to_delete %in% all_clusters]
  if (length(invalid_clusters) > 0) {
    showNotification(paste("Some selected clusters do not exist:", paste(invalid_clusters, collapse=", ")), type = "error")
    return()
  }
  
  # Check if we would delete all clusters
  if (length(clusters_to_delete) == length(all_clusters)) {
    showNotification("Cannot delete all clusters. At least one cluster must remain.", type = "error")
    return()
  }
  
  # Begin the deletion process
  withProgress(message = paste("Deleting", length(clusters_to_delete), "clusters"), detail = "Please wait...", value = 0, {
    tryCatch({
      # Count cells before deletion
      original_cell_count <- ncol(rv$sample)
      
      # Get all cells NOT in the clusters to delete (negative selection)
      cells_to_keep <- WhichCells(rv$sample, idents = clusters_to_delete, invert = TRUE)
      
      # Check if we would delete all cells
      if (length(cells_to_keep) == 0) {
        showNotification("Cannot delete all cells in the dataset. At least one cell must remain.", type = "error")
        return()
      }
      
      incProgress(0.3, detail = "Subsetting data...")
      
      # Create a new Seurat object with the cells to keep
      rv$sample <- subset(rv$sample, cells = cells_to_keep)
      
      # Count cells after deletion
      new_cell_count <- ncol(rv$sample)
      cells_removed <- original_cell_count - new_cell_count
      
      incProgress(0.7, detail = "Updating UI...")
      
      # Update the cluster names
      rv$cluster_names <- levels(Idents(rv$sample))
      
      # Update all cluster selectors
      updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
      updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
      updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
      updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
      updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
      
      # Reset the confirmation checkbox
      updateCheckboxInput(session, "confirm_delete_cluster", value = FALSE)
      
      # Show success message
      showNotification(
        paste0("Successfully deleted ", length(clusters_to_delete), " clusters (", 
               paste(clusters_to_delete, collapse=", "), "). ",
               cells_removed, " cells removed."), 
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error deleting clusters:", e$message), type = "error")
    })
  })
})
