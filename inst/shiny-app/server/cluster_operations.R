# Cluster operations module

### rename cluster ###

observeEvent(input$rename_button, {
  # Validate input
  if (input$new_cluster_name == "") {
    showNotification("New cluster name cannot be empty", type = "error", duration = 10)
    return()
  }
  
  # Get selected clusters
  selected_clusters <- input$cluster_name_selector
  
  # Validate that at least one cluster is selected
  if (length(selected_clusters) == 0) {
    showNotification("Please select at least one cluster to rename", type = "warning", duration = 10)
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
      
      showNotification(paste("Cluster", selected_clusters, "renamed to", input$new_cluster_name), type = "message", duration = 10)
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

### rename orig.ident ###

observeEvent(input$rename_orig_ident_button, {
  # Validate inputs
  req(rv$sample)
  
  # Check if orig.ident column exists
  if (!("orig.ident" %in% colnames(rv$sample@meta.data))) {
    showNotification("This Seurat object does not have an orig.ident column", type = "error", duration = 10)
    return()
  }
  
  # Get selected orig.ident value
  selected_orig_ident <- input$orig_ident_selector
  
  # Get new orig.ident name
  new_orig_ident <- input$new_orig_ident
  
  # Validate new name
  if (is.null(new_orig_ident) || new_orig_ident == "") {
    showNotification("New orig.ident name cannot be empty", type = "error", duration = 10)
    return()
  }
  
  # Start the renaming process
  withProgress(message = "Renaming orig.ident", detail = "Please wait...", value = 0, {
    tryCatch({
      # Get current unique values before the change
      original_values <- unique(as.character(rv$sample$orig.ident))
      
      # Identify cells with the selected orig.ident value
      cells_to_rename <- rownames(rv$sample@meta.data)[rv$sample$orig.ident == selected_orig_ident]
      
      if (length(cells_to_rename) == 0) {
        showNotification("No cells found with the selected orig.ident value", type = "warning", duration = 10)
        return()
      }
      
      incProgress(0.3, detail = paste("Renaming", length(cells_to_rename), "cells..."))
      
      # First convert to character to avoid factor issues
      rv$sample$orig.ident <- as.character(rv$sample$orig.ident)
      
      # Replace the orig.ident values for the selected cells
      rv$sample$orig.ident[cells_to_rename] <- new_orig_ident
      
      # Get the new unique values after replacement
      new_unique_values <- unique(rv$sample$orig.ident)
      
      # Convert back to factor with explicit levels
      rv$sample$orig.ident <- factor(rv$sample$orig.ident, levels = new_unique_values)
      
      # Force regenerate any dependent reduction plots
      if ("DimPlot" %in% names(rv$sample@misc)) {
        rv$sample@misc$DimPlot <- NULL
      }
      
      incProgress(0.8, detail = "Updating UI...")
      
      # Update the orig.ident selector with the new set of unique values
      updateSelectInput(session, "orig_ident_selector", choices = levels(rv$sample$orig.ident))
      
      # Force a refresh of any visualization that might use this metadata
      rv$force_refresh <- runif(1) # Random value to trigger reactivity
      
      # Show success notification
      showNotification(
        paste("Successfully renamed orig.ident from", selected_orig_ident, "to", new_orig_ident, 
              "for", length(cells_to_rename), "cells"),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error renaming orig.ident:", e$message), type = "error", duration = 10)
    })
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
    showNotification("No backup found. Please create a backup first.", type = "error", duration = 10)
  }
})

### delete cluster ###

observeEvent(input$delete_cluster_button, {
  # Check if deletion is confirmed
  if (!input$confirm_delete_cluster) {
    showNotification("Please confirm deletion by checking the confirmation box", type = "warning", duration = 10)
    return()
  }
  
  # Get the clusters to delete
  clusters_to_delete <- input$delete_cluster_selector
  
  # Check if any clusters were selected
  if (length(clusters_to_delete) == 0) {
    showNotification("Please select at least one cluster to delete", type = "warning", duration = 10)
    return()
  }
  
  # Verify all selected clusters exist
  all_clusters <- levels(Idents(rv$sample))
  invalid_clusters <- clusters_to_delete[!clusters_to_delete %in% all_clusters]
  if (length(invalid_clusters) > 0) {
    showNotification(paste("Some selected clusters do not exist:", paste(invalid_clusters, collapse=", ")), type = "error", duration = 10)
    return()
  }
  
  # Check if we would delete all clusters
  if (length(clusters_to_delete) == length(all_clusters)) {
    showNotification("Cannot delete all clusters. At least one cluster must remain.", type = "error", duration = 10)
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
        showNotification("Cannot delete all cells in the dataset. At least one cell must remain.", type = "error", duration = 10)
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
      showNotification(paste("Error deleting clusters:", e$message), type = "error", duration = 10)
    })
  })
})

### calculate log2UMI_cell ###

observeEvent(input$add_log2umi_button, {
  # Check if a sample is loaded
  req(rv$sample)
  
  withProgress(message = "Calculating log2UMI_cell", detail = "Please wait...", value = 0, {
    tryCatch({
      incProgress(0.3, detail = "Checking nCount_RNA exists...")
      
      # Check if nCount_RNA exists
      if (!("nCount_RNA" %in% colnames(rv$sample@meta.data))) {
        showNotification("This Seurat object does not have an nCount_RNA column", type = "error", duration = 10)
        return()
      }
      
      incProgress(0.5, detail = "Calculating log2UMI values...")
      
      # Calculate log2UMI_cell
      rv$sample@meta.data$log2UMI_cell <- log2(rv$sample@meta.data$nCount_RNA + 1)
      
      incProgress(0.8, detail = "Updating UI...")
      
      # Update the metadata columns list
      rv$metadata_columns <- colnames(rv$sample@meta.data)
      
      # Update any UI elements that use metadata columns
      if("metadata_column_selector" %in% names(input)) {
        updateSelectInput(session, "metadata_column_selector", 
                        choices = c("ident", rv$metadata_columns),
                        selected = input$metadata_column_selector)
      }
      
      # Update visualization options if they exist
      if("visualization_color_by" %in% names(input)) {
        old_selection <- input$visualization_color_by
        updateSelectInput(session, "visualization_color_by", 
                        choices = c("ident", rv$metadata_columns),
                        selected = old_selection)
      }
      
      # Force visualization refresh
      rv$force_refresh <- runif(1)
      
      # Show success notification
      showNotification("Successfully calculated log2UMI_cell metadata column", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Error calculating log2UMI_cell:", e$message), type = "error", duration = 10)
    })
  })
})
