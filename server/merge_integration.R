# Merge and integration module

### select merge folder ###

observeEvent(input$merge_folder_selector, {
  # Get files from the selected merge folder
  merge_folder <- input$merge_folder_selector
  merge_files <- get_rds_files(merge_folder)
  
  # Update the merge files selector
  updatePickerInput(session, "merge_files_selector", choices = merge_files)
  
  message("Updated merge file list for folder: ", merge_folder, ", found ", length(merge_files), " files")
})

# Merge tab outputs
output$merge_status <- renderText({
  if(is.null(input$merge_files_selector) || length(input$merge_files_selector) < 2) {
    "Please select at least two samples to merge."
  } else {
    paste(c(
      paste("Ready to merge", length(input$merge_files_selector), "samples:"),
      paste(input$merge_files_selector, collapse=", "),
      "",
      rv$merge_log
    ), collapse="\n")
  }
})

# Set up a flag for merge completion
output$merge_complete <- reactive({
  return(rv$merge_complete)
})
outputOptions(output, "merge_complete", suspendWhenHidden = FALSE)

# UMAP plot after integration
output$integrated_umap_plot <- renderPlot({
  req(rv$sample)
  req(rv$merge_complete)
  
  # Get the integrated reduction name
  if(rv$integrated_reduction != "") {
    # If we have a grouping variable or metadata column, use it
    if(!is.null(input$metadata_column_selector) && 
       input$metadata_column_selector %in% colnames(rv$sample@meta.data)) {
      DimPlot(rv$sample, reduction = rv$integrated_reduction, 
              group.by = input$metadata_column_selector, label = TRUE) + 
        ggtitle(paste("Integrated UMAP -", input$metadata_column_selector))
    } else {
      # Otherwise just show clusters
      DimPlot(rv$sample, reduction = rv$integrated_reduction, label = TRUE) + 
        ggtitle("Integrated UMAP - Clusters")
    }
  } else {
    # Fallback
    plot.new()
    text(0.5, 0.5, "No integrated UMAP available", cex = 1.2)
  }
})

# Merge button event
observeEvent(input$merge_button, {
  # Validate we have at least 2 samples selected
  if(is.null(input$merge_files_selector) || length(input$merge_files_selector) < 2) {
    showNotification("Please select at least two samples to merge.", type = "error")
    return()
  }
  
  withProgress(message = "Merging samples", detail = "Loading...", value = 0, {
    rv$merge_log <- "Starting merge process...\n"
    rv$merge_complete <- FALSE
    rv$integrated_reduction <- ""
    
    # Load all the selected Seurat objects from the chosen folder
    seurat_objects <- list()
    merge_folder <- input$merge_folder_selector
    
    incProgress(0.1, detail = "Loading samples...")
    
    # Load each file
    for(i in seq_along(input$merge_files_selector)) {
      file_name <- input$merge_files_selector[i]
      rv$merge_log <- paste0(rv$merge_log, "Loading ", file_name, "...\n")
      
      seurat_obj <- load_sample(file_name, merge_folder)
      if(!is.null(seurat_obj)) {
        seurat_objects[[file_name]] <- seurat_obj
        rv$merge_log <- paste0(rv$merge_log, "Successfully loaded ", file_name, "\n")
      } else {
        rv$merge_log <- paste0(rv$merge_log, "Failed to load ", file_name, "\n")
        showNotification(paste("Failed to load", file_name), type = "error")
      }
    }
    
    # Check if we have at least 2 objects to merge
    if(length(seurat_objects) < 2) {
      rv$merge_log <- paste0(rv$merge_log, "Error: Need at least 2 valid Seurat objects to merge.\n")
      showNotification("Failed to load enough valid Seurat objects.", type = "error")
      return()
    }
    
    # Merge the objects
    incProgress(0.3, detail = "Merging samples...")
    rv$merge_log <- paste0(rv$merge_log, "Merging ", length(seurat_objects), " objects...\n")
    
    tryCatch({
      # Create a vector of sample names
      add_cell_ids <- if(input$add_sample_id_to_cell_names) names(seurat_objects) else NULL
      
      # Merge the objects
      merged_obj <- merge(x = seurat_objects[[1]], 
                         y = seurat_objects[-1], 
                         add.cell.ids = add_cell_ids)
      
      rv$merge_log <- paste0(rv$merge_log, "Successfully merged into one object with ", 
                            ncol(merged_obj), " cells and ", nrow(merged_obj), " features.\n")
      
      incProgress(0.5, detail = "Processing merged data...")
      
      # Run standard preprocessing
      rv$merge_log <- paste0(rv$merge_log, "Processing merged data (normalize, find variable features, scale)...\n")
      
      merged_obj <- NormalizeData(merged_obj)
      merged_obj <- FindVariableFeatures(merged_obj)
      merged_obj <- ScaleData(merged_obj)
      merged_obj <- RunPCA(merged_obj)
      
      # If we're doing integration
      integration_method <- input$integration_method
      dims <- input$integration_dims
      
      if(integration_method != "none") {
        incProgress(0.7, detail = paste("Running", integration_method, "..."))
        rv$merge_log <- paste0(rv$merge_log, "Running integration with method: ", integration_method, "\n")
        
        # Create a batch variable if not exists
        if(!("batch" %in% colnames(merged_obj@meta.data))) {
          if(input$add_sample_id_to_cell_names) {
            # Extract batch from cell names if we added cell IDs
            merged_obj$batch <- sapply(strsplit(colnames(merged_obj), "_"), `[`, 1)
          } else {
            # Otherwise use orig.ident
            merged_obj$batch <- merged_obj$orig.ident
          }
        }
        
        # Set new reduction name
        new_reduction <- paste0("integrated.", tolower(gsub("Integration", "", integration_method)))
        
        # Run integration based on selected method
        tryCatch({
          if(integration_method == "CCAIntegration") {
            merged_obj <- IntegrateLayers(
              object = merged_obj, method = CCAIntegration,
              orig.reduction = "pca", new.reduction = new_reduction,
              verbose = FALSE
            )
          } else if(integration_method == "RPCAIntegration") {
            merged_obj <- IntegrateLayers(
              object = merged_obj, method = RPCAIntegration,
              orig.reduction = "pca", new.reduction = new_reduction,
              verbose = FALSE
            )
          } else if(integration_method == "HarmonyIntegration") {
            merged_obj <- IntegrateLayers(
              object = merged_obj, method = HarmonyIntegration,
              orig.reduction = "pca", new.reduction = new_reduction,
              verbose = FALSE
            )
          } else if(integration_method == "FastMNNIntegration") {
            merged_obj <- IntegrateLayers(
              object = merged_obj, method = FastMNNIntegration,
              new.reduction = new_reduction,
              verbose = FALSE
            )
          }
          
          rv$integrated_reduction <- new_reduction
          rv$merge_log <- paste0(rv$merge_log, "Integration successful. Created reduction: ", 
                                new_reduction, "\n")
        }, error = function(e) {
          rv$merge_log <- paste0(rv$merge_log, "Error in integration: ", e$message, "\n")
          rv$integrated_reduction <- "pca"  # Fallback to PCA
          showNotification(paste("Integration error:", e$message), type = "error")
        })
        
        # Run UMAP if selected
        if(input$run_umap_after_integration) {
          incProgress(0.8, detail = "Running UMAP...")
          rv$merge_log <- paste0(rv$merge_log, "Running UMAP on integrated data...\n")
          
          reduction_to_use <- ifelse(rv$integrated_reduction != "", 
                                    rv$integrated_reduction, "pca")
          
          merged_obj <- RunUMAP(merged_obj, 
                               reduction = reduction_to_use, 
                               dims = 1:dims, 
                               reduction.name = paste0("umap.", 
                                                     gsub("integrated.", "", reduction_to_use)))
          
          rv$merge_log <- paste0(rv$merge_log, "UMAP completed.\n")
        }
        
        # Find clusters if selected
        if(input$find_clusters_after_integration) {
          incProgress(0.9, detail = "Finding clusters...")
          rv$merge_log <- paste0(rv$merge_log, "Finding clusters...\n")
          
          reduction_to_use <- ifelse(rv$integrated_reduction != "", 
                                    rv$integrated_reduction, "pca")
          
          merged_obj <- FindNeighbors(merged_obj, reduction = reduction_to_use, dims = 1:dims)
          merged_obj <- FindClusters(merged_obj, resolution = input$cluster_resolution)
          
          rv$merge_log <- paste0(rv$merge_log, "Clustering completed. Found ", 
                                length(unique(Idents(merged_obj))), " clusters.\n")
        }
      } else {
        # If no integration, just run UMAP and clustering on PCA
        if(input$run_umap_after_integration) {
          incProgress(0.8, detail = "Running UMAP...")
          merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:dims)
        }
        
        if(input$find_clusters_after_integration) {
          incProgress(0.9, detail = "Finding clusters...")
          merged_obj <- FindNeighbors(merged_obj, reduction = "pca", dims = 1:dims)
          merged_obj <- FindClusters(merged_obj, resolution = input$cluster_resolution)
        }
        
        # No integration, so use regular PCA/UMAP
        rv$integrated_reduction <- ""
      }
      
      # Update our sample
      rv$sample <- merged_obj
      
      # Set metadata columns
      rv$metadata_columns <- c("ident", colnames(merged_obj@meta.data))
      if("metadata_column_selector" %in% names(input)) {
        updateSelectInput(session, "metadata_column_selector", 
                        choices = rv$metadata_columns)
      }
      
      # Update cluster names
      rv$cluster_names <- levels(Idents(merged_obj))
      
      # Update cluster selectors
      updatePickerInput(session, "cluster_name_selector", choices = rv$cluster_names)
      updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
      updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
      updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
      updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
      
      # Update reduction selector if present
      if("reduction_selector" %in% names(input)) {
        reductions_list <- Reductions(merged_obj)
        updateSelectInput(session, "reduction_selector", choices = reductions_list)
        
        # Prefer umap if it's available, otherwise use the first reduction
        if("umap" %in% reductions_list) {
          rv$selected_reduction <- "umap"
          updateSelectInput(session, "reduction_selector", choices = reductions_list, selected = "umap")
        } else if(length(reductions_list) > 0) {
          rv$selected_reduction <- reductions_list[1]
        }
      }
      
      # Set sample name
      rv$sample_name <- paste0(
        "Merged_", length(seurat_objects), "_samples_", 
        format(Sys.time(), "%Y%m%d_%H%M%S")
      )
      
      rv$merge_log <- paste0(rv$merge_log, "Merge completed successfully.\n")
      rv$merge_complete <- TRUE
      
      showNotification("Samples merged successfully!", type = "message")
    }, error = function(e) {
      rv$merge_log <- paste0(rv$merge_log, "Error in merge process: ", e$message, "\n")
      showNotification(paste("Error in merge process:", e$message), type = "error")
    })
  })
})