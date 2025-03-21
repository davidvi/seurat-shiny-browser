# server

server <- function(input, output, session) {
  rv <- reactiveValues(
    sample = sample,
    sample_name = sample_name,
    cluster_names = cluster_names,
    gene = "CSF1R",
    all_rds_files = all_rds_files,
    current_folder = ".",
    all_folders = all_folders,
    markers = data.frame(),
    selected_reduction = if(length(available_reductions) > 0) available_reductions[1] else "umap",
    metadata_columns = if(!is.null(sample)) colnames(sample@meta.data) else c("ident"),
    merge_complete = FALSE,
    merge_log = "",
    integrated_reduction = ""
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
        updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
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

  ### search gene ###

  observeEvent(input$search_button, {
    rv$gene <- toupper(input$gene)
  })
  
  ### select reduction ###
  
  observeEvent(input$reduction_selector, {
    rv$selected_reduction <- input$reduction_selector
  })
  
  ### select metadata column ###
  
  observeEvent(input$metadata_column_selector, {
    # Use either the selected metadata column or "ident" (default clustering)
    if(input$metadata_column_selector == "ident") {
      # No need to do anything special, DimPlot uses ident by default
    }
  })
  
  ### select folder ###
  
  observeEvent(input$folder_selector, {
    # Update the current folder
    current_folder <- input$folder_selector
    rv$current_folder <- current_folder
    
    # Get files from the selected folder
    files <- get_rds_files(current_folder)
    rv$all_rds_files <- files
    
    # Update the sample selector
    updateSelectInput(session, "sample", choices = files)
    
    message("Updated file list for folder: ", current_folder, ", found ", length(files), " files")
  })
  
  ### select merge folder ###
  
  observeEvent(input$merge_folder_selector, {
    # Get files from the selected merge folder
    merge_folder <- input$merge_folder_selector
    merge_files <- get_rds_files(merge_folder)
    
    # Update the merge files selector
    updatePickerInput(session, "merge_files_selector", choices = merge_files)
    
    message("Updated merge file list for folder: ", merge_folder, ", found ", length(merge_files), " files")
  })

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

  ### Load sample ###

  observeEvent(input$load_button, {
    withProgress(message = "Loading Sample", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Loading Sample")
      
      message("DEBUG: Loading sample: ", input$sample, " from folder: ", rv$current_folder)
      
      # Load the sample with the current folder
      rv$sample <- load_sample(input$sample, rv$current_folder)
      
      incProgress(0.5, "Set sample name")
      rv$sample_name <- paste0(rv$current_folder, "/", input$sample)
      
      # Set cluster names if sample loaded successfully
      if (!is.null(rv$sample)) {
        rv$cluster_names <- levels(Idents(rv$sample))
      } else {
        rv$cluster_names <- c()
      }
      
      rv$markers <- NULL
      
      # Update available reductions after loading new sample
      tryCatch({
        reductions_list <- Reductions(rv$sample)
        if(length(reductions_list) > 0) {
          updateSelectInput(session, "reduction_selector", choices = reductions_list)
          rv$selected_reduction <- reductions_list[1]
        }
      }, error = function(e) {
        # Use default if error
        updateSelectInput(session, "reduction_selector", choices = c("umap"))
        rv$selected_reduction <- "umap"
      })
      
      # Update available metadata columns
      tryCatch({
        if(!is.null(rv$sample)) {
          metadata_cols <- colnames(rv$sample@meta.data)
          rv$metadata_columns <- c("ident", metadata_cols)
          updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
        }
      }, error = function(e) {
        # Use default if error
        rv$metadata_columns <- c("ident")
        updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
      })
      incProgress(0.9, "Done")
    })
  })

  ### save as new ###

  # Function to refresh the list of files in the current folder
  refresh_folder_files <- function() {
    # Get files from the current folder
    files <- get_rds_files(rv$current_folder)
    rv$all_rds_files <- files
    
    # Update the sample selector
    updateSelectInput(session, "sample", choices = files)
    
    message("Refreshed file list for folder: ", rv$current_folder, ", now contains ", length(files), " files")
  }
  
  observeEvent(input$save_button, {
    # Get the folder and determine file name
    folder <- rv$current_folder
    
    # If the text box is empty, use the old name (without the folder prefix)
    # Otherwise, use the entered name
    if (input$new_save_name == "") {
      # Extract just the filename part from the sample_name (which includes folder path)
      current_name <- basename(rv$sample_name)
      file_name <- current_name
    } else {
      file_name <- input$new_save_name
    }
    
    # Make sure the file name ends with .rds
    if (!endsWith(file_name, ".rds")) {
      file_name <- paste0(file_name, ".rds")
    }
    
    withProgress(message = "Saving Sample", detail = "Please wait...", value = 0, {
      incProgress(0.1, "Saving Sample")
      
      # Determine the save location
      folder_path <- paste0(base_folder, folder, "/")
      # Create the directory if it doesn't exist
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      save_path <- paste0(folder_path, file_name)
      
      message("DEBUG: Saving to ", save_path)
      saveRDS(rv$sample, file = save_path)
      
      incProgress(0.8, "Setting new name")
      
      # Update the sample name with folder prefix
      rv$sample_name <- paste0(folder, "/", file_name)
      
      # Refresh the file list
      refresh_folder_files()
      incProgress(0.9, "Done")
      
      # Display success notification
      showNotification(paste("File saved as", file_name), type = "message")
    })
  })

  # Update UI elements based on reactive values
  observe({
    # Update cluster selectors with current cluster names
    updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
    updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
    updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
    updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
    
    # Ensure metadata column selector is updated
    if(!is.null(rv$sample)) {
      metadata_cols <- colnames(rv$sample@meta.data)
      rv$metadata_columns <- c("ident", metadata_cols)
      updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
    }
  })

  # Outputs for all tabs
  output$sample_name <- renderText({
    paste("Loaded Sample: ", rv$sample_name)
  })
  
  # Merge tab functionality
  
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
        updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
        updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
        updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
        updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
        
        # Update reduction selector if present
        if("reduction_selector" %in% names(input)) {
          reductions_list <- Reductions(merged_obj)
          updateSelectInput(session, "reduction_selector", choices = reductions_list)
          if(length(reductions_list) > 0) {
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
  
  output$markers_table <- renderDataTable({
    rv$markers
  })
}