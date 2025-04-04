# Server functions for the auto name tab

# Main function for the auto name tab server logic
autoNameServer <- function(input, output, session, values, singleR_values) {
  
  # Initialize the folder dropdown with all available folders
  observe({
    req(values)
    # Use the same all_folders that all other components use
    updateSelectInput(session, "reference_folder", 
                      choices = c("", all_folders),
                      selected = input$reference_folder)
  })
  
  # Update sample list when folder is selected - using global get_rds_files function
  observeEvent(input$reference_folder, {
    req(input$reference_folder, input$reference_folder != "")
    
    # Use the same file loading logic used by the rest of the app
    ref_files <- get_rds_files(input$reference_folder)
    
    # Update the sample selector
    if (length(ref_files) == 0) {
      updateSelectInput(session, "reference_sample",
                        choices = c("No RDS files found" = ""),
                        selected = "")
      message("No RDS files found in folder:", input$reference_folder)
    } else {
      updateSelectInput(session, "reference_sample",
                        choices = c("", ref_files),
                        selected = "")
      message("Found ", length(ref_files), " RDS files in folder: ", input$reference_folder)
    }
  })
  
  # Update metadata columns when sample is selected
  observe({
    req(input$reference_folder, input$reference_folder != "", 
        input$reference_sample, input$reference_sample != "")
    
    tryCatch({
      # Use the same loading function used in other parts of the app
      file_path <- paste0(base_folder, input$reference_folder, "/", input$reference_sample)
      message("Loading metadata for sample: ", file_path)
      
      # Check if file exists
      if (!file.exists(file_path)) {
        message("Reference file not found: ", file_path)
        updateSelectInput(session, "reference_column",
                          choices = c("Reference file not found" = ""),
                          selected = "")
        return(NULL)
      }
      
      # Load the object
      ref_obj <- readRDS(file_path)
      
      # Get metadata columns
      if(!inherits(ref_obj, "Seurat")) {
        message("File is not a Seurat object")
        updateSelectInput(session, "reference_column",
                          choices = c("Not a Seurat object" = ""),
                          selected = "")
        return(NULL)
      }
      
      meta_cols <- colnames(ref_obj@meta.data)
      message("Found metadata columns: ", paste(meta_cols, collapse = ", "))
      
      # Find columns that might contain cell type labels
      likely_label_cols <- meta_cols[
        grepl("label|type|cell|cluster|ident|annotation", meta_cols, ignore.case = TRUE)
      ]
      
      # If no likely columns found, use all columns
      if(length(likely_label_cols) == 0) {
        likely_label_cols <- meta_cols
      }
      
      message("Likely label columns: ", paste(likely_label_cols, collapse = ", "))
      
      # Add active identity as an option
      updateSelectInput(session, "reference_column",
                        choices = c("active.ident", likely_label_cols),
                        selected = "active.ident")
      
    }, error = function(e) {
      message("Error loading reference sample metadata: ", e$message)
      showNotification(
        paste("Error loading reference sample metadata:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Function to get reference data based on user selection
  get_reference_data <- function(ref_choice, ref_type) {
    # Show waiting notification
    showNotification(
      "Loading reference dataset...", 
      type = "message", 
      duration = NULL, 
      id = "ref_loading"
    )
    
    # Try to load reference data
    result <- tryCatch({
      if (ref_type == "builtin") {
        # Using built-in reference
        if (!requireNamespace("celldex", quietly = TRUE)) {
          stop("The 'celldex' package is required. Please install it with BiocManager::install('celldex')")
        }
        
        if (ref_choice == "hpca") {
          ref_data <- celldex::HumanPrimaryCellAtlasData()
        } else if (ref_choice == "blueprint_encode") {
          ref_data <- celldex::BlueprintEncodeData()
        } else if (ref_choice == "mouse_rnaseq") {
          ref_data <- celldex::MouseRNAseqData()
        } else {
          stop("Invalid reference dataset selected")
        }
        
        # Return as a list with the reference data and labels
        list(
          data = ref_data,
          labels = NULL,  # Will be determined later based on label_type
          type = "builtin"
        )
        
      } else if (ref_type == "sample") {
        # Using another Seurat object as reference
        req(input$reference_folder, input$reference_folder != "",
            input$reference_sample, input$reference_sample != "",
            input$reference_column, input$reference_column != "")
        
        # Use the same file path logic as the rest of the app
        file_path <- paste0(base_folder, input$reference_folder, "/", input$reference_sample)
        
        # Check if file exists
        if (!file.exists(file_path)) {
          stop("Reference file not found: ", file_path)
        }
        
        # Load the reference Seurat object
        message("Loading reference Seurat object from: ", file_path)
        ref_obj <- readRDS(file_path)
        
        if(!inherits(ref_obj, "Seurat")) {
          stop("Reference file is not a Seurat object")
        }
        
        # Get the labels from the reference object
        if(input$reference_column == "active.ident") {
          labels <- as.character(ref_obj@active.ident)
          
          # Check if we have valid labels
          if(length(labels) != ncol(ref_obj)) {
            stop("Active identities don't match the number of cells in the reference")
          }
        } else {
          # Get labels from metadata column
          if(!input$reference_column %in% colnames(ref_obj@meta.data)) {
            stop(paste("Column", input$reference_column, "not found in reference object metadata"))
          }
          
          labels <- as.character(ref_obj@meta.data[[input$reference_column]])
          
          # Check for NA or empty labels
          if(any(is.na(labels)) || any(labels == "")) {
            message("Warning: Some labels in the reference are NA or empty")
            # Filter out cells with NA or empty labels
            valid_cells <- !is.na(labels) & labels != ""
            ref_obj <- ref_obj[, valid_cells]
            labels <- labels[valid_cells]
          }
        }
        
        # Check if we have any labels
        if(length(unique(labels)) <= 1) {
          stop("Reference object doesn't have enough distinct labels")
        }
        
        message("Reference object has ", length(unique(labels)), " distinct labels")
        
        # Extract normalized data
        norm_data <- tryCatch({
          # Try Seurat 5 style
          if (inherits(ref_obj[["RNA"]], "Assay5")) {
            Seurat::GetAssayData(ref_obj, assay = "RNA", layer = "data")
          } else {
            # For Seurat v3
            ref_obj@assays$RNA@data
          }
        }, error = function(e) {
          # Generic method as fallback
          Seurat::GetAssayData(ref_obj, slot = "data")
        })
        
        # Return as a list with the reference data and labels
        list(
          data = norm_data,
          labels = labels,
          type = "sample"
        )
      } else {
        stop("Invalid reference type")
      }
    }, error = function(e) {
      removeNotification(id = "ref_loading")
      showNotification(
        paste("Error loading reference data:", e$message),
        type = "error",
        duration = 10
      )
      return(NULL)
    })
    
    # Remove loading notification
    removeNotification(id = "ref_loading")
    
    if (!is.null(result)) {
      showNotification(
        "Reference dataset loaded successfully",
        type = "message",
        duration = 3
      )
    }
    
    return(result)
  }
  
  # Run SingleR when the button is clicked
  observeEvent(input$run_singleR, {
    req(values$sample)
    
    # Set status to running
    singleR_values$running <- TRUE
    singleR_values$completed <- FALSE
    singleR_values$error <- NULL
    
    # Get reference type and details
    ref_type <- input$reference_type
    ref_choice <- if(ref_type == "builtin") input$reference_dataset else NULL
    
    # Try to run SingleR in a separate process to avoid blocking the UI
    withProgress(message = 'Running SingleR annotation...', value = 0.1, {
      tryCatch({
        # Check if SingleR package is available
        if (!requireNamespace("SingleR", quietly = TRUE)) {
          stop("The 'SingleR' package is required. Please install it with BiocManager::install('SingleR')")
        }
        
        # Get reference data
        ref_result <- get_reference_data(ref_choice, ref_type)
        req(ref_result)
        
        # Get reference data and labels based on reference type
        if (ref_result$type == "builtin") {
          ref_data <- ref_result$data
          
          # Determine label column based on user selection
          label_type <- input$label_type
          label_column <- if (label_type == "main") "label.main" else "label.fine"
          labels <- ref_data[[label_column]]
        } else {
          # For sample reference, use the normalized data and provided labels
          ref_data <- ref_result$data
          labels <- ref_result$labels
        }
        
        # Simplify SingleR processing to work at the cluster level
        # First, get cluster identities
        clusters <- values$sample@meta.data$seurat_clusters
        if (is.null(clusters)) {
          clusters <- values$sample@active.ident
        }
        
        # Ensure clusters are properly formatted
        clusters <- as.character(clusters)
        cluster_levels <- unique(clusters)
        
        # Extract normalized expression data
        # First try standard methods to get normalized data
        norm_data <- tryCatch({
          # Try Seurat 5 style
          if (inherits(values$sample[["RNA"]], "Assay5")) {
            Seurat::GetAssayData(values$sample, assay = "RNA", layer = "data")
          } else {
            # For Seurat v3
            values$sample@assays$RNA@data
          }
        }, error = function(e) {
          # Generic method as fallback
          Seurat::GetAssayData(values$sample, slot = "data")
        })
        
        # Verify the data is valid for processing
        if (nrow(norm_data) == 0) {
          stop("No normalized data found in the Seurat object.")
        }
        
        # Get cluster averages (pseudo-bulk) instead of processing all cells
        # This is more efficient and avoids memory issues
        cluster_averages <- matrix(0, nrow = nrow(norm_data), ncol = length(cluster_levels))
        rownames(cluster_averages) <- rownames(norm_data)
        colnames(cluster_averages) <- cluster_levels
        
        # Calculate average expression per cluster
        for (cluster in cluster_levels) {
          cells_in_cluster <- colnames(norm_data)[clusters == cluster]
          if (length(cells_in_cluster) > 0) {
            cluster_averages[, cluster] <- rowMeans(norm_data[, cells_in_cluster, drop = FALSE])
          }
        }
        
        # Run SingleR on cluster averages
        incProgress(0.3, detail = "Running annotation...")
        
        # Get common genes between reference and test data
        common_genes <- intersect(rownames(cluster_averages), rownames(ref_data))
        if (length(common_genes) < 200) {
          stop(paste0("Only ", length(common_genes), " genes in common between datasets. At least 200 are needed."))
        }
        
        # Subset matrices to common genes
        cluster_data <- cluster_averages[common_genes, , drop = FALSE]
        ref_subset <- ref_data[common_genes, ]
        
        # Determine the differential expression method
        de_method <- input$de_method
        
        # Run SingleR on cluster averages with appropriate method
        singler_results <- SingleR::SingleR(
          test = cluster_data,
          ref = ref_subset,
          labels = labels,
          de.method = if(de_method == "wilcox") "wilcox" else NULL,
          assay.type.test = 1,  # Use the first element (logcounts) for cluster averages
          BPPARAM = BiocParallel::SerialParam() # Force serial to avoid parallel errors
        )
        
        # Create mapping of cluster to cell type
        cluster_mapping <- data.frame(
          cluster = cluster_levels,
          cell_type = singler_results$labels,
          stringsAsFactors = FALSE
        )
        
        # Store results
        singleR_values$results <- singler_results
        singleR_values$cluster_mapping <- cluster_mapping
        singleR_values$completed <- TRUE
        
        # Show success notification
        showNotification(
          "SingleR annotation completed successfully",
          type = "message",
          duration = 5
        )
        
      }, error = function(e) {
        # Handle errors
        singleR_values$error <- e$message
        singleR_values$completed <- FALSE
        singleR_values$running <- FALSE
        
        showNotification(
          paste("Error in SingleR annotation:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Set status to not running after completion
    singleR_values$running <- FALSE
  })
  
  # Helper function to get cluster column and prepare mapping
  get_cluster_mapping <- function(seurat_obj, mapping) {
    # Get the current cluster identities column
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      cluster_column <- "seurat_clusters"
    } else {
      # If no seurat_clusters column, create one from active identities
      seurat_obj$temp_clusters <- as.character(seurat_obj@active.ident)
      cluster_column <- "temp_clusters"
    }
    
    # Create a lookup table for the mapping
    cluster_to_celltype <- setNames(
      mapping$cell_type,
      mapping$cluster
    )
    
    return(list(
      seurat_obj = seurat_obj,
      cluster_column = cluster_column,
      cluster_to_celltype = cluster_to_celltype
    ))
  }
  
  # Transfer labels to Idents (cell type names only)
  observeEvent(input$transfer_labels, {
    req(values$sample, singleR_values$completed, singleR_values$cluster_mapping)
    
    tryCatch({
      # Get the Seurat object
      seurat_obj <- values$sample
      
      # Extract the mapping
      mapping <- singleR_values$cluster_mapping
      
      # Get cluster mapping info
      mapping_info <- get_cluster_mapping(seurat_obj, mapping)
      seurat_obj <- mapping_info$seurat_obj
      cluster_column <- mapping_info$cluster_column
      cluster_to_celltype <- mapping_info$cluster_to_celltype
      
      # Add a new metadata column with the cell type labels
      cell_types <- character(length = nrow(seurat_obj@meta.data))
      for (i in seq_along(cell_types)) {
        cluster_id <- as.character(seurat_obj@meta.data[[cluster_column]][i])
        if (cluster_id %in% names(cluster_to_celltype)) {
          cell_types[i] <- cluster_to_celltype[[cluster_id]]
        } else {
          cell_types[i] <- "Unknown"
        }
      }
      
      # Add the cell types to metadata
      seurat_obj@meta.data$singleR_labels <- cell_types
      
      # Set as active ident
      seurat_obj <- Seurat::SetIdent(seurat_obj, value = "singleR_labels")
      
      # Update the Seurat object
      values$sample <- seurat_obj
      
      # Update cluster names after changing identities
      values$cluster_names <- levels(Idents(seurat_obj))
      
      # Update all cluster selectors
      updatePickerInput(session, "cluster_name_selector", choices = values$cluster_names)
      updatePickerInput(session, "delete_cluster_selector", choices = values$cluster_names)
      updateSelectInput(session, "cluster_name_selector_analyze", choices = values$cluster_names)
      updateSelectInput(session, "first_cluster_name_selector", choices = values$cluster_names)
      updateSelectInput(session, "second_cluster_name_selector", choices = values$cluster_names)
      
      # Show success notification
      showNotification(
        "Cell type labels transferred to active identities",
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("Error transferring labels:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Transfer labels with cluster prefix to Idents
  observeEvent(input$transfer_prefix_labels, {
    req(values$sample, singleR_values$completed, singleR_values$cluster_mapping)
    
    tryCatch({
      # Get the Seurat object
      seurat_obj <- values$sample
      
      # Extract the mapping
      mapping <- singleR_values$cluster_mapping
      
      # Get cluster mapping info
      mapping_info <- get_cluster_mapping(seurat_obj, mapping)
      seurat_obj <- mapping_info$seurat_obj
      cluster_column <- mapping_info$cluster_column
      cluster_to_celltype <- mapping_info$cluster_to_celltype
      
      # Add a new metadata column with prefixed labels (e.g., "0_B-cells")
      prefixed_labels <- character(length = nrow(seurat_obj@meta.data))
      for (i in seq_along(prefixed_labels)) {
        cluster_id <- as.character(seurat_obj@meta.data[[cluster_column]][i])
        if (cluster_id %in% names(cluster_to_celltype)) {
          cell_type <- cluster_to_celltype[[cluster_id]]
          prefixed_labels[i] <- paste0(cluster_id, "_", cell_type)
        } else {
          prefixed_labels[i] <- paste0(cluster_id, "_Unknown")
        }
      }
      
      # Add the prefixed labels to metadata
      seurat_obj@meta.data$singleR_prefixed_labels <- prefixed_labels
      
      # Set as active ident
      seurat_obj <- Seurat::SetIdent(seurat_obj, value = "singleR_prefixed_labels")
      
      # Update the Seurat object
      values$sample <- seurat_obj
      
      # Update cluster names after changing identities
      values$cluster_names <- levels(Idents(seurat_obj))
      
      # Update all cluster selectors
      updatePickerInput(session, "cluster_name_selector", choices = values$cluster_names)
      updatePickerInput(session, "delete_cluster_selector", choices = values$cluster_names)
      updateSelectInput(session, "cluster_name_selector_analyze", choices = values$cluster_names)
      updateSelectInput(session, "first_cluster_name_selector", choices = values$cluster_names)
      updateSelectInput(session, "second_cluster_name_selector", choices = values$cluster_names)
      
      # Show success notification
      showNotification(
        "Cluster+Cell type labels transferred to active identities",
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("Error transferring prefixed labels:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Reset to original clusters
  observeEvent(input$reset_clusters, {
    req(values$sample)
    
    tryCatch({
      # Get the Seurat object
      seurat_obj <- values$sample
      
      # Check if seurat_clusters exists
      if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        # Set back to original clusters
        seurat_obj <- Seurat::SetIdent(seurat_obj, value = "seurat_clusters")
        
        # Update the Seurat object
        values$sample <- seurat_obj
        
        # Update cluster names after changing identities
        values$cluster_names <- levels(Idents(seurat_obj))
        
        # Update all cluster selectors
        updatePickerInput(session, "cluster_name_selector", choices = values$cluster_names)
        updatePickerInput(session, "delete_cluster_selector", choices = values$cluster_names)
        updateSelectInput(session, "cluster_name_selector_analyze", choices = values$cluster_names)
        updateSelectInput(session, "first_cluster_name_selector", choices = values$cluster_names)
        updateSelectInput(session, "second_cluster_name_selector", choices = values$cluster_names)
        
        # Show success notification
        showNotification(
          "Reset to original clusters",
          type = "message",
          duration = 5
        )
      } else {
        # Warning if no seurat_clusters
        showNotification(
          "No 'seurat_clusters' column found in metadata",
          type = "warning",
          duration = 5
        )
      }
      
    }, error = function(e) {
      showNotification(
        paste("Error resetting clusters:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Output for showing the running status
  output$singleR_running <- reactive({
    singleR_values$running
  })
  outputOptions(output, "singleR_running", suspendWhenHidden = FALSE)
  
  # Output for showing completion status
  output$singleR_completed <- reactive({
    singleR_values$completed
  })
  outputOptions(output, "singleR_completed", suspendWhenHidden = FALSE)
  
  # Output for showing error status
  output$singleR_error <- reactive({
    !is.null(singleR_values$error)
  })
  outputOptions(output, "singleR_error", suspendWhenHidden = FALSE)
  
  # Output for error text
  output$error_text <- renderText({
    singleR_values$error
  })
  
  # Render the annotation table
  output$annotation_table <- DT::renderDataTable({
    req(singleR_values$completed, singleR_values$cluster_mapping)
    
    mapping <- singleR_values$cluster_mapping
    colnames(mapping) <- c("Cluster", "Cell Type")
    
    DT::datatable(
      mapping,
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        searching = TRUE
      ),
      rownames = FALSE
    )
  })
}