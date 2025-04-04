# Server functions for the auto name tab

# Main function for the auto name tab server logic
autoNameServer <- function(input, output, session, values, singleR_values) {
  
  # Function to get reference data based on user selection
  get_reference_data <- function(ref_choice) {
    # Show waiting notification
    showNotification(
      "Loading reference dataset...", 
      type = "message", 
      duration = NULL, 
      id = "ref_loading"
    )
    
    # Try to load reference data
    result <- tryCatch({
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
      ref_data
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
    
    # Get the reference dataset
    ref_choice <- input$reference_dataset
    label_type <- input$label_type
    
    # Try to run SingleR in a separate process to avoid blocking the UI
    withProgress(message = 'Running SingleR annotation...', value = 0.1, {
      tryCatch({
        # Check if SingleR package is available
        if (!requireNamespace("SingleR", quietly = TRUE)) {
          stop("The 'SingleR' package is required. Please install it with BiocManager::install('SingleR')")
        }
        
        # Get reference data
        ref_data <- get_reference_data(ref_choice)
        req(ref_data)
        
        # Determine label column based on user selection
        label_column <- if (label_type == "main") "label.main" else "label.fine"
        labels <- ref_data[[label_column]]
        
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
        
        # Run SingleR on cluster averages
        singler_results <- SingleR::SingleR(
          test = cluster_data,
          ref = ref_subset,
          labels = labels,
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
  
  # Transfer labels to Idents
  observeEvent(input$transfer_labels, {
    req(values$sample, singleR_values$completed, singleR_values$cluster_mapping)
    
    tryCatch({
      # Get the Seurat object
      seurat_obj <- values$sample
      
      # Extract the mapping
      mapping <- singleR_values$cluster_mapping
      
      # First, add the SingleR labels as a new column in metadata
      # Get the current cluster identities
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
      
      # Show success notification
      showNotification(
        "SingleR labels transferred to active identities",
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
  
  # No longer needed - removed download functionality
}