# Clustering analysis module

# Add a reactive value to track the completion of each step
rv$clustering_progress <- list(
  qc_done = FALSE,
  normalization_done = FALSE,
  variable_features_done = FALSE,
  scaling_done = FALSE,
  pca_done = FALSE,
  neighbors_done = FALSE,
  clustering_done = FALSE,
  umap_done = FALSE
)

# Helper function to update completion status
update_step_status <- function(step, status) {
  rv$clustering_progress[[step]] <- status
}

# QC metrics plots
output$qc_violin_plot <- renderPlot({
  req(rv$sample)
  
  # Check if we have mt percent calculated already
  if(!"percent.mt" %in% colnames(rv$sample@meta.data)) {
    # Add mt percent if not present and if mt genes are detected
    mt_pattern <- "^MT-"
    if(any(grepl(mt_pattern, rownames(rv$sample)))) {
      rv$sample[["percent.mt"]] <- PercentageFeatureSet(rv$sample, pattern = mt_pattern)
    }
  }
  
  # Create violin plots
  VlnPlot(rv$sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3, pt.size = 0) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

output$feature_scatter_plot <- renderPlot({
  req(rv$sample)
  
  # Create scatter plots
  p1 <- FeatureScatter(rv$sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle("Features vs Counts")
  
  if("percent.mt" %in% colnames(rv$sample@meta.data)) {
    p2 <- FeatureScatter(rv$sample, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      ggtitle("Mitochondrial % vs Counts")
    
    p1 + p2
  } else {
    p1
  }
})

# Variable features plot
output$var_features_plot <- renderPlot({
  req(rv$sample)
  
  # Check if variable features has been run
  if(length(VariableFeatures(rv$sample)) > 0) {
    # Plot variable features
    LabelPoints(VariableFeaturePlot(rv$sample), 
                points = head(VariableFeatures(rv$sample), 10), 
                repel = TRUE)
  } else {
    # Show message that variable features haven't been calculated
    plot.new()
    text(0.5, 0.5, "Variable features haven't been calculated yet.\nRun Step 3 to view this plot.", cex = 1.2)
  }
})

# PCA plots
output$elbow_plot <- renderPlot({
  req(rv$sample)
  
  # Check if PCA has been run
  if("pca" %in% Reductions(rv$sample)) {
    ElbowPlot(rv$sample, ndims = min(50, length(rv$sample@reductions$pca@stdev)))
  } else {
    # Show message that PCA hasn't been run
    plot.new()
    text(0.5, 0.5, "PCA hasn't been calculated yet.\nRun Step 5 to view this plot.", cex = 1.2)
  }
})

output$pca_heatmap <- renderPlot({
  req(rv$sample)
  
  # Check if PCA has been run
  if("pca" %in% Reductions(rv$sample)) {
    dims_to_plot <- min(6, length(rv$sample@reductions$pca@stdev))
    DimHeatmap(rv$sample, dims = 1:dims_to_plot, cells = 500, balanced = TRUE, combine = TRUE)
  } else {
    # Show message that PCA hasn't been run
    plot.new()
    text(0.5, 0.5, "PCA hasn't been calculated yet.\nRun Step 5 to view this plot.", cex = 1.2)
  }
})

# UMAP plot
output$umap_plot <- renderPlot({
  req(rv$sample)
  
  # Check if UMAP has been run
  if("umap" %in% Reductions(rv$sample)) {
    DimPlot(rv$sample, reduction = "umap", label = TRUE) + 
      ggtitle(paste0("UMAP Visualization (", length(unique(Idents(rv$sample))), " clusters)"))
  } else {
    # Show message that UMAP hasn't been run
    plot.new()
    text(0.5, 0.5, "UMAP hasn't been calculated yet.\nRun Step 8 to view this plot.", cex = 1.2)
  }
})

# Cluster statistics
output$cluster_stats <- renderTable({
  req(rv$sample)
  
  # Get cluster IDs
  clusters <- as.character(sort(unique(Idents(rv$sample))))
  
  # If there are clusters, return statistics
  if(length(clusters) > 0) {
    # Count cells per cluster
    cells_per_cluster <- sapply(clusters, function(cluster) {
      sum(Idents(rv$sample) == cluster)
    })
    
    # Create data frame
    cluster_data <- data.frame(
      Cluster = clusters,
      CellCount = cells_per_cluster,
      Percentage = round(100 * cells_per_cluster / sum(cells_per_cluster), 2)
    )
    
    cluster_data
  } else {
    # If no clusters, return empty data frame
    data.frame(
      Cluster = character(0),
      CellCount = integer(0),
      Percentage = numeric(0)
    )
  }
}, rownames = FALSE, striped = TRUE, bordered = TRUE, width = "100%")

# Download UMAP plot
output$download_umap <- downloadHandler(
  filename = function() {
    paste0(rv$sample_name, "_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  },
  content = function(file) {
    png(file, width = 1000, height = 800, res = 120)
    print(DimPlot(rv$sample, reduction = "umap", label = TRUE) + 
            ggtitle(paste0("UMAP Visualization (", length(unique(Idents(rv$sample))), " clusters)")))
    dev.off()
  }
)

# Helper function for updating UI elements
update_cluster_ui <- function() {
  # Update cluster names
  rv$cluster_names <- levels(Idents(rv$sample))
  
  # Update cluster selectors
  updatePickerInput(session, "cluster_name_selector", choices = rv$cluster_names)
  updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
  updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
  updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
  updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
  
  # Update metadata columns
  if(!is.null(rv$sample)) {
    metadata_cols <- colnames(rv$sample@meta.data)
    rv$metadata_columns <- c("ident", metadata_cols)
    updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
  }
  
  # Update reduction selector if present
  if("reduction_selector" %in% names(input)) {
    reductions_list <- Reductions(rv$sample)
    updateSelectInput(session, "reduction_selector", choices = reductions_list)
    
    # Prefer umap if available
    if("umap" %in% reductions_list) {
      rv$selected_reduction <- "umap"
      updateSelectInput(session, "reduction_selector", choices = reductions_list, selected = "umap")
    }
  }
}

# Add mitochondrial percentage if not present
ensure_mt_percent <- function() {
  # Check if we have mt percent calculated already
  if(!"percent.mt" %in% colnames(rv$sample@meta.data)) {
    # Add mt percent if not present and if mt genes are detected
    mt_pattern <- "^MT-"
    if(any(grepl(mt_pattern, rownames(rv$sample)))) {
      rv$sample[["percent.mt"]] <- PercentageFeatureSet(rv$sample, pattern = mt_pattern)
      message("Added percent.mt to metadata")
    } else {
      rv$sample[["percent.mt"]] <- 0
      message("No mitochondrial genes found, setting percent.mt to 0")
    }
  }
}

# Step 1: QC Filtering
observeEvent(input$run_qc, {
  req(rv$sample)
  
  withProgress(message = "Running QC filtering", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Ensure mitochondrial percentage is calculated
      ensure_mt_percent()
      
      # Apply QC filters
      rv$sample <- subset(rv$sample, 
                          subset = nFeature_RNA > input$min_features & 
                            nFeature_RNA < input$max_features & 
                            percent.mt < input$max_mt_percent)
      
      # Check if we still have cells
      if(ncol(rv$sample) == 0) {
        showNotification("QC filtering removed all cells. Try less stringent parameters.", type = "error")
        rv$sample <- original_sample
        return()
      }
      
      # Mark step as complete
      update_step_status("qc_done", TRUE)
      
      # Update UI elements
      update_cluster_ui()
      
      showNotification(
        paste0("QC filtering complete. Kept ", ncol(rv$sample), " cells."),
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error in QC filtering:", e$message), type = "error")
    })
  })
})

# Step 2: Normalization
observeEvent(input$run_normalization, {
  req(rv$sample)
  
  # Check if QC has been done
  if(!rv$clustering_progress$qc_done) {
    showNotification("Please run QC filtering first (Step 1)", type = "warning")
    return()
  }
  
  withProgress(message = "Running normalization", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Run normalization based on the selected method
      if(input$normalization_method == "LogNormalize") {
        rv$sample <- NormalizeData(rv$sample, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = input$scale_factor)
      } else if(input$normalization_method == "SCTransform") {
        rv$sample <- SCTransform(rv$sample)
      }
      
      # Mark step as complete
      update_step_status("normalization_done", TRUE)
      
      # If SCTransform was used, also mark variable features as complete
      if(input$normalization_method == "SCTransform") {
        update_step_status("variable_features_done", TRUE)
      }
      
      showNotification(
        "Normalization complete",
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error in normalization:", e$message), type = "error")
    })
  })
})

# Step 3: Find Variable Features
observeEvent(input$run_variable_features, {
  req(rv$sample)
  
  # Check if normalization has been done
  if(!rv$clustering_progress$normalization_done) {
    showNotification("Please run normalization first (Step 2)", type = "warning")
    return()
  }
  
  # Skip if SCTransform was used as it already finds variable features
  if(input$normalization_method == "SCTransform") {
    showNotification("SCTransform already identified variable features", type = "message")
    return()
  }
  
  withProgress(message = "Finding variable features", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Find variable features
      rv$sample <- FindVariableFeatures(rv$sample, 
                                       selection.method = "vst", 
                                       nfeatures = input$n_features)
      
      # Mark step as complete
      update_step_status("variable_features_done", TRUE)
      
      showNotification(
        paste0("Found ", length(VariableFeatures(rv$sample)), " variable features"),
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error finding variable features:", e$message), type = "error")
    })
  })
})

# Step 4: Scale Data
observeEvent(input$run_scaling, {
  req(rv$sample)
  
  # Check if variable features has been done
  if(!rv$clustering_progress$variable_features_done) {
    showNotification("Please find variable features first (Step 3)", type = "warning")
    return()
  }
  
  withProgress(message = "Scaling data", detail = "This may take a while for large datasets...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Scale data
      if(input$scale_all_genes) {
        # Scale all genes
        all_genes <- rownames(rv$sample)
        rv$sample <- ScaleData(rv$sample, features = all_genes)
      } else {
        # Scale only variable features (default)
        rv$sample <- ScaleData(rv$sample)
      }
      
      # Mark step as complete
      update_step_status("scaling_done", TRUE)
      
      showNotification(
        "Data scaling complete",
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error scaling data:", e$message), type = "error")
    })
  })
})

# Step 5: Run PCA
observeEvent(input$run_pca, {
  req(rv$sample)
  
  # Check if scaling has been done
  if(!rv$clustering_progress$scaling_done) {
    showNotification("Please scale data first (Step 4)", type = "warning")
    return()
  }
  
  withProgress(message = "Running PCA", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Run PCA
      rv$sample <- RunPCA(rv$sample, npcs = input$npcs)
      
      # Mark step as complete
      update_step_status("pca_done", TRUE)
      
      # Update reduction selector
      if("reduction_selector" %in% names(input)) {
        reductions_list <- Reductions(rv$sample)
        updateSelectInput(session, "reduction_selector", choices = reductions_list)
      }
      
      showNotification(
        paste0("PCA complete. Computed ", input$npcs, " principal components."),
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error running PCA:", e$message), type = "error")
    })
  })
})

# Step 6: Find Neighbors
observeEvent(input$run_neighbors, {
  req(rv$sample)
  
  # Check if PCA has been done
  if(!rv$clustering_progress$pca_done) {
    showNotification("Please run PCA first (Step 5)", type = "warning")
    return()
  }
  
  withProgress(message = "Finding neighbors", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Find neighbors
      pc_range <- input$pcs_for_clustering[1]:input$pcs_for_clustering[2]
      rv$sample <- FindNeighbors(rv$sample, dims = pc_range)
      
      # Mark step as complete
      update_step_status("neighbors_done", TRUE)
      
      showNotification(
        "Neighbor finding complete",
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error finding neighbors:", e$message), type = "error")
    })
  })
})

# Step 7: Find Clusters
observeEvent(input$run_clustering, {
  req(rv$sample)
  
  # Check if neighbors has been done
  if(!rv$clustering_progress$neighbors_done) {
    showNotification("Please find neighbors first (Step 6)", type = "warning")
    return()
  }
  
  withProgress(message = "Finding clusters", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Find clusters
      rv$sample <- FindClusters(rv$sample, 
                               resolution = input$resolution,
                               algorithm = as.numeric(input$cluster_algorithm),
                               modularity.fxn = as.numeric(input$modularity_function),
                               n.start = input$cluster_n_start,
                               n.iter = input$cluster_n_iter,
                               group.singletons = input$group_singletons)
      
      # Mark step as complete
      update_step_status("clustering_done", TRUE)
      
      # Update UI elements
      update_cluster_ui()
      
      showNotification(
        paste0("Clustering complete. Found ", length(unique(Idents(rv$sample))), " clusters."),
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error finding clusters:", e$message), type = "error")
    })
  })
})

# Step 8: Run UMAP
observeEvent(input$run_umap, {
  req(rv$sample)
  
  # Check if clustering has been done
  if(!rv$clustering_progress$neighbors_done) {
    showNotification("Please find neighbors first (Step 6)", type = "warning")
    return()
  }
  
  withProgress(message = "Running UMAP", detail = "Please wait...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Run UMAP
      pc_range <- input$pcs_for_clustering[1]:input$pcs_for_clustering[2]
      rv$sample <- RunUMAP(rv$sample, 
                          dims = pc_range,
                          min.dist = input$umap_mindist,
                          spread = input$umap_spread)
      
      # Mark step as complete
      update_step_status("umap_done", TRUE)
      
      # Update UI elements
      update_cluster_ui()
      
      showNotification(
        "UMAP complete",
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error running UMAP:", e$message), type = "error")
    })
  })
})

# Run all steps
observeEvent(input$run_all_steps, {
  req(rv$sample)
  
  withProgress(message = "Running all steps", detail = "This may take several minutes...", value = 0, {
    tryCatch({
      # Make a copy of the original Seurat object
      original_sample <- rv$sample
      
      # Step 1: Ensure mt percentage and QC filtering
      incProgress(0.05, detail = "Running QC filtering...")
      ensure_mt_percent()
      
      rv$sample <- subset(rv$sample, 
                          subset = nFeature_RNA > input$min_features & 
                            nFeature_RNA < input$max_features & 
                            percent.mt < input$max_mt_percent)
      
      # Check if we still have cells
      if(ncol(rv$sample) == 0) {
        showNotification("QC filtering removed all cells. Try less stringent parameters.", type = "error")
        rv$sample <- original_sample
        return()
      }
      
      update_step_status("qc_done", TRUE)
      
      # Step 2: Normalization
      incProgress(0.15, detail = "Running normalization...")
      if(input$normalization_method == "LogNormalize") {
        rv$sample <- NormalizeData(rv$sample, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = input$scale_factor)
      } else if(input$normalization_method == "SCTransform") {
        rv$sample <- SCTransform(rv$sample)
      }
      
      update_step_status("normalization_done", TRUE)
      
      # Step 3: Find variable features (skip if SCTransform was used)
      if(input$normalization_method != "SCTransform") {
        incProgress(0.25, detail = "Finding variable features...")
        rv$sample <- FindVariableFeatures(rv$sample, 
                                        selection.method = "vst", 
                                        nfeatures = input$n_features)
      }
      
      update_step_status("variable_features_done", TRUE)
      
      # Step 4: Scale data
      incProgress(0.35, detail = "Scaling data...")
      if(input$scale_all_genes) {
        all_genes <- rownames(rv$sample)
        rv$sample <- ScaleData(rv$sample, features = all_genes)
      } else {
        rv$sample <- ScaleData(rv$sample)
      }
      
      update_step_status("scaling_done", TRUE)
      
      # Step 5: Run PCA
      incProgress(0.5, detail = "Running PCA...")
      rv$sample <- RunPCA(rv$sample, npcs = input$npcs)
      
      update_step_status("pca_done", TRUE)
      
      # Step 6: Find neighbors
      incProgress(0.65, detail = "Finding neighbors...")
      pc_range <- input$pcs_for_clustering[1]:input$pcs_for_clustering[2]
      rv$sample <- FindNeighbors(rv$sample, dims = pc_range)
      
      update_step_status("neighbors_done", TRUE)
      
      # Step 7: Find clusters
      incProgress(0.8, detail = "Finding clusters...")
      rv$sample <- FindClusters(rv$sample, 
                               resolution = input$resolution,
                               algorithm = as.numeric(input$cluster_algorithm),
                               modularity.fxn = as.numeric(input$modularity_function),
                               n.start = input$cluster_n_start,
                               n.iter = input$cluster_n_iter,
                               group.singletons = input$group_singletons)
      
      update_step_status("clustering_done", TRUE)
      
      # Step 8: Run UMAP
      incProgress(0.9, detail = "Running UMAP...")
      rv$sample <- RunUMAP(rv$sample, 
                          dims = pc_range,
                          min.dist = input$umap_mindist,
                          spread = input$umap_spread)
      
      update_step_status("umap_done", TRUE)
      
      # Update UI elements
      update_cluster_ui()
      
      showNotification(
        paste0("All steps completed successfully. Found ", length(unique(Idents(rv$sample))), " clusters."),
        type = "message"
      )
    }, error = function(e) {
      rv$sample <- original_sample
      showNotification(paste("Error in analysis:", e$message), type = "error")
    })
  })
})