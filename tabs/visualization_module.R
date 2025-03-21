# Reusable visualization module for Seurat
# This module creates a tabbed visualization interface with single gene, multi-gene, and metadata visualizations

visualization_module <- function(id_prefix) {
  # Create unique IDs for this instance of the module
  violin_id <- paste0(id_prefix, "_violin_plot")
  feature_id <- paste0(id_prefix, "_feature_plot")
  dot_id <- paste0(id_prefix, "_dot_plot")
  dim_id <- paste0(id_prefix, "_dim_plot")
  
  # Return the UI elements
  tabsetPanel(
    id = paste0(id_prefix, "_viz_tabs"),
    tabPanel(
      title = "Single Gene Visualization",
      h3("Gene Expression"),
      fluidRow(
        column(6, plotOutput(outputId = violin_id)),
        column(6, plotOutput(outputId = feature_id))
      )
    ),
    tabPanel(
      title = "Multi-Gene Visualization",
      h3("Dotplot Visualization"),
      p("Shows expression level (color) and percent expressed (size) across clusters for selected genes"),
      plotOutput(outputId = dot_id, height = "500px")
    ),
    tabPanel(
      title = "Metadata Visualization",
      h3("Dimensional Reduction by Metadata"),
      p("For categorical data, cells will be colored by groups. For continuous data, a gradient will be used."),
      plotOutput(outputId = dim_id, height = "500px")
    )
  )
}

# Functions to create the server-side plot rendering for the visualization module
create_violin_plot_output <- function(output, id, rv) {
  output[[id]] <- renderPlot({
    if(is.null(rv$sample) || is.null(rv$gene)) {
      plot.new()
      text(0.5, 0.5, "No sample or gene selected", cex = 1.2)
      return()
    }
    
    tryCatch({
      VlnPlot(rv$sample, features = c(rv$gene), pt.size = 0) + NoLegend() + xlab("")
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error displaying violin plot:", e$message), cex = 1.2)
    })
  })
}

create_feature_plot_output <- function(output, id, rv) {
  output[[id]] <- renderPlot({
    if(is.null(rv$sample) || is.null(rv$gene)) {
      plot.new()
      text(0.5, 0.5, "No sample or gene selected", cex = 1.2)
      return()
    }
    
    tryCatch({
      FeaturePlot(rv$sample, 
                 features = rv$gene, 
                 reduction = rv$selected_reduction, 
                 label = TRUE) + 
        NoLegend()
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Cannot display feature on", rv$selected_reduction, "\nError:", e$message), cex = 1.2)
    })
  })
}

create_dot_plot_output <- function(output, id, rv, input) {
  output[[id]] <- renderPlot({
    if(is.null(rv$sample) || is.null(rv$multiple_genes) || length(rv$multiple_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "No sample or genes selected", cex = 1.2)
      return()
    }
    
    # Validate genes exist in the dataset
    valid_genes <- rv$multiple_genes[rv$multiple_genes %in% rownames(rv$sample)]
    if(length(valid_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "None of the selected genes were found in the dataset", cex = 1.2)
      return()
    }
    
    # If some genes were invalid, show a notification
    if(length(valid_genes) < length(rv$multiple_genes)) {
      missing_genes <- setdiff(rv$multiple_genes, valid_genes)
      showNotification(
        paste("Some genes were not found in the dataset:", 
              paste(missing_genes, collapse = ", ")), 
        type = "warning"
      )
    }
    
    tryCatch({
      # Get the current cluster identities (use default ordering)
      cluster_ids <- levels(Idents(rv$sample))
      
      # Create the dotplot with ordered clusters
      DotPlot(rv$sample, 
             features = valid_genes, 
             cols = c("lightgrey", "blue"),
             dot.scale = 6,
             cluster.idents = FALSE) + 
        RotatedAxis() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right"
        ) +
        labs(title = "Gene Expression Dotplot")
        
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error creating dotplot:", e$message), cex = 1.2)
    })
  })
}

create_dim_plot_output <- function(output, id, rv, input) {
  output[[id]] <- renderPlot({
    if(is.null(rv$sample)) {
      plot.new()
      text(0.5, 0.5, "No sample loaded", cex = 1.2)
      return()
    }
    
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
}

# Function to setup all visualization outputs for a given tab
setup_visualization_outputs <- function(output, id_prefix, rv, input) {
  create_violin_plot_output(output, paste0(id_prefix, "_violin_plot"), rv)
  create_feature_plot_output(output, paste0(id_prefix, "_feature_plot"), rv)
  create_dot_plot_output(output, paste0(id_prefix, "_dot_plot"), rv, input)
  create_dim_plot_output(output, paste0(id_prefix, "_dim_plot"), rv, input)
}