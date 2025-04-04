# Visualization module
source("tabs/visualization_module.R")

# Setup visualization outputs for each tab
# This replaces all the individual plot renderers with the modular approach

# Visualize tab plots
setup_visualization_outputs(output, "visualize", rv, input)

# Edit tab plots
setup_visualization_outputs(output, "edit", rv, input)

# Differential Expression tab plots
setup_visualization_outputs(output, "analyze", rv, input)

# Update the selected tab when finding markers to show the results
observeEvent(input$find_markers, {
  updateTabsetPanel(session, "analyze_results_tabs", selected = "Differential Expression Results")
})

observeEvent(input$find_all_markers, {
  updateTabsetPanel(session, "analyze_results_tabs", selected = "Differential Expression Results")
})

observeEvent(input$compare_button, {
  updateTabsetPanel(session, "analyze_results_tabs", selected = "Differential Expression Results")
})

# Select visualization tab when a gene in the markers table is clicked
observeEvent(input$markers_table_cell_clicked, {
  tryCatch({
    info <- input$markers_table_cell_clicked
    
    # Validate the info object is not NULL and contains necessary fields
    if (is.null(info) || is.null(info$col) || is.null(info$row)) {
      return()
    }
    
    # Check if we have valid markers data
    if (is.null(rv$markers) || nrow(rv$markers) == 0) {
      return()
    }
    
    # Check column index validity
    if (info$col == 0 && info$row >= 0 && info$row < nrow(rv$markers)) {
      # Check if the gene/feature column exists
      gene_col <- NULL
      possible_gene_cols <- c("gene", "genes", "feature", "features", "Gene", "Features")
      
      for (col in possible_gene_cols) {
        if (col %in% colnames(rv$markers)) {
          gene_col <- col
          break
        }
      }
      
      # If no gene column found, warn and exit
      if (is.null(gene_col)) {
        message("Warning: No gene/feature column found in markers data.frame")
        return()
      }
      
      # Get the gene from the table - row indices are 0-based
      gene_selected <- rv$markers[info$row, gene_col]
      
      if (!is.null(gene_selected) && !is.na(gene_selected) && gene_selected != "") {
        # Update the selected gene
        rv$gene <- gene_selected
        updateTextInput(session, "gene", value = gene_selected)
        
        # Switch to the visualization tab and set view to single gene
        updateTabsetPanel(session, "analyze_results_tabs", selected = "Visualization")
        updateTabsetPanel(session, "analyze_viz_tabs", selected = "Single Gene Visualization")
        
        # Clear the row selection by proxy
        proxy <- dataTableProxy("markers_table")
        proxy %>% selectRows(NULL)
      }
    }
  }, error = function(e) {
    message("Error in markers_table_cell_clicked handler: ", e$message)
  })
})