# Gene search and visualization options module

### search gene ###

observeEvent(input$search_button, {
  rv$gene <- toupper(input$gene)
  # Switch to the single gene visualization tab for each tab that uses the module
  updateTabsetPanel(session, "visualize_viz_tabs", selected = "Single Gene Visualization")
  updateTabsetPanel(session, "edit_viz_tabs", selected = "Single Gene Visualization")
  updateTabsetPanel(session, "analyze_viz_tabs", selected = "Single Gene Visualization")
})

### search multiple genes ###

observeEvent(input$multi_gene_search_button, {
  # Parse the input text to get a list of genes
  gene_text <- input$multi_gene_input
  
  # Split by comma or newline
  genes <- unlist(strsplit(gene_text, "([,\n\r])+"))
  
  # Trim whitespace and convert to uppercase
  genes <- toupper(trimws(genes))
  
  # Remove empty entries
  genes <- genes[genes != ""]
  
  # Handle no genes entered
  if(length(genes) == 0) {
    showNotification("Please enter at least one gene", type = "warning", duration = 10)
    return()
  }
  
  # Save the list of genes to a reactive value
  rv$multiple_genes <- genes
  
  # Provide feedback on how many genes were found
  showNotification(paste(length(genes), "genes selected for visualization"), type = "message", duration = 10)
  
  # Switch to the multi-gene visualization tab in all modules
  updateTabsetPanel(session, "visualize_viz_tabs", selected = "Multi-Gene Visualization")
  updateTabsetPanel(session, "edit_viz_tabs", selected = "Multi-Gene Visualization")
  updateTabsetPanel(session, "analyze_viz_tabs", selected = "Multi-Gene Visualization")
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
  # If metadata column changed, switch to the metadata visualization tab in all modules
  updateTabsetPanel(session, "visualize_viz_tabs", selected = "Metadata Visualization")
  updateTabsetPanel(session, "edit_viz_tabs", selected = "Metadata Visualization")
  updateTabsetPanel(session, "analyze_viz_tabs", selected = "Metadata Visualization")
})


# Initialize multiple_genes in rv if needed
observe({
  if(is.null(rv$multiple_genes)) {
    rv$multiple_genes <- c("CSF1R", "MS4A1", "CD79A")
  }
})