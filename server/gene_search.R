# Gene search and visualization options module

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