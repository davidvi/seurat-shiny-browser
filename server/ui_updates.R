# UI updates module

# Update UI elements based on reactive values
observe({
  # Update cluster selectors with current cluster names
  updateSelectInput(session, "cluster_name_selector", choices = rv$cluster_names)
  updateSelectInput(session, "cluster_name_selector_analyze", choices = rv$cluster_names)
  updateSelectInput(session, "first_cluster_name_selector", choices = rv$cluster_names)
  updateSelectInput(session, "second_cluster_name_selector", choices = rv$cluster_names)
  
  # Update the delete cluster selector (using pickerInput)
  updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
  
  # Ensure metadata column selector is updated
  if(!is.null(rv$sample)) {
    metadata_cols <- colnames(rv$sample@meta.data)
    rv$metadata_columns <- c("ident", metadata_cols)
    updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
  }
})