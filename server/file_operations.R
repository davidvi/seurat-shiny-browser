# File operations module

### select folder ###

observeEvent(input$folder_selector, {
  # Update the current folder
  current_folder <- input$folder_selector
  rv$current_folder <- current_folder
  
  # Get files from the selected folder
  files <- get_rds_files(current_folder)
  rv$all_rds_files <- files
  
  # Update the sample selectors
  updateSelectInput(session, "sample", choices = files)
  updateSelectInput(session, "delete_sample_selector", choices = files)
  
  message("Updated file list for folder: ", current_folder, ", found ", length(files), " files")
})

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
      
      # Make sure delete cluster selector is updated
      updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
    } else {
      rv$cluster_names <- c()
    }
    
    rv$markers <- NULL
    
    # Update available reductions after loading new sample
    tryCatch({
      reductions_list <- Reductions(rv$sample)
      if(length(reductions_list) > 0) {
        updateSelectInput(session, "reduction_selector", choices = reductions_list)
        
        # Prefer umap if it's available, otherwise use the first reduction
        if("umap" %in% reductions_list) {
          rv$selected_reduction <- "umap"
          updateSelectInput(session, "reduction_selector", choices = reductions_list, selected = "umap")
        } else {
          rv$selected_reduction <- reductions_list[1]
        }
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
  
  # Update the sample selectors
  updateSelectInput(session, "sample", choices = files)
  updateSelectInput(session, "delete_sample_selector", choices = files)
  
  message("Refreshed file list for folder: ", rv$current_folder, ", now contains ", length(files), " files")
}

# Function to refresh the folders list
refresh_folders <- function() {
  # Get folders in data directory
  all_folders <- list.dirs(path = base_folder, full.names = FALSE, recursive = FALSE)
  
  # Always include the root folder
  all_folders <- c(".", all_folders)
  
  # Remove empty values
  all_folders <- all_folders[all_folders != ""]
  
  # Update the reactive value
  rv$all_folders <- all_folders
  
  # Update the folder selector
  updateSelectInput(session, "folder_selector", choices = all_folders, selected = rv$current_folder)
  updateSelectInput(session, "merge_folder_selector", choices = all_folders)
  
  message("Refreshed folder list, found ", length(all_folders), " folders")
}

# Create new folder button handler
observeEvent(input$create_folder_button, {
  # Get the folder name
  new_folder_name <- input$new_folder_name
  
  # Validate folder name
  if (is.null(new_folder_name) || new_folder_name == "") {
    showNotification("Folder name cannot be empty", type = "error")
    return()
  }
  
  # Check if folder name contains invalid characters
  if (grepl("[\\/:*?\"<>|]", new_folder_name)) {
    showNotification("Folder name contains invalid characters", type = "error")
    return()
  }
  
  # Create the folder
  folder_path <- file.path(base_folder, new_folder_name)
  
  # Check if folder already exists
  if (dir.exists(folder_path)) {
    showNotification(paste("Folder", new_folder_name, "already exists"), type = "warning")
    return()
  }
  
  # Create the folder
  tryCatch({
    dir.create(folder_path, recursive = TRUE)
    
    # Refresh the folders list
    refresh_folders()
    
    # Select the new folder
    updateSelectInput(session, "folder_selector", selected = new_folder_name)
    rv$current_folder <- new_folder_name
    
    # Refresh the file list
    refresh_folder_files()
    
    # Clear the folder name input
    updateTextInput(session, "new_folder_name", value = "")
    
    showNotification(paste("Folder", new_folder_name, "created successfully"), type = "message")
  }, error = function(e) {
    showNotification(paste("Error creating folder:", e$message), type = "error")
  })
})

# Delete sample button handler
observeEvent(input$delete_sample_button, {
  # Check if deletion is confirmed
  if (!input$confirm_delete_sample) {
    showNotification("Please confirm deletion by checking the confirmation box", type = "warning")
    return()
  }
  
  # Get the file to delete
  file_to_delete <- input$delete_sample_selector
  
  # Validate file name
  if (is.null(file_to_delete) || file_to_delete == "") {
    showNotification("No file selected for deletion", type = "error")
    return()
  }
  
  # Construct the file path
  file_path <- file.path(base_folder, rv$current_folder, file_to_delete)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    showNotification(paste("File", file_to_delete, "does not exist"), type = "error")
    return()
  }
  
  # Delete the file
  tryCatch({
    # Check if the file to delete is the currently loaded sample
    if (!is.null(rv$sample_name) && (paste0(rv$current_folder, "/", file_to_delete) == rv$sample_name)) {
      showNotification("Cannot delete the currently loaded sample", type = "error")
      return()
    }
    
    file.remove(file_path)
    
    # Refresh the file list
    refresh_folder_files()
    
    # Reset the confirmation checkbox
    updateCheckboxInput(session, "confirm_delete_sample", value = FALSE)
    
    showNotification(paste("File", file_to_delete, "deleted successfully"), type = "message")
  }, error = function(e) {
    showNotification(paste("Error deleting file:", e$message), type = "error")
  })
})

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