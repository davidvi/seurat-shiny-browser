# File operations module

# Add an observer to refresh the folders list when the session starts
observe({
  # Get folders in data directory and refresh the list
  refresh_folders()
  
  # Log message
  message("Refreshed folder list on app initialization")
})

### Helper functions for raw data loading ###

# Function to check if a directory is a valid 10X Genomics data directory
is_valid_10x_dir <- function(dir_path) {
  # Check for matrix.mtx file (may be .gz compressed)
  has_matrix <- any(grepl("matrix\\.mtx(\\.gz)?$", list.files(dir_path, recursive = TRUE)))
  
  # Check for barcodes.tsv file (may be .gz compressed)
  has_barcodes <- any(grepl("barcodes\\.tsv(\\.gz)?$", list.files(dir_path, recursive = TRUE)))
  
  # Check for genes.tsv or features.tsv (may be .gz compressed)
  has_genes <- any(grepl("(genes|features)\\.tsv(\\.gz)?$", list.files(dir_path, recursive = TRUE)))
  
  # Return TRUE if all required files are present
  return(has_matrix && has_barcodes && has_genes)
}

# Function to find all subdirectories in a given folder that contain 10X data
find_10x_dirs <- function(parent_dir) {
  if(!dir.exists(parent_dir)) {
    return(character(0))
  }
  
  # Get all subdirectories (just one level deep for performance)
  all_subdirs <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  
  # Filter for valid 10X directories
  valid_10x_dirs <- all_subdirs[sapply(all_subdirs, is_valid_10x_dir)]
  
  # Return relative paths
  return(basename(valid_10x_dirs))
}

# Function to calculate percentage of mitochondrial genes
calculate_mt_percent <- function(seurat_obj) {
  tryCatch({
    # Get genes starting with MT- (mitochondrial genes)
    mt_genes <- grep("^MT-", rownames(seurat_obj@assays$RNA), value = TRUE)
    
    if(length(mt_genes) > 0) {
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      message("Calculated mitochondrial percentage for ", length(mt_genes), " genes")
    } else {
      message("No mitochondrial genes found with pattern '^MT-'")
    }
    
    return(seurat_obj)
  }, error = function(e) {
    message("Error calculating mitochondrial percentage: ", e$message)
    return(seurat_obj)
  })
}

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
  updatePickerInput(session, "move_files_selector", choices = files, selected = NULL)
  
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
    
    # Update the multiple genes text area with default genes
    updateTextAreaInput(session, "multi_gene_input", 
                       value = paste(rv$multiple_genes, collapse = "\n"))
    
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

### Raw data 10X loading ###

# Dynamic UI for raw data subdirectories
output$raw_data_subdirs <- renderUI({
  req(input$raw_data_parent_folder)
  
  # Get the full path of the parent folder
  parent_path <- file.path(base_folder, input$raw_data_parent_folder)
  
  # Find subdirectories that contain 10X data
  valid_dirs <- find_10x_dirs(parent_path)
  
  if(length(valid_dirs) == 0) {
    div(
      style = "color: #856404; background-color: #fff3cd; padding: 10px; border-radius: 4px; margin-top: 10px;",
      icon("exclamation-triangle"), 
      "No valid 10X Genomics data directories found. Each directory should contain matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv files."
    )
  } else {
    selectInput(
      inputId = "raw_data_subdir",
      label = "Select 10X Data Directory:",
      choices = valid_dirs,
      selected = if(length(valid_dirs) > 0) valid_dirs[1] else NULL
    )
  }
})

# Observer for parent folder change to scan for 10X data dirs
observeEvent(input$raw_data_parent_folder, {
  # When parent folder changes, ensure the UI updates
  invalidateLater(100)
})

# Load Raw 10X Data button handler
observeEvent(input$load_raw_button, {
  # Validate required inputs
  req(input$raw_data_parent_folder, input$raw_data_subdir, input$project_name)
  
  # Full path to the 10X data directory
  data_dir <- file.path(base_folder, input$raw_data_parent_folder, input$raw_data_subdir)
  
  # Validate directory exists and contains valid data
  if(!dir.exists(data_dir) || !is_valid_10x_dir(data_dir)) {
    showNotification(
      "Invalid 10X Genomics data directory. Make sure it contains matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv files.",
      type = "error"
    )
    return()
  }
  
  # Get parameters
  project_name <- input$project_name
  min_cells <- input$min_cells
  min_features <- input$min_features
  calc_mt <- input$calc_mt_percent
  
  withProgress(message = "Loading 10X data...", value = 0, {
    tryCatch({
      # Step 1: Read in the 10X data
      incProgress(0.1, detail = "Reading 10X data matrix")
      counts_data <- Read10X(data.dir = data_dir)
      
      # Step 2: Create Seurat object
      incProgress(0.3, detail = "Creating Seurat object")
      seurat_obj <- CreateSeuratObject(
        counts = counts_data, 
        project = project_name, 
        min.cells = min_cells, 
        min.features = min_features
      )
      
      # Step 3: Add mitochondrial percentage if requested
      if(calc_mt) {
        incProgress(0.6, detail = "Calculating mitochondrial percentage")
        seurat_obj <- calculate_mt_percent(seurat_obj)
      }
      
      # Step 4: Update the sample in reactive values
      incProgress(0.8, detail = "Finalizing")
      rv$sample <- seurat_obj
      rv$sample_name <- paste0(input$raw_data_parent_folder, "/", input$raw_data_subdir, " (unsaved)")
      
      # Set cluster names to empty since the data hasn't been clustered yet
      rv$cluster_names <- c()
      updatePickerInput(session, "delete_cluster_selector", choices = rv$cluster_names)
      
      # Clear previously computed markers
      rv$markers <- NULL
      
      # Step 5: Update available reductions (likely none for raw data)
      reductions_list <- c()
      tryCatch({
        reductions_list <- Reductions(seurat_obj)
      }, error = function(e) {
        # No reductions available in raw data
      })
      
      if(length(reductions_list) > 0) {
        updateSelectInput(session, "reduction_selector", choices = reductions_list)
        rv$selected_reduction <- reductions_list[1]
      } else {
        # No dimensional reductions yet for raw data
        updateSelectInput(session, "reduction_selector", choices = c("No reductions available"))
        rv$selected_reduction <- NULL
      }
      
      # Update available metadata columns
      if(!is.null(seurat_obj)) {
        metadata_cols <- colnames(seurat_obj@meta.data)
        rv$metadata_columns <- c("ident", metadata_cols)
        updateSelectInput(session, "metadata_column_selector", choices = rv$metadata_columns)
      }
      
      incProgress(1.0, detail = "Loading complete")
      
      # Provide success notification
      cells_loaded <- ncol(seurat_obj)
      features_loaded <- nrow(seurat_obj)
      showNotification(
        paste("Successfully loaded 10X data with", cells_loaded, "cells and", features_loaded, "features."),
        type = "message"
      )
      
      # Switch to the Normalize/Cluster tab, since that's the next logical step for raw data
      updateTabsetPanel(session, "main_tabs", selected = "Normalize/Cluster")
      
    }, error = function(e) {
      showNotification(
        paste("Error loading 10X data:", e$message),
        type = "error",
        duration = 10
      )
    })
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
  updatePickerInput(session, "move_files_selector", choices = files)
  
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
  
  # Update the folder selectors
  updateSelectInput(session, "folder_selector", choices = all_folders, selected = rv$current_folder)
  updateSelectInput(session, "merge_folder_selector", choices = all_folders)
  updateSelectInput(session, "move_target_folder", choices = all_folders)
  
  message("Refreshed folder list, found ", length(all_folders), " folders")
}

# Create new folder button handler
observeEvent(input$create_folder_button, {
  # Get the folder name
  new_folder_name <- input$new_folder_name
  
  # Validate folder name
  if (is.null(new_folder_name) || new_folder_name == "") {
    showNotification("Folder name cannot be empty", type = "error", duration = 10)
    return()
  }
  
  # Check if folder name contains invalid characters
  if (grepl("[\\/:*?\"<>|]", new_folder_name)) {
    showNotification("Folder name contains invalid characters", type = "error", duration = 10)
    return()
  }
  
  # Create the folder
  folder_path <- file.path(base_folder, new_folder_name)
  
  # Check if folder already exists
  if (dir.exists(folder_path)) {
    showNotification(paste("Folder", new_folder_name, "already exists"), type = "warning", duration = 10)
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
    
    showNotification(paste("Folder", new_folder_name, "created successfully"), type = "message", duration = 10)
  }, error = function(e) {
    showNotification(paste("Error creating folder:", e$message), type = "error", duration = 10)
  })
})

# Move files button handler
observeEvent(input$move_files_button, {
  # Get the selected files
  selected_files <- input$move_files_selector
  
  # Get the target folder
  target_folder <- input$move_target_folder
  
  # Validate inputs
  if (is.null(selected_files) || length(selected_files) == 0) {
    showNotification("No files selected for moving", type = "warning", duration = 10)
    return()
  }
  
  if (is.null(target_folder) || target_folder == "") {
    showNotification("No target folder selected", type = "error", duration = 10)
    return()
  }
  
  # Check if target folder is the same as current folder
  if (target_folder == rv$current_folder) {
    showNotification("Cannot move files to the same folder", type = "warning", duration = 10)
    return()
  }
  
  # Check if target folder exists
  target_path <- file.path(base_folder, target_folder)
  if (!dir.exists(target_path)) {
    showNotification(paste("Target folder", target_folder, "does not exist"), type = "error", duration = 10)
    return()
  }
  
  # Create a progress bar
  withProgress(message = "Moving files", value = 0, {
    # Track successful and failed moves
    success_count <- 0
    failed_moves <- c()
    
    # Set the total number of steps based on the number of files
    total_steps <- length(selected_files)
    
    # Process each file
    for (i in seq_along(selected_files)) {
      file_name <- selected_files[i]
      
      # Update progress
      incProgress(1/total_steps, detail = paste("Moving", file_name))
      
      # Construct the source and destination paths
      source_path <- file.path(base_folder, rv$current_folder, file_name)
      dest_path <- file.path(base_folder, target_folder, file_name)
      
      # Check if source file exists
      if (!file.exists(source_path)) {
        failed_moves <- c(failed_moves, file_name)
        next
      }
      
      # Check if file already exists in target folder
      if (file.exists(dest_path)) {
        # Create a new name with timestamp to avoid conflicts
        new_file_name <- paste0(
          tools::file_path_sans_ext(file_name),
          "_", format(Sys.time(), "%Y%m%d%H%M%S"),
          ".", tools::file_ext(file_name)
        )
        dest_path <- file.path(base_folder, target_folder, new_file_name)
      }
      
      # Move the file
      tryCatch({
        # Check if the file to move is the currently loaded sample
        if (!is.null(rv$sample_name) && (paste0(rv$current_folder, "/", file_name) == rv$sample_name)) {
          failed_moves <- c(failed_moves, file_name)
          showNotification("Cannot move the currently loaded sample", type = "error", duration = 10)
          next
        }
        
        # Copy file to destination
        file.copy(source_path, dest_path)
        
        # Delete the source file if copy was successful
        if (file.exists(dest_path)) {
          file.remove(source_path)
          success_count <- success_count + 1
        } else {
          failed_moves <- c(failed_moves, file_name)
        }
      }, error = function(e) {
        failed_moves <- c(failed_moves, file_name)
        message("Error moving file ", file_name, ": ", e$message)
      })
    }
    
    # Refresh the file lists for both folders
    refresh_folder_files()
    
    # Display completion notification
    if (success_count > 0) {
      showNotification(paste(success_count, "files moved successfully"), type = "message", duration = 10)
    }
    
    if (length(failed_moves) > 0) {
      showNotification(paste("Failed to move", length(failed_moves), "files"), type = "warning", duration = 10)
    }
  })
})

# Delete sample button handler
observeEvent(input$delete_sample_button, {
  # Check if deletion is confirmed
  if (!input$confirm_delete_sample) {
    showNotification("Please confirm deletion by checking the confirmation box", type = "warning", duration = 10)
    return()
  }
  
  # Get the file to delete
  file_to_delete <- input$delete_sample_selector
  
  # Validate file name
  if (is.null(file_to_delete) || file_to_delete == "") {
    showNotification("No file selected for deletion", type = "error", duration = 10)
    return()
  }
  
  # Construct the file path
  file_path <- file.path(base_folder, rv$current_folder, file_to_delete)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    showNotification(paste("File", file_to_delete, "does not exist"), type = "error", duration = 10)
    return()
  }
  
  # Delete the file
  tryCatch({
    # Check if the file to delete is the currently loaded sample
    if (!is.null(rv$sample_name) && (paste0(rv$current_folder, "/", file_to_delete) == rv$sample_name)) {
      showNotification("Cannot delete the currently loaded sample", type = "error", duration = 10)
      return()
    }
    
    file.remove(file_path)
    
    # Refresh the file list
    refresh_folder_files()
    
    # Reset the confirmation checkbox
    updateCheckboxInput(session, "confirm_delete_sample", value = FALSE)
    
    showNotification(paste("File", file_to_delete, "deleted successfully"), type = "message", duration = 10)
  }, error = function(e) {
    showNotification(paste("Error deleting file:", e$message), type = "error", duration = 10)
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
    write_rds(rv$sample, file = save_path)
    
    incProgress(0.8, "Setting new name")
    
    # Update the sample name with folder prefix
    rv$sample_name <- paste0(folder, "/", file_name)
    
    # Refresh the file list
    refresh_folder_files()
    incProgress(0.9, "Done")
    
    # Display success notification
    showNotification(paste("File saved as", file_name), type = "message", duration = 10)
  })
})

### Quick Save ###

# Observer for the quick save button
observeEvent(input$quick_save_button, {
  req(rv$sample, rv$sample_name) # Ensure an object (rv$sample) and sample name exist

  # Get the relative path from rv$sample_name
  relative_path <- rv$sample_name

  # Construct the full path using base_folder
  full_path <- file.path(base_folder, relative_path)

  # Check if the path is valid and contains a directory and filename
  if (is.null(full_path) || full_path == "" || is.null(dirname(full_path)) || is.null(basename(full_path))) {
    showNotification("Invalid sample name or path for saving.", type = "error")
    return()
  }

  # Extract directory and filename for messages
  folder_path <- dirname(full_path)
  file_name <- basename(full_path)

  # Make sure the file name ends with .rds (should already be the case if loaded correctly)
  if (!endsWith(file_name, ".rds")) {
     showNotification(paste("Invalid file extension for:", file_name, "- expected .rds"), type = "warning")
     # Optionally, handle this case differently
     return() # Stop if extension is wrong
  }

  withProgress(message = "Quick Saving Sample", detail = "Please wait...", value = 0, {
    incProgress(0.1, "Saving Sample")

    # Create the directory if it doesn't exist (safety check)
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
      showNotification(paste("Created directory:", folder_path), type = "message")
    }

    # Save the Seurat object (using rv$sample)
    tryCatch({
      saveRDS(rv$sample, file = full_path) # Use rv$sample and the full path
      incProgress(0.9, "Finalizing")
      showNotification(
        paste("Sample saved successfully as", file_name, "in", basename(dirname(relative_path))), # Show relative folder
        type = "message"
      )
    }, error = function(e) {
      showNotification(
        paste("Error saving sample:", e$message),
        type = "error",
        duration = 10
      )
    })
  })

  # No need to refresh file list as we are overwriting
})