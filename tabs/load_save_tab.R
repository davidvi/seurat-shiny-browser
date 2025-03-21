load_save_tab <- function() {
  tabPanel(
    title = "Load/Save",
    fluidRow(
      column(4,
        wellPanel(
          h3("Folder Management"),
          
          # Select folder
          selectInput(
            inputId = "folder_selector",
            label = "Target Folder:",
            choices = all_folders,
            selected = "."
          ),
          
          # Create new folder
          h4("Create New Folder"),
          div(
            style = "display: flex; align-items: flex-end;",
            div(
              style = "flex-grow: 1; margin-right: 10px;",
              textInput(
                inputId = "new_folder_name",
                label = "New Folder Name:",
                value = ""
              )
            ),
            div(
              actionButton(
                inputId = "create_folder_button",
                label = "Create Folder",
                class = "btn-info"
              )
            )
          ),
          
          hr(),
          
          h3("Save Current Sample"),
          textInput(
            inputId = "new_save_name",
            label = "Save file as:",
            value = "",
            placeholder = "Leave empty to use current name"
          ),
          p("Sample will be saved to the selected folder"),
          actionButton("save_button", "Save file", class = "btn-success")
        )
      ),
      column(8,
        wellPanel(
          h3("Load Sample"),
          p("Select a file from the current folder to load"),
          selectInput(
            inputId = "sample",
            label = "Select Sample:",
            choices = all_rds_files,
            selected = if(length(all_rds_files) > 0) all_rds_files[1] else NULL
          ),
          actionButton("load_button", "Load selected file", class = "btn-primary"),
          
          hr(),
          
          h3("Move Files"),
          div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 4px; border-left: 4px solid #17a2b8;",
            p("Select files to move to another folder"),
            
            # Multi-select files to move
            pickerInput(
              inputId = "move_files_selector",
              label = "Select Files to Move:",
              choices = all_rds_files,
              multiple = TRUE,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 2"
              )
            ),
            
            # Target folder selector
            selectInput(
              inputId = "move_target_folder",
              label = "Target Folder:",
              choices = all_folders
            ),
            
            actionButton(
              inputId = "move_files_button",
              label = "Move Selected Files",
              class = "btn-info"
            )
          ),
          
          hr(),
          
          h3("Delete Sample"),
          div(
            style = "background-color: #ffeeee; padding: 15px; border-radius: 4px; border-left: 4px solid #dc3545;",
            selectInput(
              inputId = "delete_sample_selector",
              label = "Select Sample to Delete:",
              choices = all_rds_files
            ),
            div(
              style = "margin-bottom: 15px;",
              checkboxInput(
                inputId = "confirm_delete_sample",
                label = "I understand this will permanently delete the selected file",
                value = FALSE
              )
            ),
            actionButton(
              inputId = "delete_sample_button",
              label = "Delete Sample",
              class = "btn-danger"
            )
          )
        )
      )
    )
  )
}