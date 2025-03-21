load_save_tab <- function() {
  tabPanel(
    title = "Load/Save",
    fluidRow(
      column(4,
        h3("Select Folder"),
        selectInput(
          inputId = "folder_selector",
          label = "Target Folder:",
          choices = all_folders,
          selected = "."
        ),
        hr(),
        h3("Save Current Sample"),
        textInput(
          inputId = "new_save_name",
          label = "Save file as:",
          value = ""
        ),
        p("Sample will be saved to the selected folder"),
        actionButton("save_button", "Save file")
      ),
      column(8,
        h3("Load Sample"),
        p("Files from all folders are shown with their folder prefix"),
        selectInput(
          inputId = "sample",
          label = "Select Sample:",
          choices = all_rds_files,
          selected = if(length(all_rds_files) > 0) all_rds_files[1] else NULL
        ),
        actionButton("load_button", "Load selected file")
      )
    )
  )
}