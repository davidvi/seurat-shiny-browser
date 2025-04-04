source("tabs/visualization_module.R")

edit_tab <- function() {
  tabPanel(
    title = "Edit/Rename Clusters",
    fluidRow(
      column(4,
        wellPanel(
          h3("Rename clusters"),
          pickerInput(
            inputId = "cluster_name_selector",
            label = "Select clusters to rename:",
            choices = cluster_names,
            multiple = TRUE,
            options = list(
              `actions-box` = TRUE,
              `selected-text-format` = "count > 2"
            )
          ),
          textInput(
            inputId = "new_cluster_name",
            label = "Rename selected cluster(s) to:",
            value = ""
          ),
          div(
            style = "margin: 15px 0;",
            checkboxInput(
              inputId = "add_number_suffix",
              label = "Add numeric suffix when renaming multiple clusters",
              value = TRUE
            )
          ),
          actionButton("rename_button", "Rename cluster(s)", class = "btn-primary")
        ),
        wellPanel(
          h3("Backup/restore cluster names"),
          fluidRow(
            column(6, actionButton("backup_cluster_button", "Backup cluster names", class = "btn-info")),
            column(6, actionButton("restore_cluster_button", "Restore cluster names", class = "btn-warning"))
          ),
          p("Backup saves identities to 'idents.backup'. Restore sets identities from 'idents.backup'.")
        ),
        wellPanel(
          h3("Delete clusters"),
          div(
            style = "background-color: #ffeeee; padding: 15px; border-radius: 4px; border-left: 4px solid #dc3545;",
            pickerInput(
              inputId = "delete_cluster_selector",
              label = "Select clusters to delete:",
              choices = cluster_names,
              multiple = TRUE,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 2"
              )
            ),
            div(
              style = "margin-bottom: 15px;",
              checkboxInput(
                inputId = "confirm_delete_cluster",
                label = "I understand this will permanently remove these cells",
                value = FALSE
              )
            ),
            actionButton("delete_cluster_button", "Delete selected clusters", 
                        class = "btn-danger")
          )
        ),
        wellPanel(
          h3("Edit orig.ident"),
          selectInput(
            inputId = "orig_ident_selector",
            label = "Select orig.ident value to rename:",
            choices = NULL
          ),
          textInput(
            inputId = "new_orig_ident",
            label = "Rename to:",
            value = ""
          ),
          actionButton("rename_orig_ident_button", "Rename orig.ident", class = "btn-primary")
        ),
        wellPanel(
          h3("Add metadata"),
          div(
            style = "margin-bottom: 15px;",
            p("Add log2-transformed UMI counts as a new metadata column."),
            actionButton("add_log2umi_button", "Calculate log2UMI_cell", class = "btn-info")
          )
        )
      ),
      column(8,
        # Use the modular visualization component with prefix "edit"
        visualization_module("edit")
      )
    )
  )
}