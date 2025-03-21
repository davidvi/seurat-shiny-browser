edit_tab <- function() {
  tabPanel(
    title = "Edit/Rename Clusters",
    fluidRow(
      column(4,
        h3("Rename clusters"),
        selectInput(
          inputId = "cluster_name_selector",
          label = "Select cluster:",
          choices = cluster_names
        ),
        textInput(
          inputId = "new_cluster_name",
          label = "Rename cluster to:",
          value = ""
        ),
        actionButton("rename_button", "Rename cluster"),
        hr(),
        h3("Backup/restore cluster names"),
        fluidRow(
          column(6, actionButton("backup_cluster_button", "Backup cluster names")),
          column(6, actionButton("restore_cluster_button", "Restore cluster names"))
        ),
        p("Backup saves identities to 'idents.backup'. Restore sets identities from 'idents.backup'."),
        hr(),
        h3("Delete clusters"),
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
      ),
      column(8,
        h3("Cluster visualization"),
        fluidRow(
          column(6, plotOutput(outputId = "violin_plot_edit")),
          column(6, plotOutput(outputId = "feature_plot_edit"))
        )
      )
    )
  )
}