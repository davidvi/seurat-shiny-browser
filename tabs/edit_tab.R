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
        h3("Other cluster operations"),
        fluidRow(
          column(6, actionButton("backup_cluster_button", "Backup cluster names")),
          column(6, actionButton("restore_cluster_button", "Restore cluster names"))
        ),
        p("Backup saves identities to 'idents.backup'. Restore sets identities from 'idents.backup'.")
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