source("tabs/visualization_module.R")

analyze_tab <- function() {
  tabPanel(
    title = "Differential Expression",
    fluidRow(
      column(4,
        wellPanel(
          h3("Find Marker Genes"),
          selectInput(
            inputId = "cluster_name_selector_analyze",
            label = "Select cluster:",
            choices = cluster_names
          ),
          div(
            style = "display: flex; justify-content: space-between; margin-top: 15px;",
            actionButton("find_markers", "Find markers for selected cluster", class = "btn-primary"),
            actionButton("find_all_markers", "Find markers for all clusters", class = "btn-info")
          )
        ),
        wellPanel(
          h3("Differential Expression Analysis"),
          selectInput(
            inputId = "first_cluster_name_selector",
            label = "First cluster:",
            choices = cluster_names
          ),
          selectInput(
            inputId = "second_cluster_name_selector",
            label = "Second cluster:",
            choices = cluster_names
          ),
          actionButton("compare_button", "Run differential analysis", class = "btn-success")
        )
      ),
      column(8,
        tabsetPanel(
          id = "analyze_results_tabs",
          tabPanel(
            title = "Visualization",
            # Use the modular visualization component with prefix "analyze"
            visualization_module("analyze")
          ),
          tabPanel(
            title = "Differential Expression Results",
            h3("Marker Genes Table"),
            p("Click on any gene in the table to visualize its expression"),
            dataTableOutput("markers_table"),
            downloadButton("download_data", "Download table data", class = "btn-primary")
          )
        )
      )
    )
  )
}