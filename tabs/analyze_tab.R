analyze_tab <- function() {
  tabPanel(
    title = "Differential Expression",
    fluidRow(
      column(4,
        h3("Find Marker Genes"),
        selectInput(
          inputId = "cluster_name_selector_analyze",
          label = "Select cluster:",
          choices = cluster_names
        ),
        actionButton("find_markers", "Find markers for selected cluster"),
        actionButton("find_all_markers", "Find markers for all clusters"),
        hr(),
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
        actionButton("compare_button", "Run differential analysis")
      ),
      column(8,
        h3("Gene Expression Visualization"),
        fluidRow(
          column(6, plotOutput(outputId = "violin_plot_analyze")),
          column(6, plotOutput(outputId = "feature_plot_analyze"))
        ),
        hr(),
        h3("Results"),
        dataTableOutput("markers_table"),
        downloadButton("download_data", "Download table data")
      )
    )
  )
}