visualize_tab <- function() {
  tabPanel(
    title = "Visualize",
    fluidRow(
      column(4,
        h3("Search for a gene"),
        textInput(
          inputId = "gene",
          label = "Gene Name:",
          value = "CSF1R"
        ),
        actionButton("search_button", "Search gene"),
        hr(),
        h3("Dimensional Reduction"),
        selectInput(
          inputId = "reduction_selector",
          label = "Select reduction:",
          choices = available_reductions,
          selected = available_reductions[1]
        ),
        hr(),
        h3("Metadata Visualization"),
        selectInput(
          inputId = "metadata_column_selector",
          label = "Color cells by:",
          choices = c("ident")
        )
      ),
      column(8,
        h3("Gene Expression Visualization"),
        fluidRow(
          column(6, plotOutput(outputId = "violin_plot_visualize")),
          column(6, plotOutput(outputId = "feature_plot_visualize"))
        ),
        hr(),
        h3("Dimensional Reduction by Metadata"),
        p("For categorical data, cells will be colored by groups. For continuous data, a gradient will be used."),
        plotOutput(outputId = "dim_plot_visualize", height = "300px")
      )
    )
  )
}