# UI

jsCode <- "
$(document).on('keypress', '#gene', function(e) {
  if (e.which == 13) {  // Enter key = keycode 13
    $('#search_button').click();  // trigger search_button click event
  }
});
"

ui <- fluidPage(
  useShinyjs(),
  tags$script(jsCode),
  titlePanel("Single Cell Data Explorer"),
  sidebarLayout(
    sidebarPanel(
      h3("Search for a gene"),
      textInput(
        inputId = "gene",
        label = "Gene Name:",
        value = "CSF1R"
      ),
      actionButton("search_button", "Search gene"),
      h3("Load another sample"),
      selectInput(
        inputId = "sample",
        label = "Select Sample:",
        choices = all_rds_files
      ),
      actionButton("load_button", "Load file"),
      h3("Actions on cell clusters"),
      selectInput(
        inputId = "cluster_name_selector",
        label = "Select cluster:",
        choices = cluster_names
      ),
      actionButton("find_markers", "Find markers for selected cluster"),
      actionButton("find_all_markers", "Find markers for all clusters"),
      textInput(
        inputId = "new_cluster_name",
        label = "Rename cluster to:",
        value = ""
      ),
      actionButton("rename_button", "Rename cluster"),
      selectInput(
        inputId = "second_cluster_name_selector",
        label = "Select second cluster:",
        choices = cluster_names
      ),
      actionButton("compare_button", "Differential analysis two cluster"),
      h3("Save the current sample"),
      textInput(
        inputId = "new_save_name",
        label = "Save file as:",
        value = ""
      ),
      actionButton("save_button", "Save file")
    ),
    mainPanel(
      h4(textOutput("sample_name")),

      # Plots
      fluidRow(
        column(6, plotOutput(outputId = "violin_plot_sample1")),
        column(6, plotOutput(outputId = "feature_plot_sample1"))
      ),
      dataTableOutput("markers_table"),
      downloadButton("download_data", "Download table data")
    )
  )
)
