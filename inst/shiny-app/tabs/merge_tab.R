merge_tab <- function() {
  tabPanel(
    title = "Integrate Samples",
    h3("Merge and Integrate Samples"),
    
    fluidRow(
      column(6,
        h4("Select folder"),
        selectInput(
          inputId = "merge_folder_selector",
          label = "Choose folder:",
          choices = all_folders
        )
      ),
      column(6,
        h4("Available samples in folder"),
        pickerInput(
          inputId = "merge_files_selector",
          label = "Select samples to merge:",
          choices = NULL,
          multiple = TRUE,
          options = list(
            `actions-box` = TRUE,
            `selected-text-format` = "count > 2"
          )
        )
      )
    ),
    
    hr(),
    
    fluidRow(
      column(12,
        h4("Merge and Integration Options"),
        fluidRow(
          column(4,
            checkboxInput(
              inputId = "add_sample_id_to_cell_names",
              label = "Add sample ID to cell names",
              value = TRUE
            )
          ),
          column(4,
            selectInput(
              inputId = "integration_method",
              label = "Integration Method:",
              choices = c(
                "No Integration" = "none",
                "CCA Integration" = "CCAIntegration",
                "RPCA Integration" = "RPCAIntegration",
                "Harmony Integration" = "HarmonyIntegration",
                "Joint PCA Integration" = "JointPCAIntegration"
                # Note: scVI integration requires additional setup/dependencies,
                # so we'll exclude it for simplicity
              ),
              selected = "RPCAIntegration"
            )
          ),
          column(4,
            numericInput(
              inputId = "integration_dims",
              label = "Number of dimensions:",
              value = 30,
              min = 5,
              max = 50,
              step = 5
            )
          )
        ),
        fluidRow(
          column(4,
            numericInput(
              inputId = "cluster_resolution",
              label = "Clustering resolution:",
              value = 0.8,
              min = 0.1,
              max = 3.0,
              step = 0.1
            )
          ),
          column(4,
            checkboxInput(
              inputId = "run_umap_after_integration",
              label = "Run UMAP after integration",
              value = TRUE
            )
          ),
          column(4,
            checkboxInput(
              inputId = "find_clusters_after_integration",
              label = "Find clusters after integration",
              value = TRUE
            )
          )
        ),
        actionButton(
          inputId = "merge_button",
          label = "Merge and Integrate Samples",
          class = "btn-primary"
        )
      )
    ),
    
    hr(),
    
    fluidRow(
      column(12, 
        h4("Merge Status"),
        verbatimTextOutput("merge_status")
      )
    ),
    
    hr(),
    
    fluidRow(
      column(12,
        conditionalPanel(
          condition = "input.run_umap_after_integration && output.merge_complete",
          h4("Integrated UMAP Visualization"),
          plotOutput("integrated_umap_plot", height = "500px")
        )
      )
    )
  )
}