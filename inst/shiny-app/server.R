# Main server file that sources all modules

server <- function(input, output, session) {
  # Initialize reactive values
  rv <- reactiveValues(
    sample = sample,
    sample_name = sample_name,
    cluster_names = cluster_names,
    gene = "CSF1R",
    multiple_genes = c("CSF1R", "MS4A1", "CD79A"),
    all_rds_files = all_rds_files,
    current_folder = ".",
    all_folders = all_folders,
    markers = data.frame(),
    selected_reduction = if(length(available_reductions) > 0) available_reductions[1] else "umap",
    metadata_columns = if(!is.null(sample)) colnames(sample@meta.data) else c("ident"),
    merge_complete = FALSE,
    merge_log = "",
    integrated_reduction = "",
    force_refresh = 0, # Used to force visualization updates
    singleR_results = NULL # For storing SingleR results
  )
  
  # Source all server modules
  source("server/cluster_operations.R", local = TRUE)
  source("server/gene_search.R", local = TRUE)
  source("server/file_operations.R", local = TRUE)
  source("server/differential_expression.R", local = TRUE)
  source("server/merge_integration.R", local = TRUE)
  source("server/visualizations.R", local = TRUE)
  source("server/ui_updates.R", local = TRUE)
  source("server/clustering.R", local = TRUE)
  source("server/auto_naming.R", local = TRUE)
  
  # Sample name output - shared across tabs
  output$sample_name <- renderText({
    paste("Loaded Sample: ", rv$sample_name)
  })
  
  # Auto Name Tab functionality
  
  # Reactive values to store SingleR results
  singleR_values <- reactiveValues(
    results = NULL,
    running = FALSE,
    completed = FALSE,
    error = NULL,
    cluster_mapping = NULL
  )
  
  # Output for showing the running status
  output$singleR_running <- reactive({
    singleR_values$running
  })
  outputOptions(output, "singleR_running", suspendWhenHidden = FALSE)
  
  # Output for showing completion status
  output$singleR_completed <- reactive({
    singleR_values$completed
  })
  outputOptions(output, "singleR_completed", suspendWhenHidden = FALSE)
  
  # Source the auto naming server functions
  autoNameServer(input, output, session, rv, singleR_values)
}