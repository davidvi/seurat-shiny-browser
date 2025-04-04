cluster_tab <- function() {
  tabPanel(
    title = "Normalize/Cluster",
    fluidRow(
      column(4,
        wellPanel(
          h3("Step 1: Quality Control"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            sliderInput(
              inputId = "min_features",
              label = "Min Features per Cell:",
              min = 50, 
              max = 5000, 
              value = 200, 
              step = 50
            ),
            sliderInput(
              inputId = "max_features",
              label = "Max Features per Cell:",
              min = 500, 
              max = 10000, 
              value = 2500, 
              step = 100
            ),
            sliderInput(
              inputId = "max_mt_percent",
              label = "Max Mitochondrial %:",
              min = 0, 
              max = 50, 
              value = 5, 
              step = 1
            ),
            actionButton(
              inputId = "run_qc",
              label = "1. Run Quality Control Filtering",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 2: Normalization"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            selectInput(
              inputId = "normalization_method",
              label = "Normalization Method:",
              choices = c(
                "LogNormalize" = "LogNormalize",
                "SCTransform" = "SCTransform"
              ),
              selected = "LogNormalize"
            ),
            numericInput(
              inputId = "scale_factor",
              label = "Scale Factor:",
              value = 10000,
              min = 1000,
              max = 100000,
              step = 1000
            ),
            actionButton(
              inputId = "run_normalization",
              label = "2. Run Normalization",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 3: Find Variable Features"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            numericInput(
              inputId = "n_features",
              label = "Number of Variable Features:",
              value = 2000,
              min = 500,
              max = 5000,
              step = 100
            ),
            actionButton(
              inputId = "run_variable_features",
              label = "3. Find Variable Features",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 4: Scale Data"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            checkboxInput(
              inputId = "scale_all_genes",
              label = "Scale All Genes (slower but more accurate)",
              value = FALSE
            ),
            actionButton(
              inputId = "run_scaling",
              label = "4. Scale Data",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 5: Principal Component Analysis"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            numericInput(
              inputId = "npcs",
              label = "Number of PCs to Compute:",
              value = 30,
              min = 10,
              max = 100,
              step = 5
            ),
            actionButton(
              inputId = "run_pca",
              label = "5. Run PCA",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 6: Find Neighbors"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            sliderInput(
              inputId = "pcs_for_clustering",
              label = "PCs to Use for Clustering:",
              min = 1,
              max = 50,
              value = c(1, 10),
              step = 1
            ),
            actionButton(
              inputId = "run_neighbors",
              label = "6. Find Neighbors",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 7: Find Clusters"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            selectInput(
              inputId = "cluster_algorithm",
              label = "Clustering Algorithm:",
              choices = c(
                "Original Louvain" = 1,
                "Louvain with Multilevel Refinement" = 2,
                "SLM Algorithm" = 3,
                "Leiden Algorithm" = 4
              ),
              selected = 1
            ),
            sliderInput(
              inputId = "resolution",
              label = "Resolution (Higher = More Clusters):",
              min = 0.1,
              max = 2.0,
              value = 0.8,
              step = 0.1
            ),
            selectInput(
              inputId = "modularity_function",
              label = "Modularity Function:",
              choices = c(
                "Standard" = 1,
                "Alternative" = 2
              ),
              selected = 1
            ),
            numericInput(
              inputId = "cluster_n_start",
              label = "Number of Random Starts:",
              value = 10,
              min = 1,
              max = 50,
              step = 1
            ),
            numericInput(
              inputId = "cluster_n_iter",
              label = "Max Iterations per Start:",
              value = 10,
              min = 1,
              max = 100,
              step = 1
            ),
            checkboxInput(
              inputId = "group_singletons",
              label = "Group Singletons into Nearest Cluster",
              value = TRUE
            ),
            actionButton(
              inputId = "run_clustering",
              label = "7. Find Clusters",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          h3("Step 8: Run UMAP"),
          div(
            style = "border-left: 4px solid #ddd; padding-left: 10px; margin-bottom: 15px;",
            numericInput(
              inputId = "umap_mindist",
              label = "Min Distance:",
              value = 0.3,
              min = 0.01,
              max = 0.99,
              step = 0.05
            ),
            numericInput(
              inputId = "umap_spread",
              label = "Spread:",
              value = 1.0,
              min = 0.1,
              max = 5.0,
              step = 0.1
            ),
            actionButton(
              inputId = "run_umap",
              label = "8. Run UMAP",
              class = "btn-info",
              style = "width: 100%;"
            )
          ),
          
          # Run all steps button
          div(
            style = "margin-top: 20px;",
            actionButton(
              inputId = "run_all_steps",
              label = "Run All Steps",
              class = "btn-primary",
              style = "width: 100%;"
            )
          )
        )
      ),
      column(8,
        tabsetPanel(
          id = "clustering_results_tabs",
          
          # QC metrics visualization
          tabPanel(
            title = "QC Metrics",
            plotOutput("qc_violin_plot", height = "300px"),
            plotOutput("feature_scatter_plot", height = "300px")
          ),
          
          # Variable features visualization
          tabPanel(
            title = "Variable Features",
            plotOutput("var_features_plot", height = "400px")
          ),
          
          # Elbow plot for PCA
          tabPanel(
            title = "PCA",
            plotOutput("elbow_plot", height = "300px"),
            plotOutput("pca_heatmap", height = "400px")
          ),
          
          # Clustering results
          tabPanel(
            title = "Clusters",
            plotOutput("umap_plot", height = "500px"),
            downloadButton("download_umap", "Download UMAP Plot"),
            hr(),
            tableOutput("cluster_stats")
          )
        )
      )
    )
  )
}