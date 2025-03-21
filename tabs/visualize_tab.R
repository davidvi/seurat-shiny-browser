source("tabs/visualization_module.R")

visualize_tab <- function() {
  tabPanel(
    title = "Visualize",
    fluidRow(
      column(4,
        wellPanel(
          h3("Gene Selection"),
          # Single gene search
          h4("Search for a Single Gene"),
          div(
            style = "display: flex; align-items: flex-end;",
            div(
              style = "flex-grow: 1; margin-right: 10px;",
              textInput(
                inputId = "gene",
                label = "Gene Name:",
                value = "CSF1R"
              )
            ),
            div(
              actionButton("search_button", "Search", class = "btn-primary")
            )
          ),
          hr(),
          
          # Multiple gene selection
          h4("Multiple Gene Selection"),
          div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 4px; border-left: 4px solid #28a745;",
            p("Select multiple genes for dotplot visualization"),
            textAreaInput(
              inputId = "multi_gene_input",
              label = "Enter genes (one per line or comma-separated):",
              height = "100px",
              placeholder = "e.g., MS4A1, CD79A, CD79B"
            ),
            actionButton("multi_gene_search_button", "Display Selected Genes", class = "btn-success")
          )
        ),
        wellPanel(
          h3("Display Options"),
          selectInput(
            inputId = "reduction_selector",
            label = "Dimensional Reduction Method:",
            choices = available_reductions,
            selected = available_reductions[1]
          ),
          selectInput(
            inputId = "metadata_column_selector",
            label = "Color cells by:",
            choices = c("ident")
          )
        )
      ),
      column(8,
        # Use the modular visualization component with prefix "visualize"
        visualization_module("visualize")
      )
    )
  )
}