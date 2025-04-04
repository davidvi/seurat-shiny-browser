# UI for the auto name tab using SingleR for cell type annotation

auto_name_tab <- function() {
  tabPanel(
    title = "Auto Name",
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Automatic Cell Type Annotation"),
          p("This tab uses SingleR to automatically annotate cell clusters based on reference datasets."),
          
          fluidRow(
            column(
              width = 4,
              selectInput(
                "reference_dataset",
                "Reference Dataset:",
                choices = c(
                  "Human Primary Cell Atlas" = "hpca",
                  "Blueprint/ENCODE" = "blueprint_encode",
                  "Mouse RNA-seq" = "mouse_rnaseq"
                ),
                selected = "hpca"
              )
            ),
            
            column(
              width = 4,
              selectInput(
                "label_type",
                "Label Type:",
                choices = c(
                  "Main Labels" = "main",
                  "Fine Labels" = "fine"
                ),
                selected = "main"
              )
            ),
            
            column(
              width = 4,
              actionButton(
                "run_singleR",
                "Run SingleR",
                icon = icon("play"),
                class = "btn-primary"
              )
            )
          ),
          
          tags$hr(),
          
          conditionalPanel(
            condition = "output.singleR_completed",
            actionButton(
              "transfer_labels",
              "Transfer Labels to Active Idents",
              icon = icon("exchange-alt"),
              class = "btn-success",
              style = "width: 100%;"
            )
          )
        ),
        
        # Status panels
        conditionalPanel(
          condition = "output.singleR_running",
          wellPanel(
            div(
              class = "text-center",
              h4("Running SingleR annotation..."),
              p("This may take a few minutes depending on the size of your dataset.")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.singleR_completed",
          wellPanel(
            h4("SingleR Results"),
            p("The table below shows the mapping between clusters and cell types assigned by SingleR."),
            p("Click 'Transfer Labels to Active Idents' to set these labels as the active identities in your Seurat object."),
            DT::dataTableOutput("annotation_table")
          )
        ),
        
        conditionalPanel(
          condition = "output.singleR_error",
          wellPanel(
            class = "bg-danger",
            h4("Error in SingleR Annotation"),
            verbatimTextOutput("error_text"),
            p("Possible solutions:"),
            tags$ul(
              tags$li("Make sure the SingleR and celldex packages are installed"),
              tags$li("Try a different reference dataset"),
              tags$li("Check if your Seurat object has normalized data")
            )
          )
        )
      )
    )
  )
}