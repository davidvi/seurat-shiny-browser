# UI for the auto name tab using SingleR for cell type annotation

auto_name_tab <- function() {
  tabPanel(
    title = "Auto Name",
    fluidRow(
      column(
        width = 12,
        wellPanel(
          h4("Automatic Cell Type Annotation"),
          p("This tab uses SingleR to automatically annotate cell clusters based on either built-in reference datasets or another Seurat sample."),
          
          tags$div(class = "alert alert-info",
            tags$p("You can either use built-in reference databases or another annotated Seurat object as a reference."),
            tags$p("When using another sample, make sure it has cell type annotations in its metadata.")
          ),
          
          fluidRow(
            column(
              width = 4,
              selectInput(
                "reference_type",
                "Reference Type:",
                choices = c(
                  "Built-in Reference Database" = "builtin",
                  "Use Another Sample as Reference" = "sample"
                ),
                selected = "builtin"
              ),
              
              conditionalPanel(
                condition = "input.reference_type == 'builtin'",
                selectInput(
                  "reference_dataset",
                  "Built-in Reference Dataset:",
                  choices = c(
                    "Human Primary Cell Atlas" = "hpca",
                    "Blueprint/ENCODE" = "blueprint_encode",
                    "Mouse RNA-seq" = "mouse_rnaseq"
                  ),
                  selected = "hpca"
                )
              ),
              
              conditionalPanel(
                condition = "input.reference_type == 'sample'",
                selectInput(
                  "reference_folder",
                  "Reference Sample Folder:",
                  choices = c("Loading folders..." = ""),
                  selected = ""
                ),
                selectInput(
                  "reference_sample",
                  "Reference Sample:",
                  choices = c("Select folder first" = ""),
                  selected = ""
                ),
                selectInput(
                  "reference_column",
                  "Reference Label Column:",
                  choices = c("Select sample first" = ""),
                  selected = ""
                )
              )
            ),
            
            column(
              width = 4,
              conditionalPanel(
                condition = "input.reference_type == 'builtin'",
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
              
              selectInput(
                "de_method",
                "Marker Detection Method:",
                choices = c(
                  "Classic" = "classic",
                  "Wilcoxon Test (for single-cell)" = "wilcox"
                ),
                selected = "classic"
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
            fluidRow(
              column(
                width = 12,
                actionButton(
                  "transfer_labels",
                  "Transfer Cell Type Labels to Active Idents",
                  icon = icon("exchange-alt"),
                  class = "btn-success",
                  style = "width: 100%; margin-bottom: 10px;"
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                actionButton(
                  "transfer_prefix_labels",
                  "Transfer Cluster+Cell Type to Active Idents",
                  icon = icon("tags"),
                  class = "btn-info",
                  style = "width: 100%; margin-bottom: 10px;"
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                actionButton(
                  "reset_clusters",
                  "Reset to Original Clusters",
                  icon = icon("undo"),
                  class = "btn-warning",
                  style = "width: 100%;"
                )
              )
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