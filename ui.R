# Load tab modules
source("tabs/load_save_tab.R")
source("tabs/visualize_tab.R")
source("tabs/edit_tab.R")
source("tabs/analyze_tab.R")
source("tabs/merge_tab.R")
source("tabs/cluster_tab.R")

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
  h4(textOutput("sample_name")),
  tabsetPanel(
    id = "tabs",
    load_save_tab(),
    visualize_tab(),
    edit_tab(),
    analyze_tab(),
    merge_tab(),
    cluster_tab()
  )
)
