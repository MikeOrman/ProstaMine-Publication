shiny_panel_input <- fluidPage(
  h3("Select Background Subtype"),
  DT::dataTableOutput("Subtype_table")
)