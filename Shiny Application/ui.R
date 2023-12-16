library(shiny)
library(shinydashboard)
library(shinyalert)
library(shinyjs)
source("ui_about.R", local = TRUE)
source("ui_input.R",local = TRUE)
source("ui_hits.R",local = TRUE)
#source("ui_oncoprints.R", local = TRUE)
source("ui_integrated_analysis.R",local = TRUE)

shinyUI(fluidPage(
  shinybusy::add_busy_spinner(),
  shinyjs::useShinyjs(),
  dashboardPage(
    dashboardHeader(title = "ProstaMine"),
    dashboardSidebar(disable = TRUE),
    dashboardBody(fluidRow(tabBox(
      tabPanel("About ProstaMine", shiny_panel_about),
      tabPanel("Select Subtype", shiny_panel_input),
      tabPanel("Analyze Hits", shiny_panel_hits),
      #tabPanel("Oncoprints", shiny_panel_oncoprints),
      tabPanel("Integrative Analysis", shiny_panel_integrated),
      width=12))))
))