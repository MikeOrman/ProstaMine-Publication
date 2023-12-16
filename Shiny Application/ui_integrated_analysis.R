shiny_panel_integrated <- fluidPage(
  h3("Genomic & Transcriptomic Analysis"),
  fluidRow(column(5, plotOutput("genomic_plot", height = 500)), column(7, plotOutput("gex_plot", height = 500))),
  h3("Gleason Grade Analysis"),
  fluidRow(column(6, plotOutput("gleason_wt_plot", height = 475)), column(6, plotOutput("gleason_st_plot", height = 475))),
  h3("Survival Analysis"),
  fluidRow(column(6, plotOutput("survival_wt_plot", height = 550)), column(6, plotOutput("survival_st_plot", height = 550)))
)