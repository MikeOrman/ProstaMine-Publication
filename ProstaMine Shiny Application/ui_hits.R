shiny_panel_hits <- fluidPage(
  h2(textOutput("Subtype")),
  h3("Filter Hits"),
  fluidRow(
    column(3, sliderInput("cutoff1", h3("Primary Co-Alteration Frequency Difference"), min = 0.02, max = .2, value = 0.02),
           sliderInput("cutoff2", h3("Primary Co-Alteration FDR"), min = 0, max = 1, value = 1)),
    column(3, sliderInput("cutoff3", h3("Metastatic Co-Alteration Frequency Difference"), min = 0.02, max = .2, value = 0.02),
           sliderInput("cutoff4", h3("Metastatic Co-Alteration FDR"), min = 0, max = 1, value = 1)),
    column(3, sliderInput("cutoff5", h3("Primary DGE FDR"), min = 0, max = 1, value = 1),
           sliderInput("cutoff6", h3("Metastatic DGE FDR"), min = 0, max = 1, value = 1)),
    column(3, sliderInput("cutoff7", h3("Survival p-val"), min = 0, max = .2, value = 0.2))),
  # tags$hr(),
  h3("Select Hit"),
  DT::dataTableOutput("hits_table"),
  hr(),
  # ProstaMine Score
  h3("Rank"),
  plotOutput("rank", width = "50%"),
  hr(),
  # Karyoplot
  h3("Location"),
  plotOutput("karyoplot", width = "50%"),
)