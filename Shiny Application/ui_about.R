shiny_panel_about<- fluidPage(
  h2("ProstaMine identifies subtype-specific co-alterations associated with aggressiveness in prostate cancer."),
  h3("How to use ProstaMine:"),
  h4("1. Start by selecting an input molecular subtype from the table on the `Select Subtype` tab."),
  h4("2. Observe ProstaMine hits in the table on the `Analyze Hits` tab and apply desired filtering criteria."),
  h4("3. Click on a hit. Observe the rank and location of the hit right under the table on the `Analyze Hits` tab."),
  h4("4. Naviagte to the `Integrated Analysis` tab to observe genomic, transcriptomic, and clinical analysis."),
  h3("ProstaMine Workflow:"),
  tags$figure(
    class = "centerFigure",
    tags$img(
      src = "Figure 3.jpg",
      width = 800,
    ),
  )
)