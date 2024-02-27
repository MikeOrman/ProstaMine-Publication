# Input = hits table from ProtstaMine Algorithm
prostamine_prioritization <- function(input){
  Score <- c()
  Coalteration_order <- input[order(input$`Primary Alteration Ratio Difference`, decreasing = FALSE),]
  Metastatic_order <- input[order(input$`Metastatic Alteration Ratio Difference`, decreasing = FALSE),]
  Surv_order <- input[order(input$`TCGA ST PFS pval`, decreasing = TRUE),]
  for (i in 1:nrow(input)){
    Coalteration_Rank <- match(input$Gene[i], Coalteration_order$Gene)
    Metastatic_Rank <- match(input$Gene[i], Metastatic_order$Gene)
    Surv_Rank <- match(input$Gene[i], Surv_order$Gene)
    Score[i] <- ((Metastatic_Rank/nrow(input))*(.3) +
                   (Coalteration_Rank/nrow(input))*(.3) +
                   (Surv_Rank/nrow(input))*(.4))
  }
  input$`ProstaMine Score` <- Score
  input <- input[order(input$`ProstaMine Score`, decreasing = TRUE),]
  rownames(input) <- 1:nrow(input)
  return(input)
}