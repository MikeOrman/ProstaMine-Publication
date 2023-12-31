# Computes alteration enrichment between two patient groups (Set 1 and Set 2). 
# Set 1 and Set 2 must be a binary matrix of alterations. Rows = genes | Col 1 = Hugo Symbol | Col 2 = sample name
source("alteration rank.R")
alteration.enrichment <- function(set1, set2){
  if (nrow(set1) == nrow(set2)){
    output <- data.frame()
    for (i in 1 :nrow(set1)){
      gene.set1 <- set1[i,]
      set1.altered <- sum(gene.set1 > 0)
      set1.unaltered <- sum(gene.set1 == 0)
      gene.set2 <- set2[i,]
      set2.altered <- sum(gene.set2 > 0)
      set2.unaltered <- sum(gene.set2 == 0)
      output[i,1] <- rownames(set1)[i]
      output[i,2] <- set1.altered
      output[i,3] <- set1.unaltered
      output[i,4] <- set2.altered
      output[i,5] <- set2.unaltered
      if (is.na(set1.altered) | is.na(set1.unaltered) | is.na(set2.altered) | is.na(set2.unaltered)) {
        output[i,6] <- NA
        output[i,7] <- NA
        output[i,8] <- NA
      }
      else{
        contingency.table <- matrix(data = c(set1.altered, set1.unaltered, set2.altered, set2.unaltered),
                                    nrow = 2, ncol = 2, byrow = TRUE)
        fishers.test <- fisher.test(contingency.table)
        output[i,6] <- fishers.test$estimate
        output[i,7] <- fishers.test$p.value}
    }
    output[,8] <- p.adjust(output[,7], method = "fdr")
    colnames(output) <- c("Hugo Symbol", "Set 1 altered count", "Set 1 unaltered count",
                        "Set 2 altered count", "Set 2 unaltered count", "Odds Ratio", 
                        "Enrichment p-val", "Enrichment FDR")
    return(output)
  }
  if (nrow(set1) != nrow(set2)) {
    return(print("Datasets do not have identical Column 1. Cannot compute enrichment"))
  }
}