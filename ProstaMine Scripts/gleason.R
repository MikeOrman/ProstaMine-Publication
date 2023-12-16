# Computes enrichment of high-risk gleason scores in primary patients stratified by gene expression
# Genes in binary matrix and genes in mRNA data must be harmonized
# Input = binary alteration matrix
# mRNA = gene expression matrix
gleason <- function(input, mRNA, clinical.sample){
  output <- data.frame()
# For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- rownames(input)[k]
    output[k,1] <- gene
    # Obtain sample names
    sample.names <- colnames(input)
    # Create data frame of rows = samples | columns = gleason score, risk group, Gene A expression
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- match(sample.names[i], clinical.sample$`Sample ID`)
      if (is.na(index) == FALSE & is.na(clinical.sample$`Grade Group`[index]) == FALSE
          & length(clinical.sample$`Grade Group`[index]) > 0){
        data[i,1] <- clinical.sample$`Grade Group`[index]
        if (data[i,1] == ">=8") {data[i,2] <- "High"}
        if (data[i,1] == "4+3" ) {data[i,2] <- "Intermediate"}
        if (data[i,1] == "3+4" ) {data[i,2] <- "Intermediate"}
        if (data[i,1] <= "<=6") {data[i,2] <- "Low"}
        gene.index <- match(gene, rownames(mRNA))
        sample.index <- match(sample.names[i], colnames(mRNA))
        data[i,3] <- mRNA[gene.index, sample.index]
      }
    }
    data <- na.omit(data)
    colnames(data) <- c("Score", "Risk Group", "Expression")
# Order dataset by gene expression. Subset into high and low expressers.
    data <- data[order(data[,3], decreasing = TRUE),]
    # Stratify patients into n=4 quantiles.
    quantiles <- quantile(data$Expression, c(0.5))
    high.expressing <- data[data$Expression >= quantiles[1],]
    quantiles <- quantile(data$Expression, c(0.5))
    low.expressing <- data[data$Expression <= quantiles[1],]
# Compute high-risk enrichment
    # Construct contingency table
    high.expressing.high.risk <- nrow(high.expressing[high.expressing$`Risk Group` == "High",])
    high.expressing.not.high.risk <- nrow(high.expressing) - high.expressing.high.risk
    low.expressing.high.risk <- nrow(low.expressing[low.expressing$`Risk Group` == "High",])
    low.expressing.not.high.risk <- nrow(low.expressing) - low.expressing.high.risk
    contignecy.table <- matrix(data = c(high.expressing.high.risk, high.expressing.not.high.risk,
                                        low.expressing.high.risk, low.expressing.not.high.risk), 
                              ncol = 2, nrow = 2, byrow = TRUE)
    # Compute high-risk gleason enrichment between high and low expressers
    fishers.test <- fisher.test(contignecy.table)
    pval <- fishers.test$p.value
# Add high risk OR and pval to summary table
    output[k,2] <- fishers.test$estimate
    output[k,3] <- pval
  }
  output$`fdr` <- p.adjust(output[,3], method = "fdr")
  colnames(output) <- c("Hugo Symbol", "Odds Ratio", "Tumor Grade pval", "Tumor Grade FDR")
  return(output)
}