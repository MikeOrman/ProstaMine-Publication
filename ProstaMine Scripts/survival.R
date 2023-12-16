# Input = binary alteration matrix
# mRNA = gene expression matrix
# Input and mRNA features must be harmonized
# clinical.PFS = formatted PFS data from ProstaMine
library(survival)
survival <- function(input, mRNA, clinical.PFS){
  output <- data.frame()
  # For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- rownames(input)[k]
    output[k,1] <- gene
    # Obtain sample names for input
    sample.names <- colnames(input)
    # Create data frame of rows = samples | columns = Patient ID, Group, Time, Status
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- match(sample.names[i], clinical.PFS$`Sample ID`)
      data[i,1] <- clinical.PFS$`Patient ID`[index]
      data[i,2] <- clinical.PFS$Time[index]
      data[i,3] <- clinical.PFS$Status[index]
      gene.index <- match(gene, rownames(mRNA))
      sample.index <- match(sample.names[i], colnames(mRNA))
      data[i,4] <- mRNA[gene.index, sample.index]
    }
    data <- na.omit(data)
    colnames(data) <- c("Patient ID", "PFS Months", "PFS Status", "Expression")
    # Order dataset by gene expression. Subset into high and low expressors
    data <- data[order(data$Expression, decreasing = TRUE),]
    # Stratify patients into n=4 quantiles.
    quantiles <- quantile(data$Expression, 0.5)
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,5] <- "High Expressers"
    quantiles <- quantile(data$Expression, 0.5)
    colnames(high.expressing)[5] <- "Group"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,5] <- "Low Expressers"
    colnames(low.expressing)[5] <- "Group"
    # Prepare dataframe for survival analysis
    high.expressing <- high.expressing[,c(1, 5, 2, 3)]
    low.expressing <- low.expressing[,c(1, 5, 2, 3)]
    surv_data <- rbind(high.expressing, low.expressing)
    surv_data$`PFS Months` <- as.numeric(surv_data$`PFS Months`)
    surv_data$`PFS Status` <- as.numeric(surv_data$`PFS Status`)
    # Compute survival statistics
    surv_object <- Surv(time = surv_data$`PFS Months`, event = surv_data$`PFS Status`)
    test <- survdiff(surv_object ~ Group, data = surv_data, rho = 1)
    chisq <- test$chisq
    pval <- pchisq(chisq, 1, lower.tail = FALSE)
    OR <- (test$obs[1]*test$exp[2]) / (test$obs[2]*test$exp[1])
    output[k,2] <- OR
    output[k,3] <- pval
  }
  output$`fdr` <- p.adjust(output[,3], method = "fdr")
  colnames(output) <- c("Hugo Symbol", "Odds Ratio (High Expressers:Low Expressers)", "PFS pval", "PFS FDR")
  return(output)
}