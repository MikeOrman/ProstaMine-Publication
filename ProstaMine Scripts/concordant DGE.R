# dataset = binary alteration matrix
# mRNA = mRNA data
# Note: dataset and mRNA must have identical features
#--------------------------------------COMPUTE DGE FDR--------------
DGE <- function(dataset, mRNA){
  DGE <- data.frame()
  for (i in 1:nrow(mRNA)){
    gene <- rownames(mRNA)[i]
    DGE[i,1] <- gene
    # Subset patients in dataset by alteration of gene mRNA data frame
    altered <- alteration_subset(dataset, c(gene, "altered"))
    unaltered <- alteration_subset(dataset, c(gene, "unaltered"))
    altered.patients <- colnames(altered)
    unaltered.patients <- colnames(unaltered)
    altered.GEX.patients <- intersect(altered.patients, colnames(mRNA))
    altered.GEX <- mRNA[i, match(altered.GEX.patients, colnames(mRNA))]
    unaltered.GEX.patients <- intersect(unaltered.patients, colnames(mRNA))
    unaltered.GEX <- mRNA[i,match(unaltered.GEX.patients, colnames(mRNA))]
    #---------------------------COMPUTE DGE BY STUDENT'S T-TEST-----------------
    # Check if at least 3 samples in each group
    if (length(altered.GEX.patients) >= 3 & length(unaltered.GEX.patients) >= 3 &
        class(altered.GEX.patients) == "character" & class(unaltered.GEX.patients) == "character") {
      if ((abs(var(altered.GEX)) > 0) & abs(var(unaltered.GEX)) > 0) {
        ST.expression <- mean(as.numeric(altered.GEX))
        WT.expression <- mean(as.numeric(unaltered.GEX))
        # Compute FC
        FC <- ST.expression - WT.expression
        DGE[i,2] <- FC
        # Compute pval
        DGE[i,3] <- t.test(as.numeric(altered.GEX), (as.numeric(unaltered.GEX)))$p.value
      }
      }
    else {
      DGE[i,2] <- NA
      DGE[i,3] <- NA
    }
  }
  DGE$`FDR` <- p.adjust(DGE[,3], method = "fdr")
  colnames(DGE) <- c("Hugo Symbol", "FC (altered:unaltered)", "p-val", "FDR")
  return(DGE)
}
#----------ADD DGE STATS TO SUMMARY TABLE | FILTER FOR CONCORDANCE--------------
# Input = output from the above function
# Function references summary.enrichment table from ProstaMine
concordant_DGE <- function(input){
# Initialize concordant DGE data frame
  concordant.DGE <- data.frame()
  for (i in 1:nrow(summary.enrichment)){
    concordant.DGE[i,1] <- summary.enrichment$`Hugo Symbol`[i]
    concordant.DGE[i,2] <- summary.enrichment$Alteration[i]
    # Find index in DGE
    index <- match(concordant.DGE[i,1], input[,1])
    if (concordant.DGE[i,2] == "LOF" & input[index,2] < 0 & is.na(input[index,2]) == FALSE){
      concordant.DGE[i,3] <- input[index,2]
      concordant.DGE[i,4] <- input[index,4]
    } 
    if (concordant.DGE[i,2] == "GOF" & input[index,2] > 0 & is.na(input[index,2]) == FALSE){
      concordant.DGE[i,3] <- input[index,2]
      concordant.DGE[i,4] <- input[index,4]
    }
  }
  colnames(concordant.DGE) <- c("Hugo Symbol", "Alteration", "Concordant FC (altered:unaltered)", "FDR")
  return(concordant.DGE)
  }