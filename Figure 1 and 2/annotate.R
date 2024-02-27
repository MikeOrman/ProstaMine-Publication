#input = dataframe, where column 1 = Hugo Names, and column 2 = quantification metric
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
annotate <- function(input){
  main <- input
  rownames(main) <- 1:nrow(main)
  colnames(main) <- c("Hugo_Symbol", "value")
  #Annotate Entrez Gene ID
  mapped_HGNC <- as.data.frame(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
  #Map Entrez Gene ID to HGNC Symbols in main dataframe
  for (i in 1:nrow(main)){
    index <- match(main$Hugo_Symbol[i], mapped_HGNC$symbol, nomatch = NA)
    main[i,3] <- mapped_HGNC$gene_id[index]
  }
  colnames(main)[3] <- "Gene ID"
  #Annotate Chromosome Number
  mapped_CHR <- as.data.frame(org.Hs.egCHR[mappedkeys(org.Hs.egCHR)])
  #Map Entrez Gene ID to Chromosome Number in main dataframe
  for (i in 1:nrow(main)){
    index <- match(main$`Gene ID`[i], mapped_CHR$gene_id, nomatch = NA)
    main[i,4] <- mapped_CHR$chromosome[index]
  }
  colnames(main)[4] <- "Chr"
  #Annotate Starting Position
  #This for loop is adjusted to handle multiple start site entries under the same Gene ID. It takes the start site having
  #a matching Gene ID and chromosome name
  mapped_CHRLOC <- as.data.frame(org.Hs.egCHRLOC[mappedkeys(org.Hs.egCHRLOC)])
  mapped_CHRLOC$start_location <- abs(mapped_CHRLOC$start_location)
  for (i in 1:nrow(main)){
    index <- which(main$`Gene ID`[i] == mapped_CHRLOC$gene_id)
    if (length(index) == 1) {main[i,5] <- mapped_CHRLOC$start_location[index]}
    if (length(index) > 1){
      start_locations <- mapped_CHRLOC[index,]
      index2 <- match(main$Chr[i], start_locations$Chromosome)
      main[i,5] <- start_locations$start_location[index2]
    }
  }
  colnames(main)[5] <- "Starting Position"
  #Annotate Ending Position
  #This for loop is adjusted to handle multiple end site entries under the same Gene ID. It takes the end site having
  #a matching Gene ID and chromosome name
  mapped_CHRLOCEND <- as.data.frame(org.Hs.egCHRLOCEND[mappedkeys(org.Hs.egCHRLOCEND)])
  mapped_CHRLOCEND$end_location <- abs(mapped_CHRLOCEND$end_location)
  for (i in 1:nrow(main)){
    index <- which(main$`Gene ID`[i] == mapped_CHRLOCEND$gene_id)
    if (length(index) == 1) {main[i,6] <- mapped_CHRLOCEND$end_location[index]}
    if (length(index) > 1){
      end_locations <- mapped_CHRLOCEND[index,]
      index2 <- match(main$Chr[i], end_locations$Chromosome)
      main[i,6] <- end_locations$end_location[index2]
    }
  }
  colnames(main)[6] <- "Ending Position"
  #Annotate Band
  mapped_CHRBAND <- as.data.frame(org.Hs.egMAP[mappedkeys(org.Hs.egMAP)])
  for (i in 1:nrow(main)){
    index <- which(main$`Gene ID`[i] == mapped_CHRBAND$gene_id)
    if (length(index) == 1) {main[i,7] <- mapped_CHRBAND$cytogenetic_location[index]}
    }
  colnames(main)[7] <- "Band"
  return(main)
}