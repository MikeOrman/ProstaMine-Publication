#----Get Input List----
 loss.input <- read.table("Loss_Input_1.txt")
 loss.input <- loss.input$V1
 # Make .mae
 library(curatedPCaData)
 TCGA.mae <- getPCa("tcga")
 Taylor.mae <- getPCa("taylor")
 Barbieri.mae <- getPCa("barbieri")
 Abida.mae <- getPCa("abida")
#----Run ProstaMine for Loss Subtypes----
for (z in 1 :length(loss.input)){
  source("analysis scripts.R")
  #----Specify subtype and select filtering cutoffs----
  query.gene <- loss.input[z]
  # Genomic Analysis Filters
  Primary.Alteration.FDR.cutoff <- 1
  Met.Alteration.FDR.cutoff <- 1
  Primary.Alteration.Rate.cutoff <- 0.02
  Met.Alteration.Rate.cutoff <- 0.02
  # Concordant Gene Expression Filter
  DGE.FDR.cutoff <- 1
  # Clinical Analysis Filter
  ST.cutoff <- .2
  WT.cutoff <- .3
  # Output Labels
  query.gene.alteration <- paste(query.gene, "Loss")
  query.gene.wt <- paste(query.gene, "WT")
  wt <- paste(query.gene,"WT Tumors")
  st <- paste(query.gene, "Loss Tumors")
  output.filename <- paste(query.gene, "Loss.RData", sep = "")
  #----Filter coalterations----
  # Load processed data
  primary.alteration <- read.table("primary alteration.txt", header = TRUE, check.names = FALSE, sep = "\t")
  primary.LOF <- read.table("primary LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  primary.GOF <- read.table("primary GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  met.alteration <- read.table("met alteration.txt", header = TRUE, check.names = FALSE, sep = "\t")
  met.LOF <- read.table("met LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  met.GOF <- read.table("met GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  TCGA.alteration <- read.table("TCGA alteration.txt", header = TRUE, check.names = FALSE, sep = "\t")
  TCGA.LOF <- read.table("TCGA LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  TCGA.GOF <- read.table("TCGA GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Taylor.primary.alteration <- read.table("Taylor primary alteration.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Taylor.primary.LOF <- read.table("Taylor primary LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Taylor.primary.GOF <- read.table("Taylor primary GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Taylor.met.LOF <- read.table("Taylor met LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Taylor.met.GOF <- read.table("Taylor met GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Barbieri.alteration <- read.table("Barbieri alteration.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Barbieri.LOF <- read.table("Barbieri LOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  Barbieri.GOF <- read.table("Barbieri GOF.txt", header = TRUE, check.names = FALSE, sep = "\t")
  # If the user chooses "LOF" as the alteration, then primary.LOF, met.LOF, TCGA.LOF, Taylor.LOF and Barbieri.LOF are the input datasets.
  # If the user chooses "GOF" as the alteration, then primary.GOF, met.GOF, TCGA.GOF, Taylor.GOF, and Barbieri.GOF are the input datasets.
  # Primary
  primary.ST <- alteration_subset(primary.LOF, c(query.gene, "altered"))
  primary.ST.patients <- colnames(primary.ST)
  primary.WT <- alteration_subset(primary.alteration, c(query.gene, "unaltered"))
  primary.WT.patients <- colnames(primary.WT)
  # Metastatic
  met.ST <- alteration_subset(met.LOF, c(query.gene, "altered"))
  met.ST.patients <- colnames(met.ST)
  met.WT <- alteration_subset(met.alteration, c(query.gene, "unaltered"))
  met.WT.patients <- colnames(met.WT)
  # TCGA background for DGE
  TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
  TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
  # Taylor primary background for DGE
  Taylor.primary.ST <- alteration_subset(Taylor.primary.LOF, c(query.gene, "altered"))
  Taylor.primary.WT <- alteration_subset(Taylor.primary.LOF, c(query.gene, "unaltered"))
  # Taylor met background for DGE
  Taylor.met.ST <- alteration_subset(Taylor.met.LOF, c(query.gene, "altered"))
  Taylor.met.WT <- alteration_subset(Taylor.met.LOF, c(query.gene, "unaltered"))
  # Barbieri background for DGE
  Barbieri.ST <- alteration_subset(Barbieri.LOF, c(query.gene, "altered"))
  Barbieri.WT <- alteration_subset(Barbieri.alteration, c(query.gene, "unaltered"))
  # Subset primary LOF/GOF alteration profiles by genetic background
  # Primary LOF WT
  primary.LOF.WT.patients <- intersect(colnames(primary.LOF), primary.WT.patients)
  index <- c(TRUE)
  for (i in 1:ncol(primary.LOF)){
    if (sum(colnames(primary.LOF)[i] == primary.LOF.WT.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(primary.LOF)[i] == primary.LOF.WT.patients) == 0) {index[i] = FALSE}
  }
  primary.LOF.WT <- primary.LOF[,index]
  # Primary GOF WT
  primary.GOF.WT.patients <- intersect(colnames(primary.GOF), primary.WT.patients)
  index <- c(TRUE)
  for (i in 1:ncol(primary.GOF)){
    if (sum(colnames(primary.GOF)[i] == primary.GOF.WT.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(primary.GOF)[i] == primary.GOF.WT.patients) == 0) {index[i] = FALSE}
  }
  primary.GOF.WT <- primary.GOF[,index]
  # Primary LOF ST
  primary.LOF.ST.patients <- intersect(colnames(primary.LOF), primary.ST.patients)
  index <- c(TRUE)
  for (i in 1:ncol(primary.LOF)){
    if (sum(colnames(primary.LOF)[i] == primary.LOF.ST.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(primary.LOF)[i] == primary.LOF.ST.patients) == 0) {index[i] = FALSE}
  }
  primary.LOF.ST <- primary.LOF[,index]
  # Primary GOF ST
  primary.GOF.ST.patients <- intersect(colnames(primary.GOF), primary.ST.patients)
  index <- c(TRUE)
  for (i in 1:ncol(primary.GOF)){
    if (sum(colnames(primary.GOF)[i] == primary.GOF.ST.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(primary.GOF)[i] == primary.GOF.ST.patients) == 0) {index[i] = FALSE}
  }
  primary.GOF.ST <- primary.GOF[,index]
  # Subset metastatic LOF/GOF alteration profiles by genetic background
  # Metastatic LOF WT
  met.LOF.WT.patients <- intersect(colnames(met.LOF), met.WT.patients)
  index <- c(TRUE)
  for (i in 1:ncol(met.LOF)){
    if (sum(colnames(met.LOF)[i] == met.LOF.WT.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(met.LOF)[i] == met.LOF.WT.patients) == 0) {index[i] = FALSE}
  }
  met.LOF.WT <- met.LOF[,index]
  # Metastatic GOF WT
  met.GOF.WT.patients <- intersect(colnames(met.GOF), met.WT.patients)
  index <- c(TRUE)
  for (i in 1:ncol(met.GOF)){
    if (sum(colnames(met.GOF)[i] == met.GOF.WT.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(met.GOF)[i] == met.GOF.WT.patients) == 0) {index[i] = FALSE}
  }
  met.GOF.WT <- met.GOF[,index]
  # Metastatic LOF ST
  met.LOF.ST.patients <- intersect(colnames(met.LOF), met.ST.patients)
  index <- c(TRUE)
  for (i in 1:ncol(met.LOF)){
    if (sum(colnames(met.LOF)[i] == met.LOF.ST.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(met.LOF)[i] == met.LOF.ST.patients) == 0) {index[i] = FALSE}
  }
  met.LOF.ST <- met.LOF[,index]
  # Metastatic GOF ST
  met.GOF.ST.patients <- intersect(colnames(met.GOF), met.ST.patients)
  index <- c(TRUE)
  for (i in 1:ncol(met.GOF)){
    if (sum(colnames(met.GOF)[i] == met.GOF.ST.patients) == 1) {index[i] = TRUE}
    if (sum(colnames(met.GOF)[i] == met.GOF.ST.patients) == 0) {index[i] = FALSE}
  }
  met.GOF.ST <- met.GOF[,index]
  # Compute LOF alteration rates in subtype patients
  # Primary
  primary.ST.LOF.frequency <- alterationRank(primary.LOF.ST)
  primary.ST.LOF.frequency <- primary.ST.LOF.frequency[,c(1,5)]
  colnames(primary.ST.LOF.frequency)[2] <- "LOF alteration rate: Primary ST patients"
  # Metastatic
  met.ST.LOF.frequency <- alterationRank(met.LOF.ST)
  met.ST.LOF.frequency <- met.ST.LOF.frequency[,c(1,5)]
  colnames(met.ST.LOF.frequency)[2] <- "LOF alteration rate: Metastatic ST patients"
  # Compute LOF alteration rates in WT patients
  # Primary
  primary.WT.LOF.frequency <- alterationRank(primary.LOF.WT)
  primary.WT.LOF.frequency <- primary.WT.LOF.frequency[,c(1,5)]
  colnames(primary.WT.LOF.frequency)[2] <- "LOF alteration rate: Primary WT patients"
  # Metastatic
  met.WT.LOF.frequency <- alterationRank(met.LOF.WT)
  met.WT.LOF.frequency <- met.WT.LOF.frequency[,c(1,5)]
  colnames(met.WT.LOF.frequency)[2] <- "LOF alteration rate: Metastatic WT patients"
  # Compute GOF alteration rates in subtype patients
  # Primary
  primary.ST.GOF.frequency <- alterationRank(primary.GOF.ST)
  primary.ST.GOF.frequency <- primary.ST.GOF.frequency[,c(1,5)]
  colnames(primary.ST.GOF.frequency)[2] <- "GOF alteration rate: Primary ST patients"
  # Metastatic
  met.ST.GOF.frequency <- alterationRank(met.GOF.ST)
  met.ST.GOF.frequency <- met.ST.GOF.frequency[,c(1,5)]
  colnames(met.ST.GOF.frequency)[2] <- "GOF alteration rate: Metastatic ST patients"
  # Compute GOF alteration rates in WT patients
  # Primary
  primary.WT.GOF.frequency <- alterationRank(primary.GOF.WT)
  primary.WT.GOF.frequency <- primary.WT.GOF.frequency[,c(1,5)]
  colnames(primary.WT.GOF.frequency)[2] <- "GOF alteration rate: Primary WT patients"
  # Metastatic
  met.WT.GOF.frequency <- alterationRank(met.GOF.WT)
  met.WT.GOF.frequency <- met.WT.GOF.frequency[,c(1,5)]
  colnames(met.WT.GOF.frequency)[2] <- "GOF alteration rate: Metastatic WT patients"
  # Create GOF/LOF alteration tables
  LOF.alterations <- merge.data.frame(primary.WT.LOF.frequency, primary.ST.LOF.frequency)
  LOF.alterations <- merge.data.frame(LOF.alterations, met.WT.LOF.frequency)
  LOF.alterations <- merge.data.frame(LOF.alterations, met.ST.LOF.frequency)
  GOF.alterations <- merge.data.frame(primary.WT.GOF.frequency, primary.ST.GOF.frequency)
  GOF.alterations <- merge.data.frame(GOF.alterations, met.WT.GOF.frequency)
  GOF.alterations <- merge.data.frame(GOF.alterations, met.ST.GOF.frequency)
  # Select all primary co-alterations: Alteration Rate Primary ST > Alteration Rate Primary WT
  index.LOF <- (LOF.alterations$`LOF alteration rate: Primary ST patients` > 
                  LOF.alterations$`LOF alteration rate: Primary WT patients`)
  LOF.genes <- LOF.alterations$Hugo_Symbol[index.LOF]
  index.GOF <- (GOF.alterations$`GOF alteration rate: Primary ST patients` > 
                  GOF.alterations$`GOF alteration rate: Primary WT patients`)
  GOF.genes <- GOF.alterations$Hugo_Symbol[index.GOF]
  # Assign alteration status
  filtered.LOF.alterations <- LOF.genes
  filtered.GOF.alterations <- GOF.genes
  altered.genes <- c(filtered.LOF.alterations, filtered.GOF.alterations)
  altered.genes <- altered.genes[duplicated(altered.genes)==FALSE]
  # Apply LOF/GOF gene filter to alteration datasets
  # Filter Primary alteration, ST, and WT
  primary.alteration <- primary.alteration[match(altered.genes, rownames(primary.alteration)),]
  primary.ST <- primary.ST[match(altered.genes, rownames(primary.ST)),]
  primary.WT <- primary.WT[match(altered.genes, rownames(primary.WT)),]
  # Filter Primary LOF WT and ST
  primary.LOF.ST <- primary.LOF.ST[match(altered.genes, rownames(primary.LOF.ST)),]
  primary.LOF.WT <- primary.LOF.WT[match(altered.genes, rownames(primary.LOF.WT)),]
  # Filter Primary GOF WT and ST
  primary.GOF.ST <- primary.GOF.ST[match(altered.genes, rownames(primary.GOF.ST)),]
  primary.GOF.WT <- primary.GOF.WT[match(altered.genes, rownames(primary.GOF.WT)),]
  # Filter Metastatic alteration, ST, and WT
  met.alteration <- met.alteration[match(altered.genes, rownames(met.alteration)),]
  met.ST <- met.ST[match(altered.genes, rownames(met.ST)),]
  met.WT <- met.WT[match(altered.genes, rownames(met.WT)),]
  # Filter Metastatic LOF WT and ST
  met.LOF.ST <- met.LOF.ST[match(altered.genes, rownames(met.LOF.ST)),]
  met.LOF.WT <- met.LOF.WT[match(altered.genes, rownames(met.LOF.WT)),]
  # Filter Metastatic GOF WT and ST
  met.GOF.ST <- met.GOF.ST[match(altered.genes, rownames(met.GOF.ST)),]
  met.GOF.WT <- met.GOF.WT[match(altered.genes, rownames(met.GOF.WT)),]
  # Construct summary table
  # Construct LOF alteration table
  if (length(LOF.genes) > 0){
    LOF.genes <- as.data.frame(LOF.genes)
    LOF.genes[,2] <- "LOF"
    for (i in (1:nrow(LOF.genes))) {
      index1 <- which(primary.WT.LOF.frequency$Hugo_Symbol == LOF.genes[i,1])
      index2 <- which(primary.ST.LOF.frequency$Hugo_Symbol == LOF.genes[i,1])
      index4 <- which(met.ST.LOF.frequency$Hugo_Symbol == LOF.genes[i,1])
      LOF.genes[i,3] <- (primary.ST.LOF.frequency$`LOF alteration rate: Primary ST patients`[index2] 
                         - primary.WT.LOF.frequency$`LOF alteration rate: Primary WT patients`[index1])
      LOF.genes[i,4] <- (met.ST.LOF.frequency$`LOF alteration rate: Metastatic ST patients`[index4] -
                           primary.ST.LOF.frequency$`LOF alteration rate: Primary ST patients`[index2])
    }
    colnames(LOF.genes) <- c("Hugo Symbol", "Alteration", "Primary Co-alteration Ratio Difference", "Subtype Alteration Ratio Difference (Met - Primary)")
  }
  # Construct GOF alteration table
  if (length(GOF.genes) > 0){
    GOF.genes <- as.data.frame(GOF.genes)
    GOF.genes[,2] <- "GOF"
    for (i in (1:nrow(GOF.genes))) {
      index1 <- which(primary.WT.GOF.frequency$Hugo_Symbol == GOF.genes[i,1])
      index2 <- which(primary.ST.GOF.frequency$Hugo_Symbol == GOF.genes[i,1])
      index4 <- which(met.ST.GOF.frequency$Hugo_Symbol == GOF.genes[i,1])
      GOF.genes[i,3] <- (primary.ST.GOF.frequency$`GOF alteration rate: Primary ST patients`[index2] 
                         - primary.WT.GOF.frequency$`GOF alteration rate: Primary WT patients`[index1])
      GOF.genes[i,4] <- (met.ST.GOF.frequency$`GOF alteration rate: Metastatic ST patients`[index4] -
                           primary.ST.GOF.frequency$`GOF alteration rate: Primary ST patients`[index2])
    }
    colnames(GOF.genes) <- c("Hugo Symbol", "Alteration", 
                             "Primary Co-alteration Ratio Difference", 
                             "Subtype Alteration Ratio Difference (Met - Primary)")
  }
  # Combine alteration tables
  if ((length(LOF.genes) > 0) & (length(GOF.genes) > 0)){
    summary <- rbind(LOF.genes, GOF.genes)
  }
  if ((length(LOF.genes) > 0) & (length(GOF.genes) == 0)){
    summary <- (LOF.genes)
  }
  if ((length(LOF.genes) == 0) & (length(GOF.genes) > 0)){
    summary <- (GOF.genes)
  }
  # Clear NAs from summary table. Rows with NAs mean the alteration did not occur in both primary/met settings.
  summary <- na.omit(summary)
  # Alteration Enrichment
  # Compute alteration enrichment between primary WT and primary ST
  LOF.genes <- summary$`Hugo Symbol`[summary$Alteration == "LOF"]
  GOF.genes <- summary$`Hugo Symbol`[summary$Alteration == "GOF"]
  LOF.coalteration.enrichment <- alteration.enrichment(
    primary.LOF.ST[match(LOF.genes, rownames(primary.LOF.ST)),],
    primary.LOF.WT[match(LOF.genes, rownames(primary.LOF.WT)),])
  GOF.coalteration.enrichment <- alteration.enrichment(
    primary.GOF.ST[match(GOF.genes, rownames(primary.GOF.ST)),],
    primary.GOF.WT[match(GOF.genes, rownames(primary.GOF.WT)),])
  Coalteration.enrichment <- na.omit(rbind(LOF.coalteration.enrichment, GOF.coalteration.enrichment))
  # Compute alteration enrichment between primary ST and met ST
  LOF.met.enrichment <- alteration.enrichment(
    met.LOF.ST[match(LOF.genes, rownames(met.LOF.ST)),],
    primary.LOF.ST[match(LOF.genes, rownames(primary.LOF.ST)),])
  GOF.met.enrichment <- alteration.enrichment(
    met.GOF.ST[match(GOF.genes, rownames(met.GOF.ST)),],
    primary.GOF.ST[match(GOF.genes, rownames(primary.GOF.ST)),])
  Metastatic.enrichment <- na.omit(rbind(LOF.met.enrichment, GOF.met.enrichment))
  # Add alteration enrichment to summary table
  summary.enrichment <- data.frame()
  for (i in 1:nrow(summary)){
    summary.enrichment[i,1] <- summary$`Hugo Symbol`[i]
    summary.enrichment[i,2] <- summary$`Alteration`[i]
    summary.enrichment[i,3] <- summary$`Primary Co-alteration Ratio Difference`[i]
    summary.enrichment[i,4] <- Coalteration.enrichment$`Enrichment FDR`[i]
    summary.enrichment[i,5] <- summary$`Subtype Alteration Ratio Difference (Met - Primary)`[i]
    summary.enrichment[i,6] <- Metastatic.enrichment$`Enrichment FDR`[i]
  }
  colnames(summary.enrichment) <- c("Hugo Symbol", "Alteration",
                                    "Primary Alteration Ratio Difference", "Primary Alteration FDR",
                                    "Metastatic Alteration Ratio Difference", "Metastatic Alteration FDR")
  # Filter by alteration rate and FDR
  summary.enrichment <- summary.enrichment[((abs(summary.enrichment$`Primary Alteration Ratio Difference`) >= Primary.Alteration.Rate.cutoff) &
                                            (summary.enrichment$`Primary Alteration FDR` <= Primary.Alteration.FDR.cutoff) &
                                            (abs(summary.enrichment$`Metastatic Alteration Ratio Difference`) >= Met.Alteration.Rate.cutoff) &
                                            (summary.enrichment$`Metastatic Alteration FDR` <= Met.Alteration.FDR.cutoff)),]
  
  #----Gene Expression Analysis----
  met.ST <- alteration_subset(met.LOF, c(query.gene, "altered"))
  # Format gene expression datasets
  # Abida gene expression data
  Abida.mRNA <- Abida.mae@ExperimentList@listData$gex.relz
  # Transform TCGA from rsem to cpm, normalize expression, and filter low-expressing genes
  TCGA.mRNA <- TCGA.mae@ExperimentList@listData$gex.rsem.log
  TCGA.mRNA <- scale(TCGA.mRNA)
  # Taylor gene expression data
  Taylor.mRNA <- Taylor.mae@ExperimentList@listData$gex.rma
  Taylor.mRNA <- scale(Taylor.mRNA)
  # Format Barbieri gene expression data
  Barbieri.mRNA <- Barbieri.mae@ExperimentList@listData$gex.relz
  # Intersect rows and columns of ST alteration profiles and mRNA data
  # If the user chooses "LOF" as the alteration, then Taylor.primary.LOF, Taylor.met.LOF, TCGA.LOF, Barbieri.LOF, and Abida.LOF are the input datasets.
  # If the user chooses "GOF" as the alteration, then Taylor.primary.GOF, Taylor.met.GOF, TCGA.GOF, Barbieri.GOF, and Abida.GOF are the input datasets.
  # Taylor primary
  Taylor.primary.ST <- alteration_subset(Taylor.primary.LOF, c(query.gene, "altered"))
  Taylor.primary.patients <- intersect(colnames(Taylor.primary.ST), colnames(Taylor.mRNA))
  # Taylor metastatic
  Taylor.met.ST <- alteration_subset(Taylor.met.LOF, c(query.gene, "altered"))
  Taylor.met.patients <- intersect(colnames(Taylor.met.ST), colnames(Taylor.mRNA))
  # Barbieri
  Barbieri.ST <- alteration_subset(Barbieri.LOF, c(query.gene, "altered"))
  Barbieri.patients <- intersect(colnames(Barbieri.ST), colnames(Barbieri.mRNA))
  # TCGA
  TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
  TCGA.patients <- intersect(colnames(TCGA.ST), colnames(TCGA.mRNA))
  # Abida
  Abida.ST <- na.omit(alteration_subset(met.LOF, c(query.gene, "altered")))
  Abida.patients <- intersect(colnames(Abida.ST), colnames(Abida.mRNA))
  # Make inputs for concordant DGE
  # All datasets have patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) > 1) & (length(Taylor.met.patients) > 1) & 
      (length(Barbieri.patients) > 1)){
    genes <- intersect(rownames(Taylor.primary.ST), rownames(Taylor.mRNA))
    Taylor.primary.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                       match(Taylor.primary.patients, colnames(Taylor.mRNA))]
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),
                                           match(Taylor.primary.patients, colnames(Taylor.primary.ST))]
    genes <- intersect(rownames(Taylor.met.ST), rownames(Taylor.mRNA))
    Taylor.met.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                   match(Taylor.met.patients, colnames(Taylor.mRNA))]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),
                                   match(Taylor.met.patients, colnames(Taylor.met.ST))]
    genes <- intersect(rownames(Barbieri.ST), rownames(Barbieri.mRNA))
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),
                                   match(Barbieri.patients, colnames(Barbieri.mRNA))]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),
                               match(Barbieri.patients, colnames(Barbieri.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.primary.mRNA), rownames(Taylor.met.mRNA))
    genes <- intersect(genes, rownames(Barbieri.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),]
    Taylor.met.mRNA <- Taylor.met.mRNA[match(genes, rownames(Taylor.met.mRNA)),]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Taylor.primary.mRNA)
    primary.mRNA <- cbind(primary.mRNA, Barbieri.mRNA)
    met.mRNA <- cbind(Abida.mRNA, Taylor.met.mRNA)
  }
  # Barbieri has no patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) > 1) & (length(Taylor.met.patients) > 1) & 
      (length(Barbieri.patients) <= 1)){
    genes <- intersect(rownames(Taylor.primary.ST), rownames(Taylor.mRNA))
    Taylor.primary.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                       match(Taylor.primary.patients, colnames(Taylor.mRNA))]
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),
                                           match(Taylor.primary.patients, colnames(Taylor.primary.ST))]
    genes <- intersect(rownames(Taylor.met.ST), rownames(Taylor.mRNA))
    Taylor.met.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                   match(Taylor.met.patients, colnames(Taylor.mRNA))]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),
                                   match(Taylor.met.patients, colnames(Taylor.met.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.primary.mRNA), rownames(Taylor.met.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),]
    Taylor.met.mRNA <- Taylor.met.mRNA[match(genes, rownames(Taylor.met.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Taylor.primary.mRNA)
    met.mRNA <- cbind(Abida.mRNA, Taylor.met.mRNA)
  }
  # Taylor Primary has no patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) <= 1) & (length(Taylor.met.patients) > 1) & 
      (length(Barbieri.patients) > 1)){
    genes <- intersect(rownames(Taylor.met.ST), rownames(Taylor.mRNA))
    Taylor.met.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                   match(Taylor.met.patients, colnames(Taylor.mRNA))]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),
                                   match(Taylor.met.patients, colnames(Taylor.met.ST))]
    genes <- intersect(rownames(Barbieri.ST), rownames(Barbieri.mRNA))
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),
                                   match(Barbieri.patients, colnames(Barbieri.mRNA))]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),
                               match(Barbieri.patients, colnames(Barbieri.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.met.mRNA), rownames(Barbieri.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),]
    Taylor.met.mRNA <- Taylor.met.mRNA[match(genes, rownames(Taylor.met.mRNA)),]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Barbieri.mRNA)
    met.mRNA <- cbind(Abida.mRNA, Taylor.met.mRNA)
  }
  # Taylor Met has no patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) > 1) & (length(Taylor.met.patients) <= 1) & 
      (length(Barbieri.patients) > 1)){
    genes <- intersect(rownames(Taylor.primary.ST), rownames(Taylor.mRNA))
    Taylor.primary.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                       match(Taylor.primary.patients, colnames(Taylor.mRNA))]
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),
                                           match(Taylor.primary.patients, colnames(Taylor.primary.ST))]
    genes <- intersect(rownames(Barbieri.ST), rownames(Barbieri.mRNA))
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),
                                   match(Barbieri.patients, colnames(Barbieri.mRNA))]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),
                               match(Barbieri.patients, colnames(Barbieri.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.primary.mRNA), rownames(Barbieri.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Taylor.primary.mRNA)
    primary.mRNA <- cbind(primary.mRNA, Barbieri.mRNA)
    met.mRNA <- Abida.mRNA
  }
  # Taylor Primary and Taylor Met are missing patients with match alteration/mRNA data
  if ((length(Taylor.primary.patients) <= 1) & (length(Taylor.met.patients) <= 1) & 
      (length(Barbieri.patients) > 1)){
    genes <- intersect(rownames(Barbieri.ST), rownames(Barbieri.mRNA))
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),
                                   match(Barbieri.patients, colnames(Barbieri.mRNA))]
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),
                               match(Barbieri.patients, colnames(Barbieri.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Barbieri.mRNA), rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Barbieri.mRNA)
    met.mRNA <- Abida.mRNA
  }
  # Taylor Primary and Barbieri have no patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) <= 1) & (length(Taylor.met.patients) > 1) & 
      (length(Barbieri.patients) <= 1)){
    genes <- intersect(rownames(Taylor.met.ST), rownames(Taylor.mRNA))
    Taylor.met.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                   match(Taylor.met.patients, colnames(Taylor.mRNA))]
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),
                                   match(Taylor.met.patients, colnames(Taylor.met.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.met.mRNA), rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.met.ST <- Taylor.met.ST[match(genes, rownames(Taylor.met.ST)),]
    Taylor.met.mRNA <- Taylor.met.mRNA[match(genes, rownames(Taylor.met.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- TCGA.mRNA
    met.mRNA <- cbind(Abida.mRNA, Taylor.met.mRNA)
  }
  # Taylor Met and Barbieri have no patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) > 1) & (length(Taylor.met.patients) <= 1) & 
      (length(Barbieri.patients) <= 1)){
    genes <- intersect(rownames(Taylor.primary.ST), rownames(Taylor.mRNA))
    Taylor.primary.mRNA <- Taylor.mRNA[match(genes, rownames(Taylor.mRNA)),
                                       match(Taylor.primary.patients, colnames(Taylor.mRNA))]
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),
                                           match(Taylor.primary.patients, colnames(Taylor.primary.ST))]
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Taylor.primary.mRNA), rownames(TCGA.mRNA))
    genes <- intersect(genes, rownames(Abida.mRNA))
    # Harmonize across individual mRNA datasets
    Taylor.primary.ST <- Taylor.primary.ST[match(genes, rownames(Taylor.primary.ST)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- cbind(TCGA.mRNA, Taylor.primary.mRNA)
    met.mRNA <- Abida.mRNA
  }
  # Taylor Primary, Taylor Pet, and Barbieri are missing patients with matched alteration/mRNA data
  if ((length(Taylor.primary.patients) <= 1) & (length(Taylor.met.patients) <= 1) & 
      (length(Barbieri.patients) <= 1) & (length(TCGA.patients) > 1) & (length(Abida.patients) > 1)){
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),
                           match(TCGA.patients, colnames(TCGA.mRNA))]
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),
                       match(TCGA.patients, colnames(TCGA.ST))]
    genes <- intersect(rownames(Abida.ST), rownames(Abida.mRNA))
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),
                             match(Abida.patients, colnames(Abida.mRNA))]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),
                         match(Abida.patients, colnames(Abida.ST))]
    # Pre-filter genes for DGE
    genes <- intersect(rownames(Abida.mRNA), rownames(TCGA.mRNA))
    # Harmonize across individual mRNA datasets
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    Abida.ST <- Abida.ST[match(genes, rownames(Abida.ST)),]
    Abida.mRNA <- Abida.mRNA[match(genes, rownames(Abida.mRNA)),]
    # Group primary and met mRNA datasets
    primary.mRNA <- TCGA.mRNA
    met.mRNA <- Abida.mRNA
  }
  # Harmonize combined alteration profiles with combined mRNA
  genes.primary <- intersect(rownames(primary.mRNA), rownames(primary.ST))
  primary.mRNA <- primary.mRNA[match(genes.primary, rownames(primary.mRNA)),]
  primary.ST <- primary.ST[match((genes.primary), rownames(primary.ST)),]
  genes.met <- intersect(rownames(met.mRNA), rownames(met.ST))
  met.mRNA <- met.mRNA[match(genes.met, rownames(met.mRNA)),]
  met.ST <- met.ST[match((genes.met), rownames(met.ST)),]
  # Compute concordant DGE
  primary.DGE <- DGE(primary.ST, primary.mRNA)
  met.DGE <- DGE(met.ST, met.mRNA)
  primary.DGE.concordant <- concordant_DGE(primary.DGE)
  met.DGE.concordant <- concordant_DGE(met.DGE)
  # Create table holding concordant DGE results
  concordant.DGE <- cbind(primary.DGE.concordant, met.DGE.concordant[,c(3,4)])
  colnames(concordant.DGE) <- c("Gene", "Alteration",
                                "FC Expression Primary", "Primary DGE FDR",
                                "FC Expression Metastatic", "Metastatic DGE FDR")
  #  Update summary table with concordant DGE results
  summary.molecular <- cbind(summary.enrichment, concordant.DGE[,3:6])
  colnames(summary.molecular) <- c("Hugo Symbol", "Alteration",
                                   "Primary Alteration Ratio Difference", "Primary Alteration FDR",
                                   "Metastatic Alteration Ratio Difference", "Metastatic Alteration FDR",
                                   "FC Expression Primary", "Primary DGE FDR",
                                   "FC Expression Metastatic", "Metastatic DGE FDR")
  summary.molecular <- na.omit(summary.molecular)
  # GEX filtering
  summary.molecular <- summary.molecular[((summary.molecular$`Primary DGE FDR` <= DGE.FDR.cutoff) &
                                         (summary.molecular$`Metastatic DGE FDR` <= DGE.FDR.cutoff)),]
  #----Clinical Analysis----
  # Harmonize alteration profiles, mRNA data, and gleason data for gleason grade group analysis.
  # If the user chooses "LOF" as the alteration, then Taylor.primary.LOF, TCGA.LOF, and Barbieri.LOF, are the input datasets.
  # If the user chooses "GOF" as the alteration, then Taylor.primary.GOF, TCGA.GOF, and Barbieri.GOF, are the input datasets.
  # Find overlapping genes between alteration data and mRNA profiles
  # Taylor primary
  Taylor.primary.ST <- alteration_subset(Taylor.primary.LOF, c(query.gene, "altered"))
  Taylor.primary.patients <- intersect(colnames(Taylor.primary.ST), colnames(Taylor.mRNA))
  # Barbieri
  Barbieri.ST <- alteration_subset(Barbieri.LOF, c(query.gene, "altered"))
  Barbieri.patients <- intersect(colnames(Barbieri.ST), colnames(Barbieri.mRNA))
  # TCGA
  TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
  TCGA.patients <- intersect(colnames(TCGA.ST), colnames(TCGA.mRNA))
  # Make Inputs for Gleason Analysis
  # All datasets have patients with matched alteration/mRNA
  if ((length(Taylor.primary.patients) > 0) & ((length(Barbieri.patients)) > 0) &
      (length(TCGA.patients) > 0)){
    Taylor.ST <- alteration_subset(Taylor.primary.LOF, c(query.gene, "altered"))
    Taylor.WT <- alteration_subset(Taylor.primary.alteration, c(query.gene, "unaltered"))
    Barbieri.ST <- alteration_subset(Barbieri.LOF, c(query.gene, "altered"))
    Barbieri.WT <- alteration_subset(Barbieri.alteration, c(query.gene, "unaltered"))
    TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
    TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
    Taylor.mRNA <- Taylor.mae@ExperimentList@listData$gex.rma
    Taylor.primary.mRNA <- scale(Taylor.mRNA[,match(colnames(Taylor.primary.ST), colnames(Taylor.mRNA))])
    Barbieri.mRNA <- Barbieri.mae@ExperimentList@listData$gex.relz
    TCGA.mRNA <- scale(TCGA.mae@ExperimentList@listData$gex.rsem.log)
    genes <- intersect(rownames(Taylor.ST), rownames(Barbieri.ST))
    genes <- intersect(genes, rownames(TCGA.ST))
    genes <- intersect(genes, rownames(Taylor.primary.mRNA))
    genes <- intersect(genes, rownames(Barbieri.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, summary.molecular$`Hugo Symbol`)
    # Reformat Taylor data for Gleason Analysis
    Taylor.ST <- Taylor.ST[match(genes, rownames(Taylor.ST)),]
    Taylor.WT <- Taylor.WT[match(genes, rownames(Taylor.WT)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    Taylor.gleason <- data.frame(Taylor.mae$sample_name, Taylor.mae$grade_group)
    colnames(Taylor.gleason) <- c("Sample ID", "Grade Group")
    # Reformat Barbieri data for Gleason Analysis
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.WT <- Barbieri.WT[match(genes, rownames(Barbieri.WT)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    Barbieri.gleason <- data.frame(Barbieri.mae$sample_name, Barbieri.mae$grade_group)
    colnames(Barbieri.gleason) <- c("Sample ID", "Grade Group")
    # Reformat TCGA data for Gleason Analysis
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.WT <- TCGA.WT[match(genes, rownames(TCGA.WT)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    TCGA.gleason <- data.frame(TCGA.mae$sample_name, TCGA.mae$grade_group)
    colnames(TCGA.gleason) <- c("Sample ID", "Grade Group")
    # Combine alteration profiles
    combined.ST <- cbind(Taylor.ST, Barbieri.ST)
    combined.ST <- cbind(combined.ST, TCGA.ST)
    combined.WT <- cbind(Taylor.WT, Barbieri.WT)
    combined.WT <- cbind(combined.WT, TCGA.WT)
    # Combine mRNA profiles
    combined.mRNA <- cbind(Taylor.primary.mRNA, Barbieri.mRNA)
    combined.mRNA <- cbind(combined.mRNA, TCGA.mRNA)
    # Combine gleason scores
    combined.gleason <- rbind(Taylor.gleason, Barbieri.gleason)
    combined.gleason <- rbind(combined.gleason, TCGA.gleason)
  }
  # No patients from Barbieri data
  if ((length(Taylor.primary.patients) > 0) & ((length(Barbieri.patients)) <= 0) &
      (length(TCGA.patients) > 0)){
    Taylor.ST <- alteration_subset(Taylor.primary.LOF, c(query.gene, "altered"))
    Taylor.WT <- alteration_subset(Taylor.primary.alteration, c(query.gene, "unaltered"))
    TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
    TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
    Taylor.mRNA <- Taylor.mae@ExperimentList@listData$gex.rma
    Taylor.primary.mRNA <- scale(Taylor.mRNA[,match(colnames(Taylor.primary.ST), colnames(Taylor.mRNA))])
    TCGA.mRNA <- scale(TCGA.mae@ExperimentList@listData$gex.rsem.log)
    genes <- intersect(rownames(Taylor.ST), rownames(TCGA.ST))
    genes <- intersect(genes, rownames(Taylor.primary.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, summary.molecular$`Hugo Symbol`)
    # Reformat Taylor data for Gleason Analysis
    Taylor.ST <- Taylor.ST[match(genes, rownames(Taylor.ST)),]
    Taylor.WT <- Taylor.WT[match(genes, rownames(Taylor.WT)),]
    Taylor.primary.mRNA <- Taylor.primary.mRNA[match(genes, rownames(Taylor.primary.mRNA)),]
    Taylor.gleason <- data.frame(Taylor.mae$sample_name, Taylor.mae$grade_group)
    colnames(Taylor.gleason) <- c("Sample ID", "Grade Group")
    # Reformat TCGA data for Gleason Analysis
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.WT <- TCGA.WT[match(genes, rownames(TCGA.WT)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    TCGA.gleason <- data.frame(TCGA.mae$sample_name, TCGA.mae$grade_group)
    colnames(TCGA.gleason) <- c("Sample ID", "Grade Group")
    # Combine alteration profiles
    combined.ST <- cbind(Taylor.ST, TCGA.ST)
    combined.WT <- cbind(Taylor.WT, TCGA.WT)
    # Combine mRNA profiles
    combined.mRNA <- cbind(Taylor.primary.mRNA, TCGA.mRNA)
    # Combine gleason scores
    combined.gleason <- rbind(Taylor.gleason, TCGA.gleason)
  }
  # No patients from Taylor data
  if ((length(Taylor.primary.patients) <= 0) & ((length(Barbieri.patients)) > 0) &
      (length(TCGA.patients) > 0)){
    Barbieri.ST <- alteration_subset(Barbieri.LOF, c(query.gene, "altered"))
    Barbieri.WT <- alteration_subset(Barbieri.alteration, c(query.gene, "unaltered"))
    TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
    TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
    Barbieri.mRNA <- Barbieri.mae@ExperimentList@listData$gex.relz
    TCGA.mRNA <- scale(TCGA.mae@ExperimentList@listData$gex.rsem.log)
    genes <- intersect(rownames(Barbieri.ST), rownames(TCGA.ST))
    genes <- intersect(genes, rownames(Barbieri.mRNA))
    genes <- intersect(genes, rownames(TCGA.mRNA))
    genes <- intersect(genes, summary.molecular$`Hugo Symbol`)
    # Reformat Barbieri data for Gleason Analysis
    Barbieri.ST <- Barbieri.ST[match(genes, rownames(Barbieri.ST)),]
    Barbieri.WT <- Barbieri.WT[match(genes, rownames(Barbieri.WT)),]
    Barbieri.mRNA <- Barbieri.mRNA[match(genes, rownames(Barbieri.mRNA)),]
    Barbieri.gleason <- data.frame(Barbieri.mae$sample_name, Barbieri.mae$grade_group)
    colnames(Barbieri.gleason) <- c("Sample ID", "Grade Group")
    # Reformat TCGA data for Gleason Analysis
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.WT <- TCGA.WT[match(genes, rownames(TCGA.WT)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    TCGA.gleason <- data.frame(TCGA.mae$sample_name, TCGA.mae$grade_group)
    colnames(TCGA.gleason) <- c("Sample ID", "Grade Group")
    # Combine alteration profiles
    combined.ST <- cbind(Barbieri.ST, TCGA.ST)
    combined.WT <- cbind(Barbieri.WT, TCGA.WT)
    # Combine mRNA profiles
    combined.mRNA <- cbind(Barbieri.mRNA, TCGA.mRNA)
    # Combine gleason scores
    combined.gleason <- rbind(Barbieri.gleason, TCGA.gleason)
  }
  # No patients from Barbieri or Taylor data
  if ((length(Taylor.primary.patients) <= 0) & ((length(Barbieri.patients)) <= 0) &
      (length(TCGA.patients) > 0)){
    TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
    TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
    TCGA.mRNA <- scale(TCGA.mae@ExperimentList@listData$gex.rsem.log)
    genes <- intersect(rownames(TCGA.ST), rownames(TCGA.mRNA))
    genes <- intersect(genes, summary.molecular$`Hugo Symbol`)
    # Reformat TCGA data for Gleason Analysis
    TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
    TCGA.WT <- TCGA.WT[match(genes, rownames(TCGA.WT)),]
    TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
    TCGA.gleason <- data.frame(TCGA.mae$sample_name, TCGA.mae$grade_group)
    colnames(TCGA.gleason) <- c("Sample ID", "Grade Group")
    # Combine alteration profiles
    combined.ST <- TCGA.ST
    combined.WT <- TCGA.WT
    # Combine mRNA profiles
    combined.mRNA <- TCGA.mRNA
    # Combine gleason scores
    combined.gleason <- TCGA.gleason
  }
  #  Compute association of gene expression with gleason grade group in ST and WT primary tumors
  source("gleason.R")
  WT.grade <- gleason(combined.WT, combined.mRNA, combined.gleason)
  ST.grade <- gleason(combined.ST, combined.mRNA, combined.gleason)
  # If the user chooses "LOF" as the alteration, then Taylor.primary.LOF and TCGA.LOF are the input datasets.
  # If the user chooses "GOF" as the alteration, then Taylor.primary.GOF and TCGA.GOF are the input datasets.
  # Find overlapping genes between alteration data and mRNA profiles
  TCGA.ST <- alteration_subset(TCGA.LOF, c(query.gene, "altered"))
  TCGA.WT <- alteration_subset(TCGA.alteration, c(query.gene, "unaltered"))
  TCGA.mRNA <- scale(TCGA.mae@ExperimentList@listData$gex.rsem.log)
  genes <- rownames(TCGA.mRNA)
  genes <- intersect(genes, summary.molecular$`Hugo Symbol`)
  # Reformat TCGA
  TCGA.ST <- TCGA.ST[match(genes, rownames(TCGA.ST)),]
  TCGA.WT <- TCGA.WT[match(genes, rownames(TCGA.WT)),]
  TCGA.mRNA <- TCGA.mRNA[match(genes, rownames(TCGA.mRNA)),]
  TCGA.PFS <- data.frame(TCGA.mae$sample_name,
                         TCGA.mae$patient_id,
                         TCGA.mae$days_to_disease_specific_recurrence/30,
                         TCGA.mae$disease_specific_recurrence_status)
  samples <- intersect(colnames(TCGA.mRNA), TCGA.PFS[,1])
  TCGA.mRNA <- TCGA.mRNA[,match(samples, colnames(TCGA.mRNA))]
  TCGA.PFS <- TCGA.PFS[match(samples, TCGA.PFS[,1]),]
  colnames(TCGA.PFS) <- c("Sample ID", "Patient ID", "Time", "Status")
  # Compute association of gene expression with progression free survival in ST and WT primary tumors
  # TCGA PFS pvals
  source("survival.R")
  WT.PFS.TCGA <- survival(TCGA.WT, TCGA.mRNA, TCGA.PFS)
  ST.PFS.TCGA <- survival(TCGA.ST, TCGA.mRNA, TCGA.PFS)
  # Combine WT and ST gleason data
  Grade.Enrichment <- data.frame(WT.grade$`Hugo Symbol`,
                                 WT.grade$`Odds Ratio`, WT.grade$`Tumor Grade FDR`,
                                 ST.grade$`Odds Ratio`, ST.grade$`Tumor Grade FDR`)
  colnames(Grade.Enrichment) <- c("Hugo Symbol",
                                  "WT Odds Ratio", "WT Gleason Enrichment FDR",
                                  "ST Odds Ratio", "ST Gleason Enrichment FDR")
  # Combine WT and ST survival data
  survival <- data.frame(WT.PFS.TCGA$`Hugo Symbol`,
                         WT.PFS.TCGA$`Odds Ratio (High Expressers:Low Expressers)`, WT.PFS.TCGA$`PFS pval`, 
                         ST.PFS.TCGA$`Odds Ratio (High Expressers:Low Expressers)`, ST.PFS.TCGA$`PFS pval`)
  colnames(survival) <- c("Hugo Symbol",
                          "TCGA WT Hazard Ratio", "TCGA WT PFS pval",
                          "TCGA ST Hazard Ratio", "TCGA ST PFS pval")
  #Add Gleason and Survival associations
  summary.full <- summary.molecular
  # Add gleason association
  for (i in 1:nrow(summary.full)){
    index <- match(summary.full$`Hugo Symbol`[i], Grade.Enrichment$`Hugo Symbol`)
    summary.full$`WT Odds Ratio`[i] <- Grade.Enrichment$`WT Odds Ratio`[index]
    summary.full$`WT Gleason Enrichment FDR`[i] <- Grade.Enrichment$`WT Gleason Enrichment FDR`[index]
    summary.full$`ST Odds Ratio`[i] <- Grade.Enrichment$`ST Odds Ratio`[index]
    summary.full$`ST Gleason Enrichment FDR`[i] <- Grade.Enrichment$`ST Gleason Enrichment FDR`[index]
  }
  # Add PFS association
  for (i in 1:nrow(summary.full)){
    index <- match(summary.full$`Hugo Symbol`[i], survival$`Hugo Symbol`)
    summary.full$`TCGA WT Hazard Ratio`[i] <- survival$`TCGA WT Hazard Ratio`[index]
    summary.full$`TCGA WT PFS pval`[i] <- survival$`TCGA WT PFS pval`[index]
    summary.full$`TCGA ST Hazard Ratio`[i] <- survival$`TCGA ST Hazard Ratio`[index]
    summary.full$`TCGA ST PFS pval`[i] <- survival$`TCGA ST PFS pval`[index]
  }
  # Remove hits lacking Gleason and Survival Association
  summary.full <- summary.full[is.na(summary.full$`ST Odds Ratio`) == FALSE &
                                 is.na(summary.full$`TCGA ST Hazard Ratio`) == FALSE,]
  # Capture genes with concordant ST Gleason odds ratio and ST TCGA PFS concordance
  index <- ((summary.full$`ST Odds Ratio` > 1 & summary.full$`TCGA ST Hazard Ratio` > 1) |
              (summary.full$`ST Odds Ratio` < 1 & summary.full$`TCGA ST Hazard Ratio` < 1))
  summary.full <- summary.full[index,]
  # Capture TSG / Oncogene concordance between metastasis, ST Gleason odds ratio, and ST TCGA PFS
  index1 <- (summary.full$Alteration == "LOF" & summary.full$`Metastatic Alteration Ratio Difference` < 0 & 
               summary.full$`ST Odds Ratio` > 1 & summary.full$`TCGA ST Hazard Ratio` > 1)
  index2 <- (summary.full$Alteration == "LOF" & summary.full$`Metastatic Alteration Ratio Difference` > 0 & 
               summary.full$`ST Odds Ratio` < 1 & summary.full$`TCGA ST Hazard Ratio` < 1)
  index3 <- (summary.full$Alteration == "GOF" & summary.full$`Metastatic Alteration Ratio Difference` < 0 & 
               summary.full$`ST Odds Ratio` < 1 & summary.full$`TCGA ST Hazard Ratio` < 1)
  index4 <- (summary.full$Alteration == "GOF" & summary.full$`Metastatic Alteration Ratio Difference` > 0 &
               summary.full$`ST Odds Ratio` > 1 & summary.full$`TCGA ST Hazard Ratio` > 1)
  concordancy <- index1 | index2 |index3 | index4
  summary.full <- summary.full[concordancy,]
  # Filter Survival
  summary.full <- summary.full[((summary.full$`TCGA ST PFS pval` <= ST.cutoff) &
                                  (summary.full$`TCGA WT PFS pval` >= WT.cutoff)),]
  #----Finalize and Filter Hit Table----
  # Add locus
  loci <- annotate(summary.full[,c(1,2)])
  summary.full$`Cytogenetic Band` <- loci$Band
  # Finalize hit table
  hits <- data.frame(summary.full$`Hugo Symbol`, summary.full$`Cytogenetic Band`, summary.full$Alteration,
                     summary.full$`Primary Alteration Ratio Difference`,
                     summary.full$`Primary Alteration FDR`,
                     summary.full$`Metastatic Alteration Ratio Difference`,
                     summary.full$`Metastatic Alteration FDR`,
                     summary.full$`FC Expression Primary`,
                     summary.full$`Primary DGE FDR`,
                     summary.full$`FC Expression Metastatic`,
                     summary.full$`Metastatic DGE FDR`,
                     summary.full$`WT Odds Ratio`,
                     summary.full$`WT Gleason Enrichment FDR`,
                     summary.full$`ST Odds Ratio`,
                     summary.full$`ST Gleason Enrichment FDR`,
                     summary.full$`TCGA WT Hazard Ratio`,
                     summary.full$`TCGA WT PFS pval`,
                     summary.full$`TCGA ST Hazard Ratio`,
                     summary.full$`TCGA ST PFS pval`)
  colnames(hits) <- c("Gene", "Band", "Alteration",
                      "Primary Alteration Ratio Difference", "Primary Enrichment FDR",
                      "Metastatic Alteration Ratio Difference", "Metastatic Enrichment FDR",
                      "FC Expression Primary", "Primary DGE FDR",
                      "FC Expression Metastatic", "Metastatic DGE FDR",
                      "WT Odds Ratio", "WT Gleason Enrichment FDR",
                      "ST Odds Ratio", "ST Gleason Enrichment FDR",
                      "TCGA WT Hazard Ratio", "TCGA WT PFS pval",
                      "TCGA ST Hazard Ratio", "TCGA ST PFS pval")
  hits <- prostamine_prioritization(hits)
  hits <- hits[c("Gene", "ProstaMine Score", "Band", "Alteration",
                 "Primary Alteration Ratio Difference", "Primary Enrichment FDR",
                 "Metastatic Alteration Ratio Difference", "Metastatic Enrichment FDR",
                 "FC Expression Primary", "Primary DGE FDR",
                 "FC Expression Metastatic", "Metastatic DGE FDR",
                 "WT Odds Ratio", "WT Gleason Enrichment FDR",
                 "ST Odds Ratio", "ST Gleason Enrichment FDR",
                 "TCGA WT Hazard Ratio", "TCGA WT PFS pval",
                 "TCGA ST Hazard Ratio", "TCGA ST PFS pval")]
  # Output .RData file
  ws.vars <- c("query.gene", "query.gene.alteration", "query.gene.wt", "hits",
               "Coalteration.enrichment", "LOF.coalteration.enrichment", "LOF.met.enrichment", 
               "GOF.coalteration.enrichment", "GOF.met.enrichment",
               "primary.ST", "met.ST", "primary.mRNA", "met.mRNA",
               "combined.ST", "combined.WT", "combined.mRNA", "combined.gleason", "st", "wt",
               "TCGA.ST", "TCGA.WT", "TCGA.mRNA", "TCGA.PFS")
  save(list=ws.vars, file = output.filename)
}