#--------------------------LOAD cBIO CNA DATA-----------------------------------
# TCGA
TCGA.mae <- getPCa("tcga")
TCGA <- TCGA.mae[["cna.gistic"]]
TCGA <- as.data.frame(TCGA)
TCGA.sample <- data.frame(TCGA.mae@colData@listData)
nkx.altered.patients <- read.table("NKXaltered.txt")
for (i in 1:nrow(nkx.altered.patients)){
  nkx.altered.patients[i,1] <- substr(nkx.altered.patients[i,1], 30, 44)
  nkx.altered.patients[i,1] <- gsub("-", ".", nkx.altered.patients[i,1], )
}
vals <- c()
for (i in 1:ncol(TCGA)){
  logic <- (match(colnames(TCGA)[i], nkx.altered.patients$V1) >= 1)
  if (is.na(logic) == FALSE){
    vals[i] <- -1
  } else {vals[i] <- 0}
}
df <- t(vals)
colnames(df) <- colnames(TCGA)
rownames(df) <- "NKX3-1"
TCGA <- rbind(df, TCGA)
# Taylor
Taylor.mae <- getPCa("taylor")
Taylor <- Taylor.mae[["cna.gistic"]]
Taylor <- as.data.frame(Taylor)
Taylor.sample <- data.frame(Taylor.mae@colData@listData)
# Baca
Baca.mae <- getPCa("baca")
Baca <- Baca.mae[["cna.gistic"]]
Baca <- as.data.frame(Baca)
Baca.sample <- data.frame(Baca.mae@colData@listData)
# Barbieri
Barbieri.mae <- getPCa("barbieri")
Barbieri <- Barbieri.mae[["cna.gistic"]]
Barbieri <- as.data.frame(Barbieri)
Barbieri.sample <- data.frame(Barbieri.mae@colData@listData)
# Hieronymous
Hieronymus.mae <- getPCa("hieronymus")
Hieronymus <- Hieronymus.mae[["cna.gistic"]]
Hieronymus <- as.data.frame(Hieronymus)
Hieronmous.sample <- data.frame(Hieronymus.mae@colData@listData)
# Abida
Abida.mae <- getPCa("abida")
Abida <- Abida.mae[["cna.gistic"]]
Abida <- as.data.frame(Abida)
#--------------------------SEPARATE PRIMARY AND MET SAMPLES FROM TAYLOR AND BACA----------
Taylor.primary.samples <- Taylor.sample$sample_name[Taylor.sample$sample_type == "primary"]
Taylor.met.samples <- Taylor.sample$sample_name[Taylor.sample$sample_type == "metastatic"]
Taylor.primary <- Taylor[,c(1, na.omit(match(Taylor.primary.samples, colnames(Taylor))))]
Taylor.met <- Taylor[,c(1, na.omit(match(Taylor.met.samples, colnames(Taylor))))]
Baca.primary.samples <- Baca.sample$sample_name[Baca.sample$sample_type == "primary"]
Baca.met.samples <- Baca.sample$sample_name[Baca.sample$sample_type == "metastatic"]
Baca.primary <- Baca[,c(1, na.omit(match(Baca.primary.samples, colnames(Baca))))]
Baca.met <- Baca[,c(1, na.omit(match(Baca.met.samples, colnames(Baca))))]
#--------------------------COMBINE DATASETS------------------------
# Combine primary sets
genes <- intersect(row.names(TCGA), row.names(Taylor.primary))
genes <- intersect(row.names(Baca.primary), genes)
genes <- intersect(row.names(Barbieri), genes)
genes <- intersect(row.names(Hieronymus), genes)
TCGA <- TCGA[match(genes, row.names(TCGA)),]
Taylor.primary <- Taylor.primary[match(genes, row.names(Taylor.primary)),]
Baca.primary <- Baca.primary[match(genes, row.names(Baca.primary)),]
Barbieri <- Barbieri[match(genes, row.names(Barbieri)),]
Hieronymus <- Hieronymus[match(genes, row.names(Hieronymus)),]
primary <- cbind(TCGA, Taylor.primary)
primary <- cbind(primary, Baca.primary)
primary <- cbind(primary, Barbieri)
primary <- cbind(primary, Hieronymus)
# Combine met sets
genes <- intersect(row.names(Abida), row.names(Taylor.met))
genes <- intersect(row.names(Baca.met), genes)
Abida <- Abida[match(genes, row.names(Abida)),]
Taylor.met <- Taylor.met[match(genes, row.names(Taylor.met)),]
Baca.met <- Baca.met[match(genes, row.names(Baca.met)),]
met <- cbind(Abida, Taylor.met)
met <- cbind(met, Baca.met)
# Filter genes missing in both primary and metastatic settings
genes <- intersect(row.names(primary), row.names(met))
# Filter primary and metastatic genes. The number of rows for primary and metastatic datasets should now be equal.
primary <- primary[match(genes, row.names(primary)),]
met <- met[match(genes, row.names(met)),]
#--------------------------TRANSFORM CNA DATA TO ALTERATION MATRIX--------
# Primary alteration matrix
primary.altered <- primary
primary.altered[primary > 0] <- 1
primary.altered[primary < 0] <- 1
# Metastatic alteration matrix
met.altered <- met
met.altered[met > 0] <- 1
met.altered[met < 0] <- 1
# TCGA alteration matrix
TCGA.altered <- TCGA
TCGA.altered[TCGA > 0] <- 1
TCGA.altered[TCGA < 0] <- 1
# Taylor alteration matrix
Taylor.primary.altered <- Taylor.primary
Taylor.primary.altered[Taylor.primary > 0] <- 1
Taylor.primary.altered[Taylor.primary < 0] <- 1
Taylor.met.altered <- Taylor.met
Taylor.met.altered[Taylor.met > 0] <- 1
Taylor.met.altered[Taylor.met < 0] <- 1
# Barbieri alteration matrix
Barbieri.altered <- Barbieri
Barbieri.altered[Barbieri > 0] <- 1
Barbieri.altered[Barbieri < 0] <- 1
#--------------------------TRANSFORM CNA DATA TO GOF & LOF ALTTERATION MATRICES-----
# Primary LOF binary matrix
primary.LOF <- primary
primary.LOF[primary < 0] <- 1
primary.LOF[primary > 0] <- 0
# Primary GOF binary matrix
primary.GOF <- primary
primary.GOF[primary > 0] <- 1
primary.GOF[primary < 0] <- 0
# Met LOF binary matrix
met.LOF <- met
met.LOF[met < 0] <- 1
met.LOF[met > 0] <- 0
# Met GOF binary matrix
met.GOF <- met
met.GOF[met > 0] <- 1
met.GOF[met < 0] <- 0
# TCGA LOF binary matrix
TCGA.LOF <- TCGA
TCGA.LOF[TCGA < 0] <- 1
TCGA.LOF[TCGA > 0] <- 0
# TCGA GOF binary matrix
TCGA.GOF <- TCGA
TCGA.GOF[TCGA < 0] <- 0
TCGA.GOF[TCGA > 0] <- 1
# Taylor primary LOF binary matrix
Taylor.primary.LOF <- Taylor.primary
Taylor.primary.LOF[Taylor.primary < 0] <- 1
Taylor.primary.LOF[Taylor.primary > 0] <- 0
# Taylor primary GOF binary matrix
Taylor.primary.GOF <- Taylor.primary
Taylor.primary.GOF[Taylor.primary < 0] <- 0
Taylor.primary.GOF[Taylor.primary > 0] <- 1
# Taylor met LOF binary matrix
Taylor.met.LOF <- Taylor.met
Taylor.met.LOF[Taylor.met < 0] <- 1
Taylor.met.LOF[Taylor.met > 0] <- 0
# Taylor met GOF binary matrix
Taylor.met.GOF <- Taylor.met
Taylor.met.GOF[Taylor.met < 0] <- 0
Taylor.met.GOF[Taylor.met > 0] <- 1
# Barbieri LOF binary matrix
Barbieri.LOF <- Barbieri
Barbieri.LOF[Barbieri < 0] <- 1
Barbieri.LOF[Barbieri > 0] <- 0
# Barbieri GOF binary matrix
Barbieri.GOF <- Barbieri
Barbieri.GOF[Barbieri < 0] <- 0
Barbieri.GOF[Barbieri > 0] <- 1
#--------------------------ADD DELETERIOUS MUTATIONS TO ALTERATION MATRIX-------
mutations <- read.table("Mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add mutations to primary alteration matrix
# i loop: iterate through primary samples
for (i in 1:ncol(primary.altered)){
  sample <- colnames(primary.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(primary.altered))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        primary.altered[gene.index, i] <- primary.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
# Add mutations to met alteration matrix
# i loop: iterate through met samples
for (i in 1:ncol(met.altered)){
  sample <- colnames(met.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(met.altered))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        met.altered[gene.index, i] <- met.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO LOF ALTERATION MATRIX----
for (i in 1:ncol(primary.LOF)){
  sample <- colnames(primary.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(primary.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        primary.LOF[gene.index, i] <- primary.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
for (i in 1:ncol(met.LOF)){
  sample <- colnames(met.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(met.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        met.LOF[gene.index, i] <- met.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO TCGA ALTERATION MATRIX----
TCGA.mutation <- read.table("TCGA mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add TCGA mutations to TCGA alteration matrix
for (i in 1:ncol(TCGA.altered)){
  sample <- colnames(TCGA.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(TCGA.altered))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        TCGA.altered[gene.index, i] <- TCGA.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO TCGA LOF ALTERATION MATRIX----
for (i in 1:ncol(TCGA.LOF)){
  sample <- colnames(TCGA.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(TCGA.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        TCGA.LOF[gene.index, i] <- TCGA.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO Taylor LOF primary/secondary ALTERATION MATRIX----
Taylor.mutation <- read.table("Taylor mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add Taylor mutations to Taylor alteration matrix
for (i in 1:ncol(Taylor.primary.altered)){
  sample <- colnames(Taylor.primary.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(Taylor.primary.altered))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        Taylor.primary.altered[gene.index, i] <- Taylor.primary.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
# Primary LOF
for (i in 1:ncol(Taylor.primary.LOF)){
  sample <- colnames(Taylor.primary.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(Taylor.primary.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        Taylor.primary.LOF[gene.index, i] <- Taylor.primary.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
# Met LOF
for (i in 1:ncol(Taylor.met.LOF)){
  sample <- colnames(Taylor.met.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(Taylor.met.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        Taylor.met.LOF[gene.index, i] <- Taylor.met.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------ADD DELETERIOUS MUTATIONS TO Barbieri ALTERATION MATRIX----
Barbieri.mutation <- read.table("Barbieri mutations.txt", sep = "\t", header = TRUE, check.names = FALSE)
# Add Barbieri mutations to Barbieri alteration matrix
for (i in 1:ncol(Barbieri.altered)){
  sample <- colnames(Barbieri.altered)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(Barbieri.altered))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        Barbieri.altered[gene.index, i] <- Barbieri.altered[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}

#--------------------------ADD DELETERIOUS MUTATIONS TO Barbieri LOF ALTERATION MATRIX----
for (i in 1:ncol(Barbieri.LOF)){
  sample <- colnames(Barbieri.LOF)[i]
  # Find if sample name exists in mutation data
  sample.index <- match(sample, mutations$Sample)
  # If sample exists in mutation data then initiate k loop
  if (is.na(sample.index) == FALSE) {
    sample.mutations <- mutations[sample.index,]
    # k loop: iterate through the sample's mutations.
    for (k in 1:nrow(sample.mutations)){
      gene <- sample.mutations$`Hugo Symbol`[k]
      gene.index <- match(gene, rownames(Barbieri.LOF))
      # If the gene has a mutation present, then add its LOF prediction to the alteration count.
      if (is.na(gene.index) == FALSE){
        Barbieri.LOF[gene.index, i] <- Barbieri.LOF[gene.index, i] + sample.mutations$`LOF Pred`[k]
      }
    }
  }
  # If gene does not exist in mutation data then move to next sample
}
#--------------------------REVERT MATRICES TO BINARY ALTERATION------------
primary.altered[primary.altered == 2] <- 1
primary.LOF[primary.LOF == 2] <- 1
met.altered[met.altered == 2] <- 1
met.LOF[met.LOF == 2] <- 1
TCGA.altered[TCGA.altered == 2] <- 1
TCGA.LOF[TCGA.LOF == 2] <- 1
Taylor.primary.altered[Taylor.primary.altered == 2] <- 1
Taylor.primary.LOF[Taylor.primary.LOF == 2] <- 1
Taylor.met.LOF[Taylor.met.LOF == 2] <- 1
Barbieri.altered[Barbieri.altered == 2] <- 1
Barbieri.LOF[Barbieri.LOF == 2] <- 1
#--------------------------PATIENT FILTER---------------------------------------
# z <- read.table("Group 1 patients.txt")
# patients <- z$V1
# primary.altered <- primary.altered[,na.omit(match(patients, colnames(primary.altered)))]
# primary.LOF <- primary.LOF[,na.omit(match(patients, colnames(primary.LOF)))]
# primary.GOF <- primary.GOF[,na.omit(match(patients, colnames(primary.GOF)))]
# TCGA.altered <- TCGA.altered[,na.omit(match(patients, colnames(TCGA.altered)))]
# TCGA.LOF <- TCGA.LOF[,na.omit(match(patients, colnames(TCGA.LOF)))]
# TCGA.GOF <- TCGA.GOF[,na.omit(match(patients, colnames(TCGA.GOF)))]
# Taylor.primary.altered <- Taylor.primary.altered[,na.omit(match(patients, colnames(Taylor.primary.altered)))]
# Taylor.primary.LOF <- Taylor.primary.LOF[,na.omit(match(patients, colnames(Taylor.primary.LOF)))]
# Taylor.primary.GOF <- Taylor.primary.GOF[,na.omit(match(patients, colnames(Taylor.primary.GOF)))]
# Barbieri.altered <- Barbieri.altered[,na.omit(match(patients, colnames(Barbieri.altered)))]
# Barbieri.LOF <- Barbieri.LOF[,na.omit(match(patients, colnames(Barbieri.LOF)))]
# Barbieri.GOF <- Barbieri.GOF[,na.omit(match(patients, colnames(Barbieri.GOF)))]
#--------------------------WRITE FORMATTED GOF / LOF ALTERATION TABLES----------
write.table(primary.altered, "primary alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(primary.LOF, "primary LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(primary.GOF, "primary GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.altered, "met alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.LOF, "met LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(met.GOF, "met GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.altered, "TCGA alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.LOF, "TCGA LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA.GOF, "TCGA GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor.primary.altered, "Taylor primary alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor.primary.LOF, "Taylor primary LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor.primary.GOF, "Taylor primary GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor.met.LOF, "Taylor met LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor.met.GOF, "Taylor met GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Barbieri.altered, "Barbieri alteration.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Barbieri.LOF, "Barbieri LOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Barbieri.GOF, "Barbieri GOF.txt", sep = "\t", col.names = TRUE, quote = FALSE)