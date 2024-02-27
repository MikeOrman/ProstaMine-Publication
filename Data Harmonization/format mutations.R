# Create dataframe for referencing mutations by sample. 
# 1 = predicted to be a LOF alteration, 0 = not predicted to be a LOF alteration
#-------------------------------ABIDA-------------------------------------------
x <- getPCa("abida")
x <- x[["mut"]]
Hugo_Symbol <- x@assays@unlistData@ranges@NAMES
Center <- "Abida et al"
Tumor_Sample_Barcode <- x@assays@unlistData@elementMetadata@listData$Tumor_sample_barcode
SIFT_pred <- x@assays@unlistData@elementMetadata@listData[["SIFT_pred"]]
Polyphen2_HVAR_pred <- x@assays@unlistData@elementMetadata@listData[["Polyphen2_HDIV_pred"]]
Abida.mutations <- data.frame(Hugo_Symbol, Center, Tumor_Sample_Barcode, SIFT_pred, Polyphen2_HVAR_pred)
Abida <- data.frame()
for (i in 1:nrow(Abida.mutations)){
    Abida[i,1] <- Abida.mutations$Hugo_Symbol[i]
    Abida[i,2] <- Abida.mutations$Center[i]
    Abida[i,3] <- Abida.mutations$Tumor_Sample_Barcode[i]
    sift_pred <- Abida.mutations$SIFT_pred[i]
    polyphen_pred <- Abida.mutations$Polyphen2_HVAR_pred[i]
    if (sift_pred == "T") {Abida[i,4] <- 0}
    if (sift_pred == "D") {Abida[i,4] <- 1}
    if (sift_pred == ".") {Abida[i,4] <- 0}
    if (polyphen_pred == "B") {Abida[i,5] <- 0}
    if (polyphen_pred == "D") {Abida[i,5] <- 1}
    if (polyphen_pred == "P") {Abida[i,5] <- 1}
    if (polyphen_pred == ".") {Abida[i,5] <- 0}
    if (sum((Abida[i,4] + Abida[i,5])) > 0) {Abida[i,6] <- 1}
    if (sum((Abida[i,4] + Abida[i,5])) == 0) {Abida[i,6] <- 0}
  }
colnames(Abida) <- c("Hugo Symbol", "Study", "Sample", "SIFT LOF Pred", "PolyPhen LOF Pred", "LOF Pred")
#-------------------------------TCGA----------------------------------------
x <- getPCa("tcga")
x <- x[["mut"]]
Hugo_Symbol <- x@assays@unlistData@ranges@NAMES
Center <- "TCGA"
Tumor_Sample_Barcode <- x@assays@unlistData@elementMetadata$Tumor_sample_barcode
Tumor_Sample_Barcode <- gsub("-", ".", Tumor_Sample_Barcode)
SIFT <- x@assays@unlistData@elementMetadata@listData$SIFT
PolyPhen <- x@assays@unlistData@elementMetadata@listData$PolyPhen
TCGA.mutations <- data.frame(Hugo_Symbol, Center, Tumor_Sample_Barcode, SIFT, PolyPhen)
TCGA.mutations[TCGA.mutations == "."] <- NA
TCGA <- data.frame()
for (i in 1:nrow(TCGA.mutations)){
  TCGA[i,1] <- TCGA.mutations$Hugo_Symbol[i]
  TCGA[i,2] <- TCGA.mutations$Center[i]
  TCGA[i,3] <- TCGA.mutations$Tumor_Sample_Barcode[i]
  sift_pred <- strsplit(TCGA.mutations$SIFT[i], split = "(", fixed = TRUE)[[1]][1]
  polyphen_pred <- strsplit(TCGA.mutations$PolyPhen[i], split = "(", fixed = TRUE)[[1]][1]
  if (is.na(sift_pred)) {TCGA[i,4] <- 0}
  else{
    if (sift_pred == "deleterious") {TCGA[i,4] <- 1}
    if (sift_pred == "deleterious_low_confidence") {TCGA[i,4] <- 0}
    if (sift_pred == "tolerated") {TCGA[i,4] <- 0}
    if (sift_pred == "tolerated_low_confidence") {TCGA[i,4] <- 0}
  }
  if (is.na(polyphen_pred)) {TCGA[i,5] <- 0}
  else{
    if (polyphen_pred == "benign") {TCGA[i,5] <- 0}
    if (polyphen_pred == "possibly_damaging") {TCGA[i,5] <- 0}
    if (polyphen_pred == "probably_damaging") {TCGA[i,5] <- 1}
    if (polyphen_pred == "unknown") {TCGA[i,5] <- 0}
  }
  if (sum((TCGA[i,4] + TCGA[i,5])) > 0) {TCGA[i,6] <- 1}
  if (sum((TCGA[i,4] + TCGA[i,5])) == 0) {TCGA[i,6] <- 0}
}
colnames(TCGA) <- c("Hugo Symbol", "Study", "Sample", "SIFT LOF Pred", "PolyPhen LOF Pred", "LOF Pred")
#-------------------------------BACA--------------------------------------------
x <- getPCa("baca")
x <- x[["mut"]]
Hugo_Symbol <- x@assays@unlistData@ranges@NAMES
Center <- "Baca et al"
Tumor_Sample_Barcode <- x@assays@unlistData@elementMetadata$Tumor_sample_barcode
SIFT <- x@assays@unlistData@elementMetadata@listData$SIFT
PolyPhen <- x@assays@unlistData@elementMetadata@listData$PolyPhen
Baca.mutations <- data.frame(Hugo_Symbol, Center, Tumor_Sample_Barcode, SIFT, PolyPhen)
Baca.mutations[Baca.mutations == "."] <- NA
Baca <- data.frame()
for (i in 1:nrow(Baca.mutations)){
  Baca[i,1] <- Baca.mutations$Hugo_Symbol[i]
  Baca[i,2] <- Baca.mutations$Center[i]
  Baca[i,3] <- Baca.mutations$Tumor_Sample_Barcode[i]
  sift_pred <- strsplit(Baca.mutations$SIFT[i], split = "(", fixed = TRUE)[[1]][1]
  polyphen_pred <- strsplit(Baca.mutations$PolyPhen[i], split = "(", fixed = TRUE)[[1]][1]
  if (is.na(sift_pred)) {Baca[i,4] <- 0}
  else{
    if (sift_pred == "deleterious") {Baca[i,4] <- 1}
    if (sift_pred == "deleterious_low_confidence") {Baca[i,4] <- 0}
    if (sift_pred == "tolerated") {Baca[i,4] <- 0}
    if (sift_pred == "tolerated_low_confidence") {Baca[i,4] <- 0}
  }
  if (is.na(polyphen_pred)) {Baca[i,5] <- 0}
  else{
    if (polyphen_pred == "benign") {Baca[i,5] <- 0}
    if (polyphen_pred == "possibly_damaging") {Baca[i,5] <- 0}
    if (polyphen_pred == "probably_damaging") {Baca[i,5] <- 1}
    if (polyphen_pred == "unknown") {Baca[i,5] <- 0}
  }
  if (sum((Baca[i,4] + Baca[i,5])) > 0) {Baca[i,6] <- 1}
  if (sum((Baca[i,4] + Baca[i,5])) == 0) {Baca[i,6] <- 0}
}
colnames(Baca) <- c("Hugo Symbol", "Study", "Sample", "SIFT LOF Pred", "PolyPhen LOF Pred", "LOF Pred")
#-------------------------------TAYLOR--------------------------------------------
x <- getPCa("taylor")
x <- x[["mut"]]
Hugo_Symbol <- x@assays@unlistData@ranges@NAMES
Center <- "Taylor et al"
Tumor_Sample_Barcode <- x@assays@unlistData@elementMetadata$Tumor_sample_barcode
SIFT <- x@assays@unlistData@elementMetadata@listData$SIFT
PolyPhen <- x@assays@unlistData@elementMetadata@listData$PolyPhen
Taylor.mutations <- data.frame(Hugo_Symbol, Center, Tumor_Sample_Barcode, SIFT, PolyPhen)
Taylor.mutations[Taylor.mutations == "."] <- NA
Taylor <- data.frame()
for (i in 1:nrow(Taylor.mutations)){
  Taylor[i,1] <- Taylor.mutations$Hugo_Symbol[i]
  Taylor[i,2] <- Taylor.mutations$Center[i]
  Taylor[i,3] <- Taylor.mutations$Tumor_Sample_Barcode[i]
  sift_pred <- strsplit(Taylor.mutations$SIFT[i], split = "(", fixed = TRUE)[[1]][1]
  polyphen_pred <- strsplit(Taylor.mutations$PolyPhen[i], split = "(", fixed = TRUE)[[1]][1]
  if (is.na(sift_pred)) {Taylor[i,4] <- 0}
  else{
    if (sift_pred == "deleterious") {Taylor[i,4] <- 1}
    if (sift_pred == "deleterious_low_confidence") {Taylor[i,4] <- 0}
    if (sift_pred == "tolerated") {Taylor[i,4] <- 0}
    if (sift_pred == "tolerated_low_confidence") {Taylor[i,4] <- 0}
  }
  if (is.na(polyphen_pred)) {Taylor[i,5] <- 0}
  else{
    if (polyphen_pred == "benign") {Taylor[i,5] <- 0}
    if (polyphen_pred == "possibly_damaging") {Taylor[i,5] <- 0}
    if (polyphen_pred == "probably_damaging") {Taylor[i,5] <- 1}
    if (polyphen_pred == "unknown") {Taylor[i,5] <- 0}
  }
  if (sum((Taylor[i,4] + Taylor[i,5])) > 0) {Taylor[i,6] <- 1}
  if (sum((Taylor[i,4] + Taylor[i,5])) == 0) {Taylor[i,6] <- 0}
}
colnames(Taylor) <- c("Hugo Symbol", "Study", "Sample", "SIFT LOF Pred", "PolyPhen LOF Pred", "LOF Pred")
#-------------------------------BARBIERI--------------------------------------------
x <- getPCa("barbieri")
x <- x[["mut"]]
Hugo_Symbol <- x@assays@unlistData@ranges@NAMES
Center <- "Barbieri et al"
Tumor_Sample_Barcode <- x@assays@unlistData@elementMetadata$Tumor_sample_barcode
SIFT <- x@assays@unlistData@elementMetadata@listData$SIFT
PolyPhen <- x@assays@unlistData@elementMetadata@listData$PolyPhen
Barbieri.mutations <- data.frame(Hugo_Symbol, Center, Tumor_Sample_Barcode, SIFT, PolyPhen)
Barbieri.mutations[Barbieri.mutations == "."] <- NA
Barbieri <- data.frame()
for (i in 1:nrow(Barbieri.mutations)){
  Barbieri[i,1] <- Barbieri.mutations$Hugo_Symbol[i]
  Barbieri[i,2] <- Barbieri.mutations$Center[i]
  Barbieri[i,3] <- Barbieri.mutations$Tumor_Sample_Barcode[i]
  sift_pred <- strsplit(Barbieri.mutations$SIFT[i], split = "(", fixed = TRUE)[[1]][1]
  polyphen_pred <- strsplit(Barbieri.mutations$PolyPhen[i], split = "(", fixed = TRUE)[[1]][1]
  if (is.na(sift_pred)) {Barbieri[i,4] <- 0}
  else{
    if (sift_pred == "deleterious") {Barbieri[i,4] <- 1}
    if (sift_pred == "deleterious_low_confidence") {Barbieri[i,4] <- 0}
    if (sift_pred == "tolerated") {Barbieri[i,4] <- 0}
    if (sift_pred == "tolerated_low_confidence") {Barbieri[i,4] <- 0}
  }
  if (is.na(polyphen_pred)) {Barbieri[i,5] <- 0}
  else{
    if (polyphen_pred == "benign") {Barbieri[i,5] <- 0}
    if (polyphen_pred == "possibly_damaging") {Barbieri[i,5] <- 0}
    if (polyphen_pred == "probably_damaging") {Barbieri[i,5] <- 1}
    if (polyphen_pred == "unknown") {Barbieri[i,5] <- 0}
  }
  if (sum((Barbieri[i,4] + Barbieri[i,5])) > 0) {Barbieri[i,6] <- 1}
  if (sum((Barbieri[i,4] + Barbieri[i,5])) == 0) {Barbieri[i,6] <- 0}
}
colnames(Barbieri) <- c("Hugo Symbol", "Study", "Sample", "SIFT LOF Pred", "PolyPhen LOF Pred", "LOF Pred")
#-------------------------------COMBINED----------------------------------------
all_mutations <- rbind(Abida, TCGA)
all_mutations <- rbind(all_mutations, Baca)
all_mutations <- rbind(all_mutations, Taylor)
all_mutations <- rbind(all_mutations, Barbieri)
write.table(all_mutations, file = "Mutations.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Taylor, file = "Taylor mutations.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(Barbieri, file = "Barbieri mutations.txt", sep = "\t", col.names = TRUE, quote = FALSE)
write.table(TCGA, file = "TCGA mutations.txt", sep = "\t", col.names = TRUE, quote = FALSE)