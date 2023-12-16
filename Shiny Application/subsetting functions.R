# Subset CNA dataset based on copy number loss, WT/gain genotypes
# dataset = CNA matrix with row names as genes and col names as samples
# x = c("Hugo_Symbol", "type of loss or gain")
subtype_subset <- function (dataset, x) {
  if (x[2] == "loss"){
    samples = t(subset(dataset, rownames(dataset) == x[1]))
    index = which(samples < 0)
    subset = dataset[, c(1,index)]
    return(subset)
  }
  if (x[2]=="gain"){
    samples = t(subset(dataset, rownames(dataset) == x[1]))
    index = which(samples > 0)
    subset = dataset[, c(1,index)]
    return(subset)
  }
  if (x[2] == "WT"){
    samples = t(subset(dataset, rownames(dataset) == x[1]))
    samples <- as.numeric(samples)
    index = which(samples == 0)
    subset = dataset[, c(1,index)]
    return(subset)
  }
}

# Subsets a binary alteration matrix by alteration of a given gene
# Dataset is a binary alteration matrix. Row names are Hugo Symbols and column names are sample IDs.
# Matrix must be clear of NAs
# x = c("Hugo Symbol of alteration background", "unaltered or altered")
alteration_subset <- function (dataset, x) {
  if (x[2] == "altered") {
    index <- match(x[1], rownames(dataset))
    background.alteration.profile <- dataset[index,]
    subset <- dataset[,background.alteration.profile==1]
    return(subset)
  }
  if (x[2] == "unaltered") {
    index <- match(x[1], rownames(dataset))
    background.alteration.profile <- dataset[index,]
    subset <- dataset[,background.alteration.profile==0]
    return(subset)
  }
}


# For oncoprints
primary.loss_func<-function(primary,gene1,gene2){
  subtype_subset(subtype_subset(primary, c(gene1, "loss")), c(gene2, "loss"))
}

primary.wt_func<-function(primary,gene1,gene2){
  subtype_subset(subtype_subset(primary, c(gene1, "WT")), c(gene2, "WT"))
} 

met.loss_func<-function(gene1,gene2){
  subtype_subset(subtype_subset(met, c(gene1, "loss")), c(gene2, "loss"))
}

met.wt_func<-function(gene1,gene2){
  subtype_subset(subtype_subset(met, c(gene1, "WT")), c(gene2, "WT"))
}