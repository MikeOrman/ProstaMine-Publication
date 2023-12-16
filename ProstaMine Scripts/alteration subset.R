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