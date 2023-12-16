# Dataset = Binary alteration matrix. Row names are Hugo symbols and column names are sample IDs.
alterationRank <- function(dataset){
  gene.names = rownames(dataset)
  alteration.freq = data.frame()
  logical = (dataset > 0)
  for (i in 1:nrow(logical)) {
    alteration.freq[i,1] = gene.names[i]
    alteration.freq[i,2] = sum(logical[i,1:ncol(logical)], na.rm = TRUE)
    alteration.freq[i,3] = (ncol(logical)) - sum(logical[i,1:ncol(logical)], na.rm = TRUE)
    alteration.freq[i,4] = sum(logical[i,1:ncol(logical)], na.rm = TRUE)/(ncol(logical))
  }
  ordered = order(alteration.freq[,4])
  x = c()
  y = c()
  z = c()
  a = c()
  b = c()
  for (i in 1:length(ordered)) {
    x[i] = i
    y[i] = alteration.freq[ordered[i], 4]
    z[i] = gene.names[ordered[i]]
    a[i] = alteration.freq[ordered[i], 2]
    b[i] = alteration.freq[ordered[i], 3]
  }
  output = data.frame(z, x, a, b, y)
  colnames(output) = c("Hugo_Symbol", "Gene Rank", "total altered", 
                       "total unaltered", "alteration freq")
  return(output)
}