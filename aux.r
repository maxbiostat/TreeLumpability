load_clade_matrix <- function(mapfile){
  ## Source: https://github.com/maxbiostat/coalescent_validation/blob/main/aux_clade_efficiency.r
  raw <- read.table(mapfile, header = TRUE) 
  inds <- raw[-1, -1]
  #
  freqfile <- gsub("cladematrix", "cladetable", mapfile)
  freqfile <- gsub(".cmap", ".txt", freqfile)
  ctable <- read.table(freqfile, header = TRUE)
  freqs <- subset(ctable, Cred < 1)
  pos <- match(freqs$No., ctable$No.)
  #
  final.inds <- inds[, pos]
  colnames(final.inds)  <- freqs$Members
  return(list(
    clade_frequencies = freqs,
    clade_indicators = final.inds
  ))
}