library(PhyloMarkovChains)
library(phangorn)
n <- 5
all.splits <- sapply(allSplits(n), paste, collapse = "")
J <- length(all.splits)
## Running the Markov chain
M <- 50000
# set.seed(666)
X <- vector(M + 1, mode = "list")
independent <- TRUE
if(independent){
  X[[1]] <- rcoal(n) #rtree(n, rooted = TRUE, equiprob = TRUE)
  for(i in 2:(M+1)){
      X[[i]] <- rcoal(n) # rtree(n, rooted = TRUE, equiprob = TRUE)
  }
}else{
  X[[1]] <- rtree(n, rooted = TRUE, equiprob = TRUE)
  for(i in 2:(M+1)){
    cand <- rSPR(X[[i-1]])
    while(!is.rooted(cand)){
      cand <- rSPR(X[[i-1]])
    }
    X[[i]] <- cand
  }
}
class(X) <- "multiPhylo"
## Processing output
Splits <- parallel::mclapply(X, get_clades, mc.cores = 6)

Indexes <- parallel::mclapply(Splits,
                              function(x) charmatch(x, all.splits),
                              mc.cores = 6)

clade.matrix <- do.call(rbind,
                        parallel::mclapply(Indexes,
                                           clade_indicators, L = J,
                                           mc.cores = 6))
colnames(clade.matrix) <- all.splits


plot(clade.matrix[, 6], type = "l")
colMeans(clade.matrix)
