library(PhyloMarkovChains)
library(phangorn)
n <- 4
all.splits <- sapply(allSplits(n), paste, collapse = "")
J <- length(all.splits)

all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)

## Running the Markov chain
M <- 50000 ## Number of iterations
# set.seed(666)
X <- vector(M + 1, mode = "list")
independent <- FALSE

if (independent) { ## IID sampling 
  X[[1]] <- rcoal(n) #rtree(n, rooted = TRUE, equiprob = TRUE)
  for (i in 2:(M + 1)) {
    X[[i]] <- rcoal(n) # rtree(n, rooted = TRUE, equiprob = TRUE)
  }
} else{ ## Proper Markov 
  
  compute <- FALSE
  
  if (compute) {
    comp.time <- system.time(SPR.mat <-
                               as.matrix(rspr_matrix(all.trees) > 0) + 0)
    save(SPR.mat,
         file = paste0("saved_data/SPR_matrix_n=", n, ".RData"))
    write.csv(x = SPR.mat,
              file = paste0("saved_data/SPR_matrix_n=", n, ".csv"))
  } else{
    load(paste0("saved_data/SPR_matrix_n=", n, ".RData"))
  }
  ####
  # Neighbourhood stuff
  incidence.mat <- SPR.mat
  neighbourhood.sizes <- colSums(incidence.mat)
  
  neigh.ratios <- matrix(NA, nrow = K, ncol = K)
  diag(neigh.ratios) <- 0
  for (i in 1:K) {
    neigh.ratios[i, ] <-
      sapply(neighbourhood.sizes[i] / neighbourhood.sizes,
             function(x)
               min(1, x)) /
      neighbourhood.sizes[i]
  }
  
  trans.mat.MH <- neigh.ratios * incidence.mat
  diag(trans.mat.MH) <- 1 - rowSums(trans.mat.MH)
  diag(trans.mat.MH)
  
  States <- paste0("tree_", 1:K)
  
  phylo.MH <- new(
    "markovchain",
    states = States,
    transitionMatrix = trans.mat.MH,
    name = "rootedMH"
  )
  
  ini.state <- sample(States, 1)
  
  TreeInds <- markovchain::markovchainSequence(n = M,
                                               markovchain = phylo.MH,
                                               t0 = ini.state)
  X <- all.trees[match(TreeInds, States)]
}
class(X) <- "multiPhylo"

## Processing output
Splits <- parallel::mclapply(X, get_clades, mc.cores = 6)

Indexes <- parallel::mclapply(Splits,
                              function(x)
                                charmatch(x, all.splits),
                              mc.cores = 6)

clade.matrix <- do.call(rbind,
                        parallel::mclapply(Indexes,
                                           clade_indicators, L = J,
                                           mc.cores = 6))
colnames(clade.matrix) <- all.splits


plot(clade.matrix[, 6], type = "l")
colMeans(clade.matrix)
