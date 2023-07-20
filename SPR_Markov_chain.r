library(PhyloMarkovChains)
library(phangorn)
########
n <- 4
all.splits <- sapply(allSplits(n), paste, collapse = "")
J <- length(all.splits)

all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)
States <- paste0("tree_", 1:K)

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
  load(paste0("saved_data/MH_matrix_n=", n, ".RData"))
  TreeInds <- run_one_MC(Niter = M, MH_mat = trans.mat.MH)
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

