n <- 4
all.trees <- phangorn::allTrees(n, rooted = TRUE)
K <- length(all.trees)

load(paste0("saved_data/SPR_matrix_n=", n, ".RData"))

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

save(trans.mat.MH,
     file = paste0("saved_data/MH_matrix_n=", n, ".RData"))