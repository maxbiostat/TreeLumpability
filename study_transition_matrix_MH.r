compute_all_errors <- function(clade){
  
  csize <- length(strsplit(clade, "")[[1]])
  clade.ind <- clade_index(clade, all.splits)
  
  S1.c.pos <- which(complete.clade.matrix[, clade.ind] == 1)
  S0.c.pos <- setdiff(1:K, S1.c.pos)
  
  grid.pos <- subset(expand.grid(S1.c.pos, S1.c.pos),
                     Var1 < Var2)
  
  all.errors <- do.call(rbind,
                        apply(grid.pos, 1,
                              function(row){
                                x.pos <- row[1]
                                y.pos <- row[2]
                                sum_x <- sum(trans.mat.MH[x.pos, S0.c.pos])
                                sum_y <- sum(trans.mat.MH[y.pos, S0.c.pos])
                                return(data.frame(
                                  index.x = x.pos,
                                  index.y = y.pos,
                                  Nx = neighbourhood.sizes[x.pos],
                                  Ny = neighbourhood.sizes[y.pos],
                                  A0cx = A0s.sizes[x.pos],
                                  A0cy = A0s.sizes[y.pos],
                                  sum.x = sum_x,
                                  sum.y = sum_y,
                                  error = sum_x - sum_y
                                ))
                              })
  )
  all.errors$ntaxa <- n
  all.errors$clade_size <- csize
  all.errors$clade <- clade
  all.errors <- tibble::tibble(all.errors)
  return(all.errors)
}

compute_all_dists <- function(clade){
  clade.ind <- clade_index(clade, all.splits)
  S1.c.pos <- which(complete.clade.matrix[, clade.ind] == 1)
  grid.pos <- subset(expand.grid(S1.c.pos, S1.c.pos),
                     Var1 < Var2)
  all.dists <- do.call(rbind,
                        apply(grid.pos, 1,
                              function(row){
                                x.pos <- row[1]
                                y.pos <- row[2]
                                return(data.frame(
                                  index.x = x.pos,
                                  index.y = y.pos,
                                  rspr_dist = rspr(
                                    all.trees[[x.pos]],
                                    all.trees[[y.pos]]
                                  )
                                ))
                              })
  )
  all.dists$ntaxa <- n
  all.dists$clade <- clade
  all.dists <- tibble::tibble(all.dists)
  return(all.dists)
}
#########################################
library(PhyloMarkovChains)
library(phangorn)
library(tidyverse)

n <- 4
all.splits <- sapply(allSplits(n + 1),  paste, collapse = "")
J <- length(all.splits)

all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)

compute <- FALSE

if(compute){
  comp.time <- system.time(
    SPR.mat <- as.matrix(rspr_matrix(all.trees) > 0) + 0
  )
  save(SPR.mat,
       file = paste0("saved_data/SPR_matrix_n=", n, ".RData"))
}else{
  load(paste0("saved_data/SPR_matrix_n=", n, ".RData"))
}

####
# Neighbourhood stuff
incidence.mat <- SPR.mat
neighbourhood.sizes <- colSums(incidence.mat)

neigh.ratios <- matrix(NA, nrow = K, ncol = K)
diag(neigh.ratios) <- 0
for(i in 1:K){
  neigh.ratios[i, ] <- sapply(neighbourhood.sizes[i]/neighbourhood.sizes,
                                          function(x) min(1, x))/
    neighbourhood.sizes[i] 
}

trans.mat.MH <- neigh.ratios * incidence.mat
diag(trans.mat.MH) <- 1-rowSums(trans.mat.MH)
diag(trans.mat.MH)

#####

complete.Splits <- parallel::mclapply(all.trees,
                                      get_clades, mc.cores = 12)

complete.indexes <- parallel::mclapply(complete.Splits,
                                       function(x) charmatch(x, all.splits),
                                       mc.cores = 12)

complete.clade.matrix <- do.call(rbind, parallel::mclapply(complete.indexes,
                                                           clade_indicators,
                                                           L = J,
                                                           mc.cores = 12))
colnames(complete.clade.matrix) <- all.splits
colMeans(complete.clade.matrix)

## Computing |A0c|
clade <- "12"
clade.nice <- paste0("t", strsplit(clade, "")[[1]])
clade.ind <- clade_index(x = clade, all.splits = all.splits)
hasclade <- which(complete.clade.matrix[, clade.ind] == 1)
parts <- lapply(all.trees[hasclade],
                safe_f_c, clade = clade.nice)

Nxprimes <- unlist(lapply(parts,
                          function(x) neighbourhood_size(x$x_prime)$n_size))
Nphixcs <- unlist(lapply(parts,
                         function(x) neighbourhood_size(x$phi_x_c)$n_size))

A0s.sizes <- rep(NA, K)
A0s.sizes[hasclade] <- neighbourhood.sizes[hasclade] - Nxprimes + Nphixcs

all.errors <- compute_all_errors(clade)
dists <- compute_all_dists(clade)

with.dists <- tibble::tibble(
  merge(all.errors, dists,
        by = c("index.x", "index.y",
               "ntaxa", "clade") )
)
with.dists <- with.dists %>%
  mutate(bound = A0cx/Nx - (5*A0cy)/(6*Ny))
with.dists

## Plotting
boxplot(error ~ rspr_dist, with.dists)
boxplot(bound ~ rspr_dist, with.dists)

plot(error~bound, with.dists)
abline(a = 0, b = 1, lwd = 2)
plot(abs(error)~bound, with.dists)
abline(a = 0, b = 1, lwd = 2)

range(with.dists$error)
range(with.dists$bound)

with.dists[which.min(with.dists$error), ]
with.dists[which.max(with.dists$error), ]

with.dists[which.min(with.dists$bound), ]
with.dists[which.max(with.dists$bound), ]