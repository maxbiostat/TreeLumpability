get_csize <- function(clade){
  length(strsplit(clade, "")[[1]])
}
get_stuff <- function(clade){
  
  compute_all_errors <- function(clade){
    
    csize <- get_csize(clade)
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
  ###
  
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
  
  all.errors.dt <- tibble::tibble(
    merge(all.errors, dists,
          by = c("index.x", "index.y",
                 "ntaxa", "clade") )
  )
  all.errors.dt <- all.errors.dt %>%
    mutate(bound = A0cx/Nx - (5*A0cy)/(6*Ny))
  return(all.errors.dt)
}
#########################################
library(PhyloMarkovChains)
library(phangorn)
library(tidyverse)

n <- 5
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
clade.sizes <- unlist(lapply(all.splits, get_csize))
nontrivial <- intersect(which(clade.sizes > 1), which(clade.sizes < n))
all.ntclades <- all.splits[nontrivial]
all.errors <- parallel::mclapply(all.ntclades,
                     function(clade) get_stuff(clade), mc.cores = 8)

all.errors.dt <- do.call(rbind, all.errors)

write.csv(all.errors.dt,
          file = paste0("csv/errors_n=", n, ".csv"),
          row.names = FALSE)

## Plotting
boxplot(error ~ rspr_dist, all.errors.dt)
boxplot(bound ~ rspr_dist, all.errors.dt)

plot(error~bound, all.errors.dt)
abline(a = 0, b = 1, lwd = 2)

range(all.errors.dt$error)
range(all.errors.dt$bound)

all.errors.dt[which.min(all.errors.dt$error), ]
all.errors.dt[which.max(all.errors.dt$error), ]

all.errors.dt[which.min(all.errors.dt$bound), ]
all.errors.dt[which.max(all.errors.dt$bound), ]

library(ggplot2)

vsbound <- ggplot(all.errors.dt,
       aes(x = bound, y = error,
           colour = rspr_dist)) +
   geom_point() +
  facet_wrap(clade~., scales = "free") +
  geom_abline(intercept = 0, slope = 1, linetype = "longdash") + 
  theme_bw(base_size = 20)

vsbound

ggsave(
  plot = vsbound,
  filename = paste0("plots/errors_vs_bound_n=", n, ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300 
)
