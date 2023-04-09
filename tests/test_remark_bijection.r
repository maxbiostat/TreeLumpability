library(PhyloMarkovChains)
library(phangorn)

n <- 7
all.splits <- sapply(allSplits(n + 1),  paste, collapse = "")
J <- length(all.splits)

all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)

compute <- FALSE

if (compute) {
  comp.time <- system.time(SPR.mat <-
                             as.matrix(rspr_matrix(all.trees) > 0) + 0)
  save(SPR.mat,
       file = paste0("../saved_data/SPR_matrix_n=", n, ".RData"))
} else{
  load(paste0("../saved_data/SPR_matrix_n=", n, ".RData"))
}

####
# Neighbourhood stuff
incidence.mat <- as.matrix(apply(SPR.mat, 2,
                                 function(x)
                                   as.numeric(x == 1)))

n.sizes.inc <- colSums(incidence.mat)
n.sizes.decomp <- unlist(lapply(all.trees,
                                function(t)
                                  neighbourhood_size(t)$n_size))

all.equal(n.sizes.inc, n.sizes.decomp)

complete.Splits <- parallel::mclapply(all.trees, get_clades,
                                      mc.cores = 6)

complete.indexes <- parallel::mclapply(complete.Splits,
                                       function(x)
                                         charmatch(x, all.splits),
                                       mc.cores = 6)

complete.clade.matrix <-
  do.call(rbind,
          parallel::mclapply(
            complete.indexes,
            clade_indicators,
            L = J,
            mc.cores = 6
          ))
colnames(complete.clade.matrix) <- all.splits
colMeans(complete.clade.matrix)
# colSums(complete.clade.matrix)
# cladeCorrelation::pclade(2, n)

clade <- "14"
clade.ind <- clade_index(x = clade, all.splits = all.splits)

S1.c.pos <- which(complete.clade.matrix[, clade.ind] == 1)

in.neighbourhood.pos <- apply(incidence.mat, 2,
                              function(x)
                                which(x == 1))

A1.c.pos <- lapply(in.neighbourhood.pos,
                   function(x)
                     intersect(x, S1.c.pos))

A1.c.sizes <- unlist(lapply(A1.c.pos, length))

sizes.dt <- data.frame(
  tree = paste0("tree_", 1:K),
  neighbourhood_size = n.sizes.inc,
  A1c_size  = A1.c.sizes
)

hasclade <-  S1.c.pos

clade.nice <- paste0("t", strsplit(clade, "")[[1]])

parts <- lapply(all.trees[hasclade],
                safe_f_c, clade = clade.nice)

Nxprimes <- unlist(lapply(parts,
                          function(x)
                            neighbourhood_size(x$x_prime)$n_size))
Ntauxcs <- unlist(lapply(parts,
                         function(x)
                           neighbourhood_size(x$tau_x_c)$n_size))

cbind(sizes.dt[hasclade,], Nxprimes + Ntauxcs)