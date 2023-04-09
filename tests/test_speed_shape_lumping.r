library(phangorn)
library(PhyloMarkovChains)

##############
n <- 5
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

### Building the MH transition matrix
incidence.mat <- SPR.mat
neighbourhood.sizes <- colSums(incidence.mat)

neigh.ratios <- matrix(NA, nrow = K, ncol = K)
diag(neigh.ratios) <- 0
for (i in 1:K) {
  neigh.ratios[i,] <-
    sapply(neighbourhood.sizes[i] / neighbourhood.sizes,
           function(x)
             min(1, x)) /
    neighbourhood.sizes[i]
}

trans.mat.MH <- neigh.ratios * incidence.mat
diag(trans.mat.MH) <- 1 - rowSums(trans.mat.MH)

get_shapes_slow <- function(all_trees) {
  K <- length(all_trees)
  
  Shape <- matrix(NA, nrow = K, ncol = K)
  
  for (i in 1:K) {
    for (j in 1:K) {
      Shape[i, j] <- all.equal.phylo(
        all_trees[[i]],
        all_trees[[j]],
        use.edge.length = FALSE,
        use.tip.label = FALSE
      )
    }
  }
  Shape_num <- matrix(as.numeric(as.logical(Shape)),
                      nrow = K,
                      ncol = K)
  return(Shape_num)
}

vec2symMat <- function (x, diag = TRUE, byrow = FALSE) {
  # stolen from https://github.com/mikewlcheung/metasem/blob/17542243fc3f69c1b6d3c53ec68c00f4f8bbb81f/R/vec2symMat.R
  m <- length(x)
  d <- if (diag)
    1
  else
    - 1
  n <- floor((sqrt(1 + 8 * m) - d) / 2)
  if (m != n * (n + d) / 2)
    stop("Cannot make a square matrix as the length of \"x\" is incorrect.")
  mat <- diag(n)
  
  ## Row major
  if (byrow) {
    mat[upper.tri(mat, diag = diag)] <- x
    index <- lower.tri(mat)
    mat[index] <- t(mat)[index]
  } else {
    ## Column major: default behavior
    mat[lower.tri(mat, diag = diag)] <- x
    # Just mirroring the matrix, exclude the diagonals
    ## mat[upper.tri(mat, diag=FALSE)] <- mat[lower.tri(mat, diag=FALSE)]
    ## Corrected a bug
    index <- upper.tri(mat)
    mat[index] <- t(mat)[index]
  }
  mat
}

get_shapes_fast <- function(all_trees, ncores = 6) {
  K <- length(all_trees)
  grid <- subset(expand.grid(1:K, 1:K), Var1 < Var2)
  N <- nrow(grid)
  comps <- unlist(parallel::mclapply(1:N,
                                     function(k) {
                                       i <- grid[k, 1]
                                       j <- grid[k, 2]
                                       all.equal.phylo(all_trees[[i]],
                                                       all_trees[[j]],
                                                       use.edge.length = FALSE,
                                                       use.tip.label = FALSE)
                                     }, mc.cores = ncores))
  comps <- as.numeric(comps)
  Shape_num <- vec2symMat(x = comps, diag = FALSE, byrow = TRUE)
  diag(Shape_num) <- 1
  return(Shape_num)
}

SM1 <- get_shapes_slow(all.trees)
colSums(SM1)
unique(colSums(SM1))

SM2 <- get_shapes_fast(all.trees)
colSums(SM2)
unique(colSums(SM2))

bench::mark(get_shapes_slow(all.trees),
            get_shapes_fast(all.trees),
            memory = FALSE)
