library(phangorn)
library(PhyloMarkovChains)
##############
###$$$ Aux functions

vec2symMat <- function(x, diag = TRUE, byrow = FALSE) {
  # stolen from https://github.com/mikewlcheung/metasem/blob/17542243fc3f69c1b6d3c53ec68c00f4f8bbb81f/R/vec2symMat.R
  m <- length(x)
  d <- if (diag)
    1
  else-1
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
##############
n <- 7
all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)

compute <- FALSE

if (compute) {
  comp.time <- system.time(SPR.mat <-
                             as.matrix(rspr_matrix(all.trees) > 0) + 0)
  save(SPR.mat,
       file = paste0("saved_data/SPR_matrix_n=", n, ".RData"))
} else{
  load(paste0("saved_data/SPR_matrix_n=", n, ".RData"))
}

### Building the MH transition matrix
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


Shape_num <- get_shapes_fast(all.trees)

dif_shapes <- unique.array(Shape_num)
DS <- nrow(dif_shapes)
shapes_each_total <- rowSums(dif_shapes)
Max_shape <- max(shapes_each_total)
Mat_shapes <- matrix(0, nrow = Max_shape, ncol = DS)


for (i in 1:DS) {
  Mat_shapes[(1:shapes_each_total[i]), i] <- which(apply(Shape_num, 1,
                                                         function(x)
                                                           all.equal(x[],
                                                                     dif_shapes[i,]) == TRUE))
  
}

one_shape_each <- array(NA, dim = DS)

for (i in 1:DS) {
  one_shape_each[i] <- Mat_shapes[1, i]
}

Mat_lump <- matrix(NA, nrow = DS, ncol = DS)

for (i in 1:DS) {
  x_S <- array(NA, dim = K)
  x_S <- trans.mat.MH[one_shape_each[i],]
  for (j in 1:DS) {
    Mat_lump[i, j] <- sum(x_S[Mat_shapes[, j]])
  }
}

Mat_lump
