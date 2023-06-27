library(PhyloMarkovChains)
library(phangorn)
library(tidyverse)
library(markovchain)
library(rwty)
source("aux_asdsf.r")
######
n <- 4
load(paste0("saved_data/MH_matrix_n=", n, ".RData"))
all.trees <- allTrees(n, rooted = TRUE)
K <- length(all.trees)
States <- paste0("tree_", 1:K)

phylo.MH <- new(
  "markovchain",
  states = States,
  transitionMatrix = trans.mat.MH,
  name = "rootedMH"
)
Niter <- 1E3
Nchains <- 100
ncores <- 4

run_once <- function(i){
  ini.state <- sample(States, 1)
  raw.chain <- markovchain::markovchainSequence(n = Niter,
                                                      markovchain = phylo.MH,
                                                      t0 = ini.state)
  return(raw.chain)
}

raw.chains <- parallel::mclapply(1:Nchains,
                                 run_once,
                                 mc.cores = min(ncores, Nchains))

process_chain <- function(i){
  chain <- all.trees[match(raw.chains[[i]], States)]
  class(chain) <- "multiPhylo"
  return(chain)
}

chains <- parallel::mclapply(1:Nchains,
                             process_chain,
                             mc.cores = min(ncores, Nchains))

save(chains, file = paste0("saved_data/MH_n=", n,
                    "_Nchains=", Nchains,
                    ".RData"))

asdsf.stuff <- compute_asdsf(chains)

asdsf.df <- asdsf.stuff$asdsf_table

write.csv(asdsf.df,
          file = paste0("saved_data/MH_n=", n,
                           "_Nchains=", Nchains,
                           ".csv"),
          row.names = FALSE)

tail(asdsf.df)