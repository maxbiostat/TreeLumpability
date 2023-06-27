library(PhyloMarkovChains)
library(phangorn)
library(tidyverse)
library(markovchain)
library(rwty)
source("aux_asdsf.r")
######
n <- 4
load(paste0("saved_data/MH_matrix_n=", n, ".RData"))

all.trees <- phangorn::allTrees(n, rooted = TRUE)
K <- length(all.trees)
States <- paste0("tree_", 1:K)

phylo.MH <- new(
  "markovchain",
  states = States,
  transitionMatrix = trans.mat.MH,
  name = "rootedMH"
)
Niter <- 1E3
Nchains <- 4
raw.chains <- vector(Nchains, mode = "list")
chains <- vector(Nchains, mode = "list")

for (i in 1:Nchains) {
  ini.state <- sample(States, 1)
  raw.chains[[i]] <- markovchain::markovchainSequence(n = Niter,
                                                      markovchain = phylo.MH,
                                                      t0 = ini.state)
  chains[[i]] <- all.trees[match(raw.chains[[i]], States)]
  class(chains[[i]]) <- "multiPhylo"
}

asdsf.stuff <- compute_asdsf(chains)

asdsf.df <- asdsf.stuff$asdsf_table
tail(asdsf.df)

plot(ASDSF ~ Generation, asdsf.df, xlim = c(0, 500))

the.plt <- plot.asdsf(dat = asdsf.df,
                      log.y = FALSE)
ggsave(
  plot = the.plt$asdsf.plot,
  filename = paste0("asdsf_n=", n, ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
