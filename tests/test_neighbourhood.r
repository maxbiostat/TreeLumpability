library(PhyloMarkovChains)

ns <- 3:20
plot(ns, Nsize_ratio(ns), type = "l")

ladT <- makeLadder(ntaxa = 5)
neighbourhood_size(ladT)
min_size(5)
