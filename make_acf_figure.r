### Necessary data (.tar.gz'ed) can be downloaded from 10.5281/zenodo.8168348

library(BinaryMarkovChains)
library(ggplot2)

source("aux.r")
get_data <- function(ntaxa, lazyr, index) {
  fpath1 <-
    paste0("~/DUMP/data_repo/derived_data/chains_n=", ntaxa, "/")
  fpath2 <- paste0(
    "cladematrix_trees_n=",
    ntaxa,
    "_rho=",
    lazyr,
    "_Niter=10000_replicate=",
    index,
    ".cmap"
  )
  
  obs.Tree.Data <- load_clade_matrix(paste0(fpath1, fpath2))
  obs.cmat <-  obs.Tree.Data$clade_indicators
  return(obs.cmat)
}

get_theo_acf <- function(y, lags) {
  alpha.hat <- BinaryMarkovChains::estimate_alpha(sample = y)
  p.hat <- mean(y)
  beta.hat <- BinaryMarkovChains::get_beta(alpha = alpha.hat,
                                           p = p.hat)
  rho.hat <- function(k) {
    BinaryMarkovChains::MC_autocorr(lag = k,
                                    alpha = alpha.hat,
                                    beta = beta.hat)
  }
  out <- tibble::tibble(lag = lags,
                        theo_rho = rho.hat(lags))
  return(out)
}
######
rho <- 0.9
rep <- 42
ns <- c(4, 5, 6, 7)

simus <- lapply(ns, function(nt) {
  get_data(ntaxa = nt,
           lazyr = rho,
           index = rep)
})

c1 <- "{t1,t2}" 
c2 <- "{t1,t2,t3}"
c1s <- lapply(simus, function(X)
  X[, c1])
c2s <- lapply(simus, function(X)
  X[, c2])

K <- length(c1s)
max_lag <- 50

acfs.c1 <- lapply(1:K, function(k) {
  x <- c1s[[k]]
  bacf <- acf(x, plot = FALSE, lag.max = max_lag)
  bacfdf <- tibble::tibble(with(bacf, data.frame(lag, acf)))
  bacfdf$n_taxa <- paste0("n = ", ns[k])
  bacfdf$clade <- c1
  return(bacfdf)
})

acfs.c2 <- lapply(1:K, function(k) {
  x <- c2s[[k]]
  bacf <- acf(x, plot = FALSE, lag.max = max_lag)
  bacfdf <- tibble::tibble(with(bacf, data.frame(lag, acf)))
  bacfdf$n_taxa <- paste0("n = ", ns[k])
  bacfdf$clade <- c2
  return(bacfdf)
})

complete.df <- rbind(do.call(rbind, acfs.c1),
                     do.call(rbind, acfs.c2))


theos.c1 <- lapply(1:K, function(k) {
  x <- c1s[[k]]
  res <- get_theo_acf(y = x, lags = 0:max_lag)
  res$n_taxa <- paste0("n = ", ns[k])
  res$clade <- c1
  return(res)
})

theos.c2 <- lapply(1:K, function(k) {
  x <- c2s[[k]]
  res <- get_theo_acf(y = x, lags = 0:max_lag)
  res$n_taxa <- paste0("n = ", ns[k])
  res$clade <- c2
  return(res)
})

theoACF <-  rbind(do.call(rbind, theos.c1),
                  do.call(rbind, theos.c2))

complete.df$clade <- factor(complete.df$clade,
                            levels = c(c1, c2))
theoACF$clade <- factor(theoACF$clade,
                            levels = c(c1, c2))

ciline <- qnorm((1 + .95) / 2) / sqrt(acf(c1s[[2]], plot = FALSE)$n.used)


############# Plotting

final.plot <- ggplot(data = complete.df,
                     mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_line(aes(x = lag, y = theo_rho),
            data = theoACF,
            colour = "red", linetype = "longdash") +
  facet_grid(n_taxa ~ clade, switch = "y") +
  scale_x_continuous("lag (k)") +
  scale_y_continuous("Autocorrelation") + 
  geom_hline(aes(yintercept = ciline),
             linetype = 3,
             color = 'darkblue') +
  geom_hline(aes(yintercept = -ciline),
             linetype = 3,
             color = 'darkblue') +
  theme_bw(base_size = 20) +
  theme(strip.text.y = element_text(angle = 0),
        strip.text.y.left = element_text(angle = 0))


ggsave(
  plot = final.plot,
  filename = paste0("plots/LazyMH_ACF_rho=", rho, ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
