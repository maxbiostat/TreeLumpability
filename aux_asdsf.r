#' Compute the average standard deviation of split (clade) frequencies from a set
#' of trees
#'
#' @param x a 'multiPhylo' object 
#' @param bnin the burnin in number of iterations
#' @param wsize the window size at which to take averages
#' @param minf minimum frequency for inclusion
#' @param ncores number of cores
#'
#' @return a list containing a frequency table and an ASDSF table
compute_asdsf <- function(x,
                          bnin = 0,
                          wsize = 20,
                          minf = .01,
                          ncores = 6) {
  sf.list <- parallel::mclapply(
    x,
    rwty:::get.slide.freq.table,
    burnin = bnin,
    window.size = wsize,
    gens.per.tree = 1,
    mc.cores = ncores
  )
  delta <- tibble::as_tibble(rwty:::get.asdsfs(sf.list, min.freq = minf))
  return(list(freq_table = sf.list,
              asdsf_table = delta))
}

plot.asdsf <- function(dat, log.y = TRUE) {
  require(ggplot2)
  ## copied from RWTY code base
  asdsf.plot <- ggplot(dat,
                       aes(x = as.numeric(as.character(Generation)))) +
    geom_line(aes(color = 14, y = min), linetype = 3) +
    geom_line(aes(color = 13, y = lower.95), linetype = 2) +
    geom_line(aes(color = 12, y = lower.75), linetype = 7) +
    geom_line(aes(color = 12, y = upper.75), linetype = 7) +
    geom_line(aes(color = 13, y = upper.95), linetype = 2) +
    geom_line(aes(color = 14, y = max), linetype = 3) +
    geom_ribbon(aes(
      ymin = min,
      ymax = lower.95,
      fill = 14
    ), alpha = 0.50) +
    geom_ribbon(aes(
      ymin = lower.95,
      ymax = lower.75,
      fill = 13
    ), alpha = 0.50) +
    geom_ribbon(aes(
      ymin = lower.75,
      ymax = upper.75,
      fill = 12
    ), alpha = 0.50) +
    geom_ribbon(aes(
      ymin = upper.75,
      ymax = upper.95,
      fill = 13
    ), alpha = 0.50) +
    geom_ribbon(aes(
      ymin = upper.95,
      ymax = max,
      fill = 14
    ), alpha = 0.50) +
    geom_line(aes(y = ASDSF)) +
    geom_point(aes(y = ASDSF)) +
    scale_color_viridis_b(begin = 0.2,
                          end = .9,
                          option = "D") +
    scale_fill_viridis_b(begin = 0.2,
                         end = .9,
                         option = "D") +
    expand_limits(y = 0) +
    theme(legend.position = "none") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    xlab("Iteration") +
    ggtitle("Average Standard Deviation of Split Frequencies")
  
  if (log.y == TRUE) {
    asdsf.plot <-
      asdsf.plot + scale_y_log10() + ylab("Log Standard Deviation of Split Frequencies")
  } else {
    asdsf.plot <-
      asdsf.plot + ylab("Log Standard Deviation of Split Frequencies")
  }
  return(list("asdsf.plot" = asdsf.plot))
}