DiagQuantile <- function (real.data, post.means, post.vcvs, at = c(0, 0.5, 1))
  {
    ### return quantile diagnostics for MCMC runs of multinormal parameters
    real.quant <- aaply (real.data, 2, quantile, probs = at)
    pops <- aaply (1:nrow (post.means), 1, function (i)
                   rmvnorm(nrow (real.data),
                           mean = post.means [i, ],
                           sigma = post.vcvs [i, , ]))
    post.quant <- aaply (pops, c(1, 3), quantile, probs = at)
    dimnames (post.quant) [-1] <- dimnames (real.quant)
    names (dimnames (post.quant)) <- c('sample', 'variable', 'quantile')
    names (dimnames (real.quant)) <- c('variable', 'quantile')
    post.quant <- melt (post.quant)
    post.quant $ variable <- as.factor (post.quant $ variable)
    post.quant $ quantile <- as.factor (post.quant $ quantile)
    real.quant <- melt (real.quant)
    real.quant $ variable <- as.factor (real.quant $ variable)
    real.quant $ quantile <- as.factor (real.quant $ quantile)
    Plot <-
      ggplot (post.quant, aes (x = value)) +
        geom_histogram (aes (fill = quantile), position = 'identity', alpha = 0.5) +
          facet_wrap (~ variable, scales = 'free_x') +
            geom_vline (aes (xintercept = value, colour = quantile), real.quant)
    return (Plot)
  }