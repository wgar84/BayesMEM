DiagQuantilePop <- function (node, extract.stanfit,
                             main.data, tree, what = 'local', at = c (0, 0.5, 1))
  {
    subtree <- extract.clade(tree, node)
    n.extant <- length (subtree $ tip.label)
    or.pops <- llply (main.data [match (subtree $ tip.label, names (main.data))],
                      function (L) L [[what]])
    or.quant <-
      laply (or.pops, function (L) aaply (L, 2, quantile, probs = at))

    or.sample.sizes <- laply (or.pops, nrow)
            
    dimnames (or.quant) [[1]] <- names (or.pops)

    sampled.means <- extract.stanfit $ Xbar [, 1:n.extant, ]
    sampled.vcv <- extract.stanfit $ Sigma

    sim.quant <-
      aaply(1:(dim (sampled.means) [2]), 1, function (i)
            aaply(1:(dim (sampled.means) [1]), 1, function (j)
                  aaply (rmvnorm(or.sample.sizes [i], sampled.means [j, i, ],
                                 sampled.vcv [j, , ]), 2, quantile, probs = at)))
    
    names (dimnames (sim.quant)) <- c ('otu', 'iteration', 'trait', 'quantile')
    names (dimnames (or.quant)) <- c ('otu', 'trait', 'quantile')
    dimnames (sim.quant) [c(1, 3:4)] <- dimnames (or.quant)

    or.quant.df <- melt (or.quant)
    sim.quant.df <- melt (sim.quant)

    Plot <-
      ggplot() +
        geom_point(mapping = aes(x = otu, y = value, colour = quantile),
                   data = or.quant.df, size = 2)
    Plot <- Plot +
      geom_point(mapping = aes(x = otu, y = value, colour = quantile),
                 data = sim.quant.df,
                 size = 1, alpha = 0.1) +
                   facet_wrap(~ trait, scales = 'free_y') +
                     scale_x_discrete(labels = abbreviate) +
                       theme(axis.text.x = element_text(angle = 90))
    Plot
  }

