require (mvtnorm)
require (phytools)

ancModel <-
  function (node, main.data, tree, what = 'local',
            prior.list, ...)
  {
    stan.data <- list()
    subtree <- extract.clade(tree, node)
    stan.data $ C <- vcvPhylo (subtree, anc.nodes = TRUE)
    subset <- match (subtree $ tip.label, names(main.data))
    if (what == 'local')
      raw.data <- laply (main.data [subset], function (L) L $ mean)
    if (what == 'ed')
      raw.data <- laply (main.data [subset], function (L) L $ ed.mean)
    stan.data <-
      within (stan.data,
              {
                k <- ncol (raw.data)
                m <- nrow (raw.data)
                Xbar <- raw.data
                #priorS <- prior.list $ vcv
                #priorX <- prior.list $ mean
              })
    start.values <- list()
    start.values $ c1 <- list()
    start.values $ c1 <- 
      within (start.values $ c1,
              {
                ancestor <- rmvnorm(stan.data $ m - 2, prior.list $ mean, prior.list $ vcv)
                root <- as.vector (rmvnorm(1, prior.list $ mean, prior.list $ vcv))
                Gamma_bm <- t (chol (cov2cor (prior.list $ vcv)))
                sigma_bm <- diag (prior.list $ vcv)
              })
    model.file <- '../Stan/Ancestor.stan'
    model.fit <- stan (file = model.file, chains = 1, 
                       data = stan.data, init = start.values, ...)
    ext <- extract (model.fit)
    ext <- llply (ext, function (L)
                  {
                    change.k <- length (which (dim (L) == stan.data $ k)) != 0
                    if (change.k)
                      for (i in which (dim (L) == stan.data $ k))
                        {
                          names (dimnames (L)) [i] <- 'trait'
                          dimnames (L) [[i]] <- colnames (raw.data)
                        }
                    change.m <- length (which (dim (L) == 2*stan.data $ m - 2)) != 0
                    if (change.m)
                      for (i in which (dim (L) == 2*stan.data $ m - 2))
                        {
                          names (dimnames (L)) [i] <- 'node'
                          dimnames (L) [[i]] <-
                            c(subtree $ tip.label,
                              as.character (stan.data $ m + (2:(stan.data $ m - 1))))
                        }
                    L
                  })
    return (list ('model' = model.fit, 'extract' = ext))
  }


