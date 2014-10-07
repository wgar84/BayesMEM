require (mvtnorm)
require (phytools)

mainModel <-
  function (node, main.data, tree, what = 'local', initial.state = 'stan', 
            prior.list, model = 'oneSigma_noAnc', ...)
  {
    stan.data <- list()
    subtree <- extract.clade(tree, node)
    stan.data $ C <- vcvPhylo (subtree,
                               anc.nodes = ifelse (grepl ('_Anc', model), TRUE, FALSE))
    subset <- match (subtree $ tip.label, names(main.data))
    if (what == 'local')
      raw.data <- llply (main.data [subset], function (L) L $ local)
    if (what == 'ed')
      raw.data <- llply (main.data [subset], function (L) L $ ed)
    stan.data <-
      within (stan.data,
              {
                k <- ncol (raw.data [[1]])
                m <- length (subtree $ tip.label)
                ni <- laply (raw.data, function (L) nrow (L))
                ni_max <- max (ni)
                X <- array (0, c(m, ni_max, k))
                for (i in 1:m)
                  X [i, 1:ni[i], ] <- raw.data [[i]]
                priorS <- prior.list $ vcv
                priorX <- prior.list $ mean
              })
    start.values <- list()
    start.values $ c1 <- list()
    start.values $ c1 <- 
      within (start.values $ c1,
              {
                Xbar <- rmvnorm(stan.data $ m, prior.list $ mean, prior.list $ vcv)
                alpha <- as.vector (rmvnorm(1, prior.list $ mean, prior.list $ vcv))
                if (grepl('multi', model))
                  {
                    Gamma <- array (t (chol (cov2cor (prior.list $ vcv))),
                                       c(stan.data$k, stan.data$k, stan.data$m))
                    ### cholesky_factor_corr[k] Gamma[m] (m x k x k)
                    Gamma <- aperm (Gamma, c(3, 1, 2))
                    sigma <- t (array (diag (prior.list $ vcv),
                                       c(stan.data$k, stan.data$m)))
                  }
                else
                  {
                    Gamma <- t (chol (cov2cor (prior.list $ vcv)))
                    sigma <- diag (prior.list $ vcv)
                  }
                if (grepl ('_Anc', model))
                  Xbar <-
                      rbind (Xbar, rmvnorm (stan.data$m - 2, prior.list $ mean,
                                               prior.list $ vcv))
                Gamma_bm <- t (chol (cov2cor (prior.list $ vcv)))
                sigma_bm <- diag (prior.list $ vcv)
              })
    model.file <- paste ('../Stan/', model, '.stan', sep = '')
    if (initial.state == 'stan')
      model.fit <- stan (file = model.file, chains = 1, 
                         data = stan.data, ...)
    else
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
                          dimnames (L) [[i]] <- colnames (raw.data [[1]])
                        }
                    if (ifelse (grepl ('_Anc', model), TRUE, FALSE))
                      {
                        change.m <- length (which (dim (L) == 2*stan.data $ m - 2)) != 0
                        if (change.m)
                          for (i in which (dim (L) == 2*stan.data $ m - 2))
                            {
                              names (dimnames (L)) [i] <- 'node'
                              dimnames (L) [[i]] <-
                                c(subtree $ tip.label,
                                  as.character (stan.data $ m + (2:(stan.data $ m - 1))))
                            }
                      }
                    else
                      {
                        change.m <- length (which (dim (L) == stan.data $ m)) != 0
                        if (change.m)
                          for (i in which (dim (L) == stan.data $ m))
                            {
                              names (dimnames (L)) [i] <- 'node'
                              dimnames (L) [[i]] <- subtree $ tip.label
                            }
                      }
                    L
                  })
    return (list ('model' = model.fit, 'extract' = ext))
  }


