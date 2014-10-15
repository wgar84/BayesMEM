require (mvtnorm)
require (phytools)

mainModel <-
  function (node, main.data, tree, what = 'local',  
            prior.list, model = 'oneSigma_noAnc', corC = TRUE, ...)
  {
    stan.data <- list()
    subtree <- extract.clade(tree, node)
    stan.data $ C <-
      vcvPhylo (subtree,
                anc.nodes = ifelse (grepl ('_Anc', model), TRUE, FALSE))
    if (corC)
      stan.data $ C <- solve (cov2cor (stan.data $ C))
    else
      stan.data $ C <- solve (stan.data $ C)
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
                                        #priorS <- prior.list $ vcv
                                        #priorX <- prior.list $ mean
              })
    start.values <- list()
    start.values $ c1 <- list()
    start.values $ c1 <- 
      within (start.values $ c1,
              {
                terminal <- rmvnorm(stan.data $ m, prior.list $ mean, prior.list $ vcv)
                root <- as.vector (rmvnorm(1, prior.list $ mean, prior.list $ vcv))
                if (grepl('multi', model))
                  {
                    Gamma <- array (t (chol (cov2cor (prior.list $ vcv))),
                                    c(stan.data$k, stan.data$k, stan.data$m))
                                        # cholesky_factor_corr[k] Gamma[m] (m x k x k)
                    GammaW <- aperm (Gamma, c(3, 1, 2))
                    sigmaW <- t (array (sqrt (diag (prior.list $ vcv)),
                                       c(stan.data$k, stan.data$m)))
                  }
                else
                  {
                    GammaW <- t (chol (cov2cor (prior.list $ vcv)))
                    sigmaW <- sqrt (diag (prior.list $ vcv))
                  }
                if (grepl ('_Anc', model))
                  ancestor <-
                    rmvnorm (stan.data$m - 2, prior.list $ mean,
                             prior.list $ vcv)
                if (grepl ('alt', model))
                  {
                    drift <- 1
                    Gamma_beta <- t (chol (cov2cor (prior.list $ vcv)))
                    sigma_beta <- sqrt (diag (prior.list $ vcv))
                  }
                if (grepl ('pcaS', model))
                  {
                    dim.flag <- ifelse (stan.data $ k >= stan.data $ m,
                                        stan.data $m - 1, stan.data $k)
                    GammaB <- t (chol (diag (dim.flag)))
                    sigmaB <- rep (1, times = dim.flag)
                  }
                else
                  {
                    GammaB <- t (chol (cov2cor (prior.list $ vcv)))
                    sigmaB <- sqrt (diag (prior.list $ vcv))
                  }
              })
    model.file <- paste ('../Stan/', model, '.stan', sep = '')
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
                    change.m <- length (which (dim (L) == stan.data $ m)) != 0
                    if (change.m)
                      for (i in which (dim (L) == stan.data $ m))
                        {
                          names (dimnames (L)) [i] <- 'node'
                          dimnames (L) [[i]] <- subtree $ tip.label
                        }
                    if (ifelse (grepl ('_Anc', model), TRUE, FALSE))
                      {
                        change.mm <- length (which (dim (L) == stan.data $ m - 2)) != 0
                        if (change.mm)
                          for (i in which (dim (L) == stan.data $ m - 2))
                           {
                             names (dimnames (L)) [i] <- 'node'
                             dimnames (L) [[i]] <-
                               as.character (stan.data $ m + (2:(stan.data $ m - 1)))
                           }
                     }
                    L
                  })
    return (list ('model' = model.fit, 'extract' = ext))
  }


