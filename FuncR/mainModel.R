require (mvtnorm)
require (phytools)

mainModel <-
  function (node, main.data, tree, what = 'local',
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
    model.fit <- stan (file = model.file, chains = 1, 
                       data = stan.data, init = start.values, ...)
    return (model.fit)
  }


