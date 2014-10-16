pcaModelAlt <-
  function (node, main.data, tree, what = 'local', corC = TRUE, ...,
            filter.sv = 1e2)
  {
    options (contrasts = c('contr.sum', 'contr.poly'))
    require (mvtnorm)
    require (Morphometrics)
    stan.data <- list()
    subtree <- extract.clade(tree, node)
    stan.data $ C <- vcvPhylo (subtree, anc.nodes = TRUE)
    if (corC)
      stan.data $ C <- solve (cov2cor (stan.data $ C))
    else
      stan.data $ C <- solve (stan.data $ C)
    subset <- match (subtree $ tip.label, names(main.data))
    if (what == 'local')
      {
        raw.data <- llply (main.data [subset], function (L) L $ local)
        raw.data.df <- ldply (main.data [subset], function (L) L $ local)
      }
    if (what == 'ed')
      {
        raw.data <- llply (main.data [subset], function (L) L $ ed)
        raw.data.df <- ldply (main.data [subset], function (L) L $ ed)
      }
    formula.mlW <- paste (colnames (raw.data.df) [-1], collapse = ', ')
    formula.mlW <- paste ('cbind(', formula.mlW, ') ~ .id', sep = '')
    formula.mlW <- as.formula(formula.mlW)
    model.mlW <- lm (formula.mlW, data = raw.data.df)
    stan.data <-
      within (stan.data,
              {
                m <- length (subtree $ tip.label)
                k <- ncol (raw.data [[1]])
                n_pc <- ifelse (k >= m, m - 1, k)
                ni <- laply (raw.data, function (L) nrow (L))
                ni_max <- max (ni)
                X <- array (0, c(m, ni_max, k))
                for (i in 1:m)
                  X [i, 1:ni[i], ] <- raw.data [[i]]
              })
    aux <- list()
    aux <- within(aux,
                  {
                    mlW <- CalculateMatrix (model.mlW)
                    eigenW <- eigen (mlW)
                    raw.rotated <-
                      llply (raw.data,
                             function (L)
                             t (t (eigenW $ vectors [,1:stan.data $ n_pc]) %*% t (L)))
                    pca.means <- laply (raw.rotated, colMeans)
                    raw.means <- laply (raw.data, colMeans)
                    raw.gmean <- colMeans (raw.means)
                    B.pca <- var(pca.means)
                    jitter.mean <- rmvnorm(stan.data $ m, sigma = filter.sv * mlW)
                  })
    stan.data $ eigen_Wml <- aux $ eigenW $ vectors [,1:stan.data $ n_pc]
    start.values <- list()
    start.values $ c1 <- list()
    start.values $ c1 <- 
      within (start.values $ c1,
              {
                terminal <- aux $ raw.means + aux $ jitter.mean
                root <- c(rmvnorm(1, aux $ raw.gmean, filter.sv * aux $ mlW))
                ancestor <- rmvnorm(stan.data $ m - 2, aux $ raw.gmean,
                                    sigma = filter.sv * aux $ mlW)
                GammaW <- diag(stan.data $ k)
                sigmaW <- sqrt(filter.sv * diag(aux $ mlW))
                GammaB <- t(chol(cov2cor (aux $ B.pca)))
                sigmaB <- sqrt(filter.sv * diag(aux $ B.pca))
              })
    model.fit <- stan (file = '../Stan/pcaSigma_Anc.stan', chains = 1, 
                       data = stan.data, init = start.values, ...)
    return (list ('model' = model.fit, 'data' = stan.data, 'start' = start.values,
                  'aux' = aux))
  }
