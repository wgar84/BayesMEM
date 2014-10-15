pcaModel <-
  function (node, main.data, tree, what = 'local', corC = TRUE, ...)
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
        formula.mlW <- paste (colnames (raw.data.df) [-1], collapse = ', ')
        formula.mlW <- paste ('cbind(', formula.mlW, ') ~ .id', sep = '')
        formula.mlW <- as.formula(formula.mlW)
        model.mlW <- lm (formula.mlW, data = raw.data.df)
        mlW <- CalculateMatrix (model.mlW)
        eigenW <- eigen (mlW)
      }
    if (what == 'ed')
      {
        raw.data <- llply (main.data [subset], function (L) L $ ed)
        raw.data.df <- ldply (main.data [subset], function (L) L $ ed)
        formula.mlW <- paste (colnames (raw.data.df) [-1], collapse = ', ')
        formula.mlW <- paste ('cbind(', formula.mlW, ') ~ .id', sep = '')
        formula.mlW <- as.formula(formula.mlW)
        model.mlW <- lm (formula.mlW, data = raw.data.df)
        mlW <- CalculateMatrix (model.mlW)
        eigenW <- eigen (mlW)
      }
    stan.data <-
      within (stan.data,
              {
                m <- length (subtree $ tip.label)
                k <- ncol (raw.data [[1]])
                k <- ifelse (k >= m, m - 1, k)
                ni <- laply (raw.data, function (L) nrow (L))
                ni_max <- max (ni)
                X <- array (0, c(m, ni_max, k))
                raw.rotated <-
                  llply (raw.data, function (L) t (eigenW $ vectors [,1:k]) %*% L)
                raw.means <- laply (raw.rotated, colMeans)
                for (i in 1:m)
                  X [i, 1:ni[i], ] <- raw.rotated [[i]]
              })
    start.values <- list()
    start.values $ c1 <- list()
    start.values $ c1 <- 
      within (start.values $ c1,
              {
                W.pca <- diag (stan.data $ eigenW $ values)
                B.pca <- var(stan.data $ raw.means)
                grand.mean <- colMeans(stan.data $ raw.means)
                BW.vcv <- var(stan.data $ raw.means)
                terminal <- stan.data $ raw.means +
                  rmvnorm(m, sigma = W.pca)
                root <- rmvnorm(1, grand.mean, W.pca)
                ancestor <- rmvnorm(m - 2, grand.mean, sigma = W.pca)
                GammaW <- diag(stan.data$k)
                sigmaW <- sqrt(diag (W.pca))
                GammaB <- t(chol(B.pca))
                sigmaB <- sqrt(diag(B.pca))
              })
    model.fit <- stan (file = '../Stan/oneSigma_Anc.stan', chains = 1, 
                       data = stan.data, init = start.values, ...)
    return (model.fit)
  }
