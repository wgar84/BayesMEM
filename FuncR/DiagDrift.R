require (psych)

singleDrift <- function (W.mat, terminal)
  {
    options(contrasts = c('contr.sum', 'contr.poly'))
    m <- nrow (terminal)
    k <- ncol (terminal)
    W.eig.dec <- eigen (W.mat)
    if (k >= m)
      {
        W.pc <- W.eig.dec $ vectors [, 1:(m-1)]
        var.W <- W.eig.dec $ values [1:(m-1)]
      }
    else
      {
        W.pc <- W.eig.dec $ vectors
        var.W <- W.eig.dec $ values
      }
    B.on.W <- terminal %*% W.pc ### B.on.W será m x m - 1 ou m x k
    ### variance
    var.B.on.W <- aaply (B.on.W, 2, var)
    reg.test <- lm (log (var.B.on.W) ~ log (var.W))
    var.df <- data.frame (1:length (var.W), log (var.B.on.W), log (var.W))
    names (var.df) <- c('PC', 'B', 'W')
    ### correlation
    cor.B.on.W <- cor (B.on.W)
    cor.test <- cortest.bartlett (cor.B.on.W, n = m)
    return (list ('reg' = reg.test, 'cor' = cor.test, 'mat' = cor.B.on.W,
                  'var.df' = var.df))
  }

matrixDrift <- function (W.mat, B.mat, n.term)
  {
    options(contrasts = c('contr.sum', 'contr.poly'))
    m <- n.term
    k <- ncol (B.mat)
    W.eig.dec <- eigen (W.mat)
    if (k >= m)
      {
        W.pc <- W.eig.dec $ vectors [, 1:(m-1)]
        var.W <- W.eig.dec $ values [1:(m-1)]
      }
    else
      {
        W.pc <- W.eig.dec $ vectors
        var.W <- W.eig.dec $ values
      }
    B.on.W <- t(W.pc) %*% B.mat %*% W.pc 
    ### variance
    var.B.on.W <- diag (B.on.W)
    reg.test <- lm (log (var.B.on.W) ~ log (var.W))
    var.df <- data.frame (1:length (var.W), log (var.B.on.W), log (var.W))
    names (var.df) <- c('PC', 'B', 'W')
    ### correlation
    cor.B.on.W <- cov2cor (B.on.W)
    cor.test <- cortest.bartlett (cor.B.on.W, n = m)
    return (list ('reg' = reg.test, 'cor' = cor.test, 'mat' = cor.B.on.W,
                  'var.df' = var.df))
  }

DiagDrift <- function (extraction, mode = 'terminal')
  {
    fisherTrans <- function (r) 0.5 * (log (1 + r) - log (1 - r))
    W.ensemble <- extraction $ SigmaW
    terminal.ensemble <- extraction $ terminal
    B.ensemble <- extraction $ SigmaB
    n.it <- dim (W.ensemble) [1]
    km <- dim (terminal.ensemble) [3:2]
    ncor <- ifelse (km [1] >= km [2], km [2] - 1, km [1])
    df.bartlett <- (ncor^2 - ncor) / 2
    results.ensemble <-
      alply (1:n.it, 1, function (i)
             {
               if (mode == 'terminal')
                 ds.model <- singleDrift (W.ensemble [i, , ], terminal.ensemble [i, , ])
               else if (mode == 'matrix')
                 ds.model <- matrixDrift (W.ensemble [i, , ], B.ensemble [i, , ], km [2])
               ds.model
             })
    regcor.ens <-
      ldply (results.ensemble, function (L) c (coef (L $ reg), L $ cor $ chisq))
    mat.ens <-
      laply (results.ensemble, function (L)
             abs (fisherTrans (L $ mat [lower.tri (L $ mat)])))
    var.df <- ldply (results.ensemble, function (L) L $ var.df)
    colnames (var.df) [1] <- 'iterations'
    s <- which (lower.tri(diag(ncor)), arr.ind = TRUE) [, 2:1]
    names.cor <- paste (s [, 1], s[, 2], sep = '-')
    colnames (mat.ens) <- names.cor
    names (dimnames (mat.ens)) <- c('Iteration', 'Correlation')
    mat.df <- melt (mat.ens)
    colnames (regcor.ens) <- c ('Iteration', 'Intercept', 'Slope',
                                'Correlation')

    var.df.mean <- ddply (var.df, .(PC), summarize,
                          meanW = mean (W),
                          meanB = mean (B))
    var.df.hull <- ddply (var.df, .(PC), summarize,
                          hullW = W [chull(W, B)],
                          hullB = B [chull(W, B)])
    
    chisq.crit.value <- qchisq(0.95, df.bartlett)
    
    Plots <- list()

    Plots $ BvsW <- ggplot (var.df.mean) +
      geom_text(aes (x = meanW, y = meanB, label = PC), size = 3, alpha = 1) +
        theme_minimal() +
          geom_polygon(aes (x = hullW, y = hullB, group = PC),
                       data = var.df.hull, alpha = 0.1) +
                         xlab(expression (paste ('log ',lambda[W]))) +
                           ylab(expression (paste ('log ',lambda[B]))) +
                             labs (title = 'Eigenvalue Regression')
    for (i in 1:nrow(regcor.ens))
      Plots $ BvsW <- Plots $ BvsW +
        geom_abline(intercept = regcor.ens[i, 'Intercept'],
                    slope = regcor.ens[i, 'Slope'], alpha = 0.05)
    
    Plots $ reg <- ggplot (regcor.ens) +
      geom_histogram(aes(x = Slope), position = 'identity') +
        theme_minimal() + geom_vline (xintercept = 1) +
          labs (title = 'Regression Test') +
            xlab (expression (paste ('Slope ', 'log', lambda[B]) %~%
                              paste ('log', lambda[W])))
    Plots $ mat <- ggplot (mat.df) + 
      geom_histogram(aes(x = value), position = 'identity') +
        facet_wrap (~ Correlation, scales = 'free') +
          theme_minimal() + geom_vline (xintercept = 0) +
            labs (title = 'Individual Correlation Test') +
              theme(axis.text = element_text(size = 4))+
                xlab ('absolute Fisher z-transformed correlation')
    Plots $ cor <- ggplot (regcor.ens) +
      geom_histogram(aes(x = Correlation), position = 'identity') +
        theme_minimal() + geom_vline (xintercept = chisq.crit.value) +
          labs (title = 'Joint Correlation Test') +
            xlab(expression(paste (chi^2, ' Value')))
    Plots
  }
