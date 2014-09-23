PostDeltaZ <- function (extracted, tree, beta = FALSE)
  {
    #### Xbar alpha
    m <- length (tree $ tip.label)
    anc.states.tree <- array (0, dim (extracted $ Xbar) + c(0, 1, 0))
    anc.states.tree [, 1:m, ] <- extracted $ Xbar[, 1:m, ]
    anc.states.tree [, m+1, ] <- extracted $ alpha
    anc.states.tree [, (m+2):(2*m - 1), ] <- extracted $ Xbar[, (m+1):(2*m - 2), ]
    
    DeltaZ <-
      aaply (tree $ edge, 1, function (index)
             {
               o <- index [1]
               t <- index [2] ### o terminal é sempre o segundo
               anc.states.tree [, t, ] - anc.states.tree [, o, ]
             })
    DeltaZ <- aperm (DeltaZ, c(2, 1, 3))
    
    if (beta)
      {
        Beta <- aaply (1:(dim (DeltaZ) [1]), 1, function (i)
                       solve (extracted $ Sigma [i, , ]) %*% t(DeltaZ [i, , ]))
        Beta <- aperm (Beta, c(1, 3, 2))
        return (list ('DeltaZ' = DeltaZ, 'Beta' = Beta))
      }
      
    return (DeltaZ)
      
  }

