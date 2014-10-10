PostDeltaZ <- function (extracted, tree, beta = FALSE)
  {
    #### 3dSandwich (alpha, Xbar)
    n.taxa <- length (tree $ tip.label)
    names.taxa <- c(tree $ tip.label, as.character (n.taxa + 1:tree $ Nnode))
    m <- length (tree $ tip.label)
    anc.states.tree <- array (0, dim (extracted $ Xbar) + c(0, 1, 0))
    anc.states.tree [, 1:m, ] <- extracted $ Xbar[, 1:m, ]
    anc.states.tree [, m+1, ] <- extracted $ alpha
    anc.states.tree [, (m+2):(2*m - 1), ] <- extracted $ Xbar[, (m+1):(2*m - 2), ]
    
    DeltaZ <-
      aaply (tree $ edge, 1, function (index)
             {
               o <- index [1] ### inicial
               t <- index [2] ### final
               anc.states.tree [, t, ] - anc.states.tree [, o, ]
             })
    DeltaZ <- aperm (DeltaZ, c(2, 1, 3))
    names (dimnames (DeltaZ)) <- c('iterations', 'edge', 'trait')
    dimnames (DeltaZ) [c(1, 3)] <- list (as.character (1:100),
                                         dimnames (extracted $ Xbar) [[3]])
    dimnames (DeltaZ) [[2]] <- paste (names.taxa [tree $ edge [, 1]],
                                      names.taxa [tree $ edge [, 2]], sep = '-')
    if (beta)
      {
        Beta <- aaply (1:(dim (DeltaZ) [1]), 1, function (i)
                       solve (extracted $ Sigma [i, , ]) %*% t(DeltaZ [i, , ]))
        Beta <- aperm (Beta, c(1, 3, 2))
        dimnames (Beta) <- dimnames (DeltaZ)
        names (dimnames (Beta)) <-  names (dimnames (DeltaZ))
        return (list ('DeltaZ' = DeltaZ, 'Beta' = Beta))
      }
      
    return (DeltaZ)
    
  }

