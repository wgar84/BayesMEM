PostDeltaZ <- function (extracted, tree, beta = FALSE)
  {
    #### 3dSandwich (alpha, Xbar)
    m <- length (tree $ tip.label)
    names.taxa <- c(tree $ tip.label, as.character (m + 1:tree $ Nnode))
    anc.states.tree <- array (0, dim (extracted $ terminal) + c(0, m - 1, 0))
    anc.states.tree [, 1:m, ] <- extracted $ terminal
    anc.states.tree [, m+1, ] <- extracted $ root
    anc.states.tree [, (m+2):(2*m - 1), ] <- extracted $ ancestor
    
    DeltaZ <-
      aaply (tree $ edge, 1, function (index)
             {
               o <- index [1] ### inicial
               t <- index [2] ### final
               anc.states.tree [, t, ] - anc.states.tree [, o, ]
             })
    DeltaZ <- aperm (DeltaZ, c(2, 1, 3))
    names (dimnames (DeltaZ)) <- c('iterations', 'edge', 'trait')
    dimnames (DeltaZ) [c(1, 3)] <- list (as.character (1:dim (DeltaZ) [1]),
                                         dimnames (extracted $ terminal) [[3]])
    dimnames (DeltaZ) [[2]] <- paste (names.taxa [tree $ edge [, 1]],
                                      names.taxa [tree $ edge [, 2]], sep = '-')
    if (beta)
      {
        Beta <- aaply (1:(dim (DeltaZ) [1]), 1, function (i)
                       solve (extracted $ SigmaW [i, , ], t(DeltaZ [i, , ])))
        Beta <- aperm (Beta, c(1, 3, 2))
        dimnames (Beta) <- dimnames (DeltaZ)
        names (dimnames (Beta)) <-  names (dimnames (DeltaZ))
        return (list ('DeltaZ' = DeltaZ, 'Beta' = Beta))
      }
      
    return (DeltaZ)
    
  }

