DiagGelman <- function (extract.list)
  {
    nchain <- length (extract.list)
    if(nchain < 2)
      stop('supply at least two chains')
    theta.names <- names (extract.list [[1]])
    theta.names <- theta.names [-length(theta.names)]
    print(theta.names)
    theta.loa <- alply (theta.names, 1, function (theta)
                        laply (extract.list, function (L) L [[theta]]))
    names (theta.loa) <- theta.names
    theta.loa <- llply (theta.loa, function (L) {names (dimnames (L)) [1] <- 'chain'; L})
    some.var <- llply(theta.loa, function (L) aaply(L, c(1, 2), var))
    llply (some.var, dim)
  }
