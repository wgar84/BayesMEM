require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 3)
#registerDoMC (cores = 8)
#registerDoMC (cores = 60)
require (reshape2)
require (ggplot2)
require (grid)
require (gridExtra)
require (rstan)
require (phytools)
require (geiger)
require (mvtnorm)
require (boa)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')

#attach ('Test.RData')

test.post <- list()
test.post $ ext <- extract(Test)
test.post $ ext2 <- extract(Test2)

test.post $ subtree <- extract.clade(Tree [[5]], 138)

test.post $ ext <- llply (test.post $ ext, function (L)
                  {
                    change.k <- length (which (dim (L) == 39)) != 0
                    if (change.k)
                      for (i in which (dim (L) == 39))
                        {
                          names (dimnames (L)) [i] <- 'trait'
                          dimnames (L) [[i]] <- colnames (OneDef [[1]] $ local)
                        }
                    change.m <- length (which (dim (L) == 22)) != 0
                    if (change.m)
                      for (i in which (dim (L) == 22))
                        {
                          names (dimnames (L)) [i] <- 'node'
                          dimnames (L) [[i]] <-
                            c(test.post $ subtree $ tip.label, as.character (14:23))
                        }
                    L
                  })

test.post $ ext2 <- llply (test.post $ ext2, function (L)
                  {
                    change.k <- length (which (dim (L) == 39)) != 0
                    if (change.k)
                      for (i in which (dim (L) == 39))
                        {
                          names (dimnames (L)) [i] <- 'trait'
                          dimnames (L) [[i]] <- colnames (OneDef [[1]] $ local)
                        }
                    change.m <- length (which (dim (L) == 22)) != 0
                    if (change.m)
                      for (i in which (dim (L) == 22))
                        {
                          names (dimnames (L)) [i] <- 'node'
                          dimnames (L) [[i]] <-
                            c(test.post $ subtree $ tip.label, as.character (14:23))
                        }
                    L
                  })

test.quantile.diags <- list ()

test.quantile.diags $ Quant1 <- DiagQuantilePop(138, test.post $ ext, OneDef, Tree [[5]])
test.quantile.diags $ Quant2 <- DiagQuantilePop(138, test.post $ ext2, OneDef, Tree [[5]])

test.quantile.diags $ stackQuant <-
  arrangeGrob(test.quantile.diags $ Quant1 + labs (title = 'start with Cebus apella'),
              test.quantile.diags $ Quant2 + labs (title = 'start with Saguinus mystax'),
              ncol = 2)
ggsave(test.quantile.diags $ stackQuant, filename = 'quantile.convergence.jpg',
       width = 30, height = 16)

test.trace.diags1 <- DiagTrace(test.post $ ext, OneDef $ 'Cebus_apella' $ mean)

test.trace.diags1 $ alpha
test.trace.diags1 $ extant
test.trace.diags1 $ ancestor

test.trace.diags2 <- DiagTrace(test.post $ ext2, OneDef $ 'Saguinus_mystax' $ mean)

test.trace.diags2 $ alpha
test.trace.diags2 $ extant
test.trace.diags2 $ ancestor

