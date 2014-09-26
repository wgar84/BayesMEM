require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 8)
#registerDoMC (cores = 60)
require (reshape2)
require (ggplot2)
require (rstan)
require (phytools)
require (geiger)
require (mvtnorm)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

source ('../FuncR/mainModel.R')

## RAxML.tree <- read.tree ('raxml.nm.tre')
## RAxML.tree <- treedata (RAxML.tree, OneDef, TRUE) $ phy
## Tree [[5]] <- RAxML.tree
## names (Tree) [5] <- 'raxml.nm'
## save (Tree, file = '../../Databases/Tree.RData')
## plot (RAxML.tree, direction = 'upwards')

plot (extract.clade (Tree [[5]], 138), direction = 'upwards')
nodelabels()

ancTest <- ancModel(138, OneDef, Tree [[5]], 'local',
                    list ('mean' = c(5.2, rep (0, times = 38)),
                          'vcv' = post.vcv $ ss.grand.mean),
                    pars = c('ancestor', 'root'),
                    warmup = 500, iter = 1000, thin = 5)

ancTest.post <- list ()
ancTest.post $ ext <- extract (ancTest)

dimnames (ancTest.post $ ext $ ancestor) <- list (1:100, as.character (14:23),
                                                  names (OneDef [[1]] $ mean))

names (dimnames (ancTest.post $ ext $ ancestor)) <- c ('iteration', 'node', 'trait')

ancTest.post $ anc.df <- melt (ancTest.post $ ext $ ancestor)
ancTest.post $ anc.df $ node <- as.character(ancTest.post $ anc.df $ node)

ggplot (ancTest.post $ anc.df, aes (y = value, x = iteration, colour = node)) +
  geom_point () + facet_wrap (~ trait) + geom_line()
