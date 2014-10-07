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
require (rstan)
require (phytools)
require (geiger)

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

## plot (extract.clade (Tree [[5]], 138), direction = 'upwards')
## nodelabels()

Test.prior <- mainModel(138, OneDef, Tree [[5]], 'local',
                  list ('mean' = OneDef [['Cebus_apella']] $ mean,
                        'vcv' = post.vcv $ ss.grand.mean),
                  model = 'oneSigma_Anc_prior',
                  pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
                  warmup = 500, iter = 1000, thin = 5)

save (Test.prior, file = 'Test.prior.RData')
