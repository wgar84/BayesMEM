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

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

failMainModel <- failwith(NULL, mainModel)

Test.NoPrior <- alply(1:8, 1, function (i)
                      failMainModel(138, OneDef, Tree [[5]], 'local',
                                list ('mean' = OneDef [[i+26]] $ mean,
                                      'vcv' = post.vcv $ ss.grand.mean),
                                model = 'oneSigma_Anc', initial.state = 'R',
                                pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
                                warmup =1000, iter = 2000, thin = 10), .parallel = TRUE)

### uma falhou
save (Test.NoPrior, file = 'TestNO.RData')
