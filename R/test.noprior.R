require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 12)
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

noTest <-
  alply(1:12, 1, function (i)
        failMainModel(138, OneDef, Tree [[1]], 'local', corC = TRUE, 
                      list ('mean' = OneDef [[i+26]] $ mean,
                            'vcv' = 10*OneDef[['Callithrix_kuhlii']] $ml.vcv),
                      model = 'oneSigma_Anc',
                      pars = c('terminal', 'root', 'ancestor',
                        'SigmaW', 'SigmaB'),
                      control = list ('chain_id' = i),
                      warmup = 500, iter = 1000, thin = 5),
        .parallel = TRUE)

save (noTest, file = 'noTest.RData')





