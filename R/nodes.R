require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
#registerDoMC (cores = 8)
registerDoMC (cores = 60)
require (reshape2)
require (ggplot2)
require (rstan)
require (phytools)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

source ('../FuncR/mainModel.R')

Nodes <- list ()

### fazer modelo oneSigma_Anc
Nodes $ intermediate <- c(111, 112, 113, 115, 117, 121, 122, 123,
                          127, 129, 131, 132, 133, 134, 137, 138,
                          139, 141, 142, 143, 144, 147, 149, 150,
                          151, 152, 153, 154, 155, 156, 157, 158,
                          159, 161, 162, 163, 164, 175, 176, 177,
                          178, 180, 181, 182, 184, 185, 188, 189,
                          198, 199, 200, 201, 202, 203, 206, 207,
                          210, 211, 216)

### fazer modelo multiSigma_noAnc
Nodes $ terminal <- c(114, 115, 117, 123, 128, 129, 133, 134,
                      141, 149, 152, 157, 164, 175, 181, 199, 210)

Nodes $ inter.model <-
  alply (Nodes $ intermediate, 1, mainModel,
         main.data = OneDef, tree = Tree [[1]], what = 'local',
         prior.list = list ('mean' = OneDef [['Cebus_apella']] $ mean,
           'vcv' = post.vcv $ ss.grand.mean),
         model = 'oneSigma_Anc',
         pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
         warmup = 500, iter = 1000, thin = 5, .parallel = TRUE)

save (Nodes, file = 'Nodes.RData')

Nodes $ term.model <-
  alply (Nodes $ terminal, 1, mainModel,
         main.data = OneDef, tree = Tree [[1]], what = 'local',
         prior.list = list ('mean' = OneDef [['Cebus_apella']] $ mean,
           'vcv' = post.vcv $ ss.grand.mean),
         model = 'multiSigma_noAnc',
         pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
         warmup = 500, iter = 1000, thin = 5, .parallel = TRUE)

save (Nodes, file = 'Nodes.RData')
