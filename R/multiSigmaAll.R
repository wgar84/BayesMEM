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

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

source ('../FuncR/mainModel.R')

fit.MN <- mainModel (110, OneDef, Tree [[1]], 'local',
                     list ('mean' = OneDef [['Cebus_apella']] $ mean,
                           'vcv' = post.vcv $ ss.grand.mean),
                     model = 'multiSigma_noAnc',
                     pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
                     warmup = 500, iter = 1000, thin = 5)

save (fit.MN, file = 'multiSigmaAll.RData')
