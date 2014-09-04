require (ape)
require (MCMCglmm)
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

All <- list ()

All $ fit.SA <- mainModel (110, OneDef, Tree [[1]], 'local',
                           list ('mean' = OneDef [['Cebus_apella']] $ mean,
                                 'vcv' = post.vcv $ ss.grand.mean),
                           model = 'oneSigma_Anc',
                           pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
                           warmup = 50, iter = 100)

All $ oneSigma <- list()
All $ oneSigma <-
  within (All $ oneSigma,
          {
            k <- 39
            m <- 109
            C <- vcvPhylo (Tree [[1]], anc.nodes = TRUE)
            ni <- Aux $ sample.size 
            ni_max <- max (ni)
            X <- array (0, c (m, ni_max, k))
            for (i in 1:m)
              X [i, 1:ni [i], ] <-
                OneDef [[i]] $ local
            priorS <- post.vcv $ ss.grand.mean
            priorX <- OneDef [['Cebus_apella']] $ mean
          })

save (All, file = 'All.RData')
