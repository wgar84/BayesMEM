require (ape)
require (MCMCglmm)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
#registerDoMC (cores = 8)
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

All <- list ()

All $ oneSigma <- list()
All $ oneSigma <-
  within (All $ oneSigma,
          {
            k <- 39
            m <- 109
            C <- vcvPhylo (Tree [[1]], anc.nodes = FALSE)
            ni <- Aux $ sample.size 
            ni_max <- max (ni)
            X <- array (0, c (m, ni_max, k))
            for (i in 1:m)
              X [i, 1:ni [i], ] <-
                OneDef [[i]] $ local
            priorS <- post.vcv $ ss.grand.mean
            priorX <- OneDef [['Cebus_apella']] $ mean
          })

All $ oneSigma.start <-
  list('c1' = list(
         'Xbar' = t(array (rep (OneDef [['Cebus_apella']] $ mean, 109), c(39, 109))),
         'alpha' = OneDef [['Cebus_apella']] $ mean,
         'Gamma_bm' = t (chol (post.vcv $ ss.grand.mean)),
         'Gamma' = t (chol (post.vcv $ ss.grand.mean)),
         'sigma_bm' = diag (post.vcv $ ss.grand.mean),
         'sigma' = diag (post.vcv $ ss.grand.mean)))
         

All $ fit.oneSigma <- stan(file = '../Stan/oneSigma_cor.stan',
                           data = All $ oneSigma,
                           pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
                           init = All $ oneSigma.start,
                           warmup = 1, iter = 2, chains = 1,
                           control = list (refresh = 1))

save (All, file = 'All.RData')
