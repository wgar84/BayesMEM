require (ape)
require (MCMCglmm)
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
Sys.setenv(MAKEFLAGS = '-j2')
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
            id <- diag (k)
            zeroes <- rep(0, times = k)
            priorS <- post.vcv $ ss.mean [, , 44]
            priorX <- OneDef [['Homo_sapiens']] $ mean
          })

All $ oneSigma.start <-
  list('c1' = list(
         'Xbar' = t(array (rep (OneDef [['Homo_sapiens']] $ mean, 109), c(39, 109))),
         'alpha' = OneDef [['Alouatta_belzebul']] $ mean,
         'Sigma_bm' = post.vcv $ ss.grand.mean,
         'Sigma' = post.vcv $ ss [, , 22, 1]),
       'c2' = list(
         'Xbar' = t(array (rep (OneDef [['Miopithecus_talapoin']] $ mean, 109),
           c(39, 109))),
         'alpha' = OneDef [['Gorilla_gorilla']] $ mean,
         'Sigma_bm' = post.vcv $ ss [, , 102, 1],
         'Sigma' = post.vcv $ ss [, , 98, 19])
       )

All $ fit.oneSigma <- stan(file = '../Stan/oneSigma_noAnc.stan',
                           data = All $ oneSigma,
                           init = All $ oneSigma.start [1],
                           warmup = 1, iter = 1, chains = 1)

All $ func.oneSigma <-
  function (chain)
  stan (fit = All $ fit.oneSigma, data = All $ oneSigma,
        init = chain, warmup = 500, iter = 1000, chains = 1)

All $ fit.oneSigma.par <-
  foreach (i = 1:2) %dopar% All $ func.oneSigma (All $ oneSigma.start [i])

save(All, file = 'All.RData')
