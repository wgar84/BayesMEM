require (ape)
require (MCMCglmm)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
#registerDoMC (cores = 8)
registerDoMC (cores = 10)
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

Callithrix <- list ()

Callithrix $ ancSigma <- list()
Callithrix $ ancSigma <-
  within (Callithrix $ ancSigma,
          {
            k <- 39
            m <- 7
            invC <- solve (vcvPhylo (extract.clade (Tree [[1]], 143), anc.nodes = TRUE))
            ni <- Aux $ sample.size [which (grepl ('Callithrix', names (OneDef)))]
            ni_max <- max (ni)
            X <- array (0, c (m, ni_max, k))
            for (i in 1:m)
              X [i, 1:ni [i], ] <-
                OneDef [[which (grepl ('Callithrix', names (OneDef))) [i]]] $ local
            id <- diag (k)
            zeroes <- rep(0, times = k)
            priorS <- post.vcv $ ss.grand.mean
            priorX <- OneDef [['Callithrix_kuhlii']] $ mean
          })

Callithrix $ ancSigma.start <-
  list('c1' = list(
         'Xbar' = t(array (rep (OneDef [['Callithrix_kuhlii']] $ mean, 12), c(39, 12))),
         'alpha' = OneDef [['Callithrix_kuhlii']] $ mean,
         'Sigma_bm' = post.vcv $ ss.grand.mean,
         'Sigma' = post.vcv $ ss [, , 32, 1]),
       'c2' = list(
         'Xbar' = t(array (rep (OneDef [['Callithrix_penicillata']] $ mean, 12),
           c(39, 12))),
         'alpha' = OneDef [['Callithrix_jacchus']] $ mean,
         'Sigma_bm' = post.vcv $ ss [, , 32, 1],
         'Sigma' = post.vcv $ ss [, , 35, 19])
       )

Callithrix $ fit.ancSigma <- stan(file = '../Stan/ancSigma.stan',
                                  data = Callithrix $ ancSigma,
                                  warmup = 1, iter = 1, chains = 1)

Callithrix $ fit.fun <- function() stan(fit = Callithrix $ fit.ancSigma,
                                        data = Callithrix $ ancSigma,
                                        iter = 2000, thin = 100, chains = 1)

Callithrix $ fit.parl <- foreach (i = 1:10) %dopar% Callithrix $ fit.fun()

save(Callithrix, file = 'Callithrix.RData')
