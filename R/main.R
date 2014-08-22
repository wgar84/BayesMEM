require (ape)
require (MCMCglmm)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 8)
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

par (mfrow = c(2, 2))

vcv.phylo <- llply (Tree, vcvPhylo, anc.nodes = FALSE)

Callithrix <- list ()

### com anc.nodes = TRUE aqui e uma modificação no código, se pá é possível tirar os
### estados de caráter ancestral do modelo direto

Callithrix $ reduced <- list()
Callithrix $ reduced <-
  within (Callithrix $ reduced,
          {
            k <- 39
            m <- 7
            C <- vcvPhylo (extract.clade (Tree [[1]], 143), anc.nodes = FALSE)
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

Callithrix $ reduced.start <-
  list('c1' = list(
         'Xbar' = t(array (rep (OneDef [['Callithrix_kuhlii']] $ mean, 7), c(39, 7))),
         'alpha' = OneDef [['Callithrix_kuhlii']] $ mean,
         'Sigma_bm' = post.vcv $ ss.grand.mean,
         'Sigma' = aperm (post.vcv $ ss [, , 32, 1:7], c(3, 1, 2))),
       'c2' = list(
         'Xbar' = t(array (rep (OneDef [['Callithrix_penicillata']] $ mean, 7), c(39, 7))),
         'alpha' = OneDef [['Callithrix_jacchus']] $ mean,
         'Sigma_bm' = post.vcv $ ss [, , 32, 1],
         'Sigma' = aperm (post.vcv $ ss [, , 35, 8:14], c(3, 1, 2)))
       )

Callithrix $ fit.reduced <- stan(file = '../Stan/reduced.stan',
                                 data = Callithrix $ reduced,
                                 init = Callithrix $ reduced.start,
                                 warmup = 500, iter = 1500, chains = 2)

