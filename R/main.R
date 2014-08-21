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

Callithrix $ phyloVCV <- vcvPhylo (extract.clade (Tree [[1]], 143), anc.nodes = FALSE)

Callithrix $ data <- list()

for (i in which (grepl ('Callithrix', names (OneDef))))
  Callithrix $ data [[length (Callithrix $ data) + 1]] <- OneDef [[i]] $ local

names (Callithrix $ data) <- rownames (Callithrix $ phyloVCV)

Callithrix $ sample.sizes <- laply (Callithrix $ data, nrow)
names (Callithrix $ sample.sizes) <- names (Callithrix $ data)

Callithrix $ stan.data <- list(
  'k' = length (Callithrix $ data),
  'p' = ncol (Callithrix $ data [[1]]),
  'C' = Callithrix $ phyloVCV,
  'geoffroyi' = Callithrix $ data [[1]],
  'kuhlii' = Callithrix $ data [[2]],
  'jacchus' = Callithrix $ data [[3]],
  'penicillata' = Callithrix $ data [[4]],
  'argentata' = Callithrix $ data [[5]],
  'humeralifera' = Callithrix $ data [[6]],
  'pygmaea' = Callithrix $ data [[7]],
  'priorV' = post.vcv $ ss.grand.mean,
  'priorB' = OneDef [['Callithrix_kuhlii']] $ mean
  )

Callithrix $ stan.start <-
  list('c1' = list (
         'B' = t(array (rep (OneDef [['Callithrix_kuhlii']] $ mean, 7), c(39, 7))),
         'alpha' = OneDef [['Callithrix_kuhlii']] $ mean,
         'V' = post.vcv $ ss.grand.mean,
         'sigma_geoffroyi' = post.vcv $ ss.grand.mean,
         'sigma_kuhlii' = post.vcv $ ss.grand.mean,
         'sigma_jacchus' = post.vcv $ ss.grand.mean,
         'sigma_penicillata' = post.vcv $ ss.grand.mean,
         'sigma_argentata' = post.vcv $ ss.grand.mean,
         'sigma_humeralifera' = post.vcv $ ss.grand.mean,
         'sigma_pygmaea' = post.vcv $ ss.grand.mean
         )
       )

Callithrix.fit <- stan(file = '../Stan/test.stan',
                       data = Callithrix $ stan.data,
                       pars = c('B', 'alpha', 'V', 'sigma_geoffroyi',
                         'sigma_kuhlii', 'sigma_jacchus', 'sigma_penicillata', 
                         'sigma_argentata', 'sigma_humeralifera', 'sigma_pygmaea'),
                       init = Callithrix $ stan.start,
                       iter = 1000, chains = 1)

Callithrix $ post <- extract(Callithrix.fit)

par (mfrow = c(1, 2))
color2D.matplot (cov2cor (Callithrix $ post $ sigma_kuhlii [1, ,]))
color2D.matplot (cov2cor (Callithrix $ post $ sigma_kuhlii [500, ,]))

X11()
color2D.matplot (cov2cor (OneDef [['Callithrix_kuhlii']] $ ml.vcv))
