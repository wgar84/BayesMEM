require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 4)
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
#attach ('../../covTensor/Work/post.vcv.RData')

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

subtree <- extract.clade (Tree [[5]], 138)

data <- list ()
data <-
  within (data,
          {
            k <- 39
            m <- 12 
            C <- solve(vcvPhylo(subtree))
            ni <- laply (OneDef [subtree $ tip.label],
                         function (L) nrow (L $ local))
            ni_max <- max(ni)
            X <- array (0, c(m, ni_max, k))
            for (i in 1:m)
              X[i, 1:ni[i], ] <- OneDef[[i+26]] $ local
            n_fac <- 4
            ###  decomposition
            a1W <- 2; b1W <- 1 # shrinkageW
            a2W <- 3; b2W <- 1
            asW <- 3; bsW <- 1
            a1S <- 2; b1S <- 1 # shrinkageS
            a2S <- 3; b2S <- 1
            asS <- 3; bsS <- 1
            niW <- max(ni); niS <- 11
          })
  
fail.model <-
  failwith(NULL, function (i)
           stan('../Stan/oneSigma_Anc_runcie_WBS.stan', data = data,
                pars = c('terminal', 'ancestor', 'root',
                  'LambdaW', 'LambdaS', 'PsiW', 'PsiS','SigmaW', 'SigmaB', 'SigmaS', 
                  'deltaW', 'deltaS', 'phiW', 'phiS', 'variance_brownian'),
                warmup = 1000, iter = 2000, thin = 10, chains = 1,
                control = list ('chain_id' = i)))

fit <- alply (1:4, 1, fail.model, .parallel = TRUE)

save (data, fit, file = 'vanillaFlat.RData')

rm (list = ls())

