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
require (geiger)
require (mvtnorm)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

runcie.aux <- list()
runcie.aux $ subtree <- extract.clade (Tree [[5]], 138)

runcie.data <- list ()
runcie.data <-
  within (runcie.data,
          {
            k <- 39
            m <- 12 
            C <- solve(vcvPhylo(runcie.aux $ subtree))
            ni <- laply (OneDef [runcie.aux $ subtree $ tip.label],
                         function (L) nrow (L $ local))
            ni_max <- max(ni)
            X <- array (0, c(m, ni_max, k))
            for (i in 1:m)
              X[i, 1:ni[i], ] <- OneDef[[i+26]] $ local
            n_fac <- 6
            ### runcie decomposition
            a1W <- 3; b1W <- 2 # shrinkageW
            a2W <- 4; b2W <- 2
            asW <- 3; bsW <- 2
            a1S <- 3; b1S <- 2 # shrinkageB
            a2S <- 4; b2S <- 2
            asS <- 3; bsS <- 2
            niW <- max(ni); niS <- 11
          })
  
fail.runcie <-
  failwith(NULL, function (i)
           stan('../Stan/oneSigma_Anc_runcie_WBS.stan', data = runcie.data,
                pars = c('terminal', 'ancestor', 'root',
                  'LambdaW', 'LambdaS', 'PsiW', 'PsiS','SigmaW', 'SigmaB', 'SigmaS',
                  'deltaW', 'deltaS', 'phiW', 'phiS', 'variance_brownian'),
              warmup = 1000, iter = 2000, thin = 10, chains = 1,
                control = list ('chain_id' = i)))

runcie.test <- alply (1:4, 1, fail.runcie, .parallel = TRUE)

save (runcie.test, file = 'testRuncieWBS.RData')

rm (list = ls())

