require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 4)
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
runcie.aux $ subtree <- extract.clade (Tree [[1]], 138)

runcie.data <- list ()
runcie.data <-
  within (runcie.data,
          {
            k <- 39
            m <- 12 
            C <- solve(cov2cor (vcvPhylo(runcie.aux $ subtree)))
            ni <- laply (OneDef [runcie.aux $ subtree $ tip.label],
                         function (L) nrow (L $ local))
            ni_max <- max(ni)
            X <- array (0, c(m, ni_max, k))
            for (i in 1:m)
              X[i, 1:ni[i], ] <- OneDef[[i+26]] $ local
            n_fac <- 8
            ### runcie decomposition
            a1W <- 5; b1W <- 2 # shrinkageW
            a2W <- 4; b2W <- 2
            asW <- 2; bsW <- 1
            a1B <- 5; b1B <- 2 # shrinkageB
            a2B <- 4; b2B <- 2
            asB <- 2; bsB <- 1
            niW <- max(ni); niB <- 11
          })
  
### precio fazer três desses:
##### um com matriz de covariância vinda da Tree[[5]] ### cov5
##### outra com a matriz de correlação vinda da Tree[[1]] ### cor1
##### e uma (ou duas) sem ancestrais


fail.runcie <-
  failwith(NULL, function (i)
           stan('../Stan/oneSigma_Anc_runcie.stan', data = runcie.data,
                pars = c('terminal', 'ancestor', 'root',
                  'LambdaW', 'LambdaB', 'PsiW', 'PsiB','SigmaW', 'SigmaB',
                  'deltaW', 'deltaB', 'phiW', 'phiB'),
              warmup = 1000, iter = 2000, thin = 10, chains = 1,
                control = list ('chain_id' = i)))

runcie.test <- alply (1:4, 1, fail.runcie, .parallel = TRUE)

save (runcie.test, file = 'testRuncieCor1.RData')

rm (list = ls())

