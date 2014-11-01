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
attach ('../../covTensor/Work/post.vcv.RData')

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
            a1W <- 10; b1W <- 2 # shrinkageW
            a2W <- 10; b2W <- 5
            asW <- 10; bsW <- 2
            niW <- max(ni)
          })
  
fail.model <-
  failwith(NULL, function (i)
           stan('../Stan/oneSigma_Anc_runcie_drift.stan', data = data,
                pars = c('terminal', 'ancestor', 'root',
                  'LambdaW', 'PsiW','SigmaW', 'variance_brownian',
                  'deltaW', 'phiW'),
                warmup = 1000, iter = 2000, thin = 10, chains = 1,
                control = list ('chain_id' = i)))

fit <- alply (1:4, 1, fail.model, .parallel = TRUE)

save (data, fit, file = 'driftFlat.RData')

rm (list = ls())


alpha = 10
beta = 2
alpha2 = 10
beta2 = 5
boxplot (cbind (
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta))), 
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(1, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(2, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(3, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(4, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(5, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(6, shape = alpha2, rate = beta2))),
  1/raply (1000, prod (rgamma(1, shape = alpha, rate = beta),
                      rgamma(7, shape = alpha2, rate = beta2)))))

  
