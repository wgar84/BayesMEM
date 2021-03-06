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

node <- 138
subtree <- extract.clade (Tree [[5]], node)

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
            n_fac <- 3
            ###  decomposition
            a1W <- 20; b1W <- 10 # shrinkageW
            a2W <- 20; b2W <- 15
            asW <- 2; bsW <- 1
            niW <- min(ni)
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

data $ subtree <- subtree
data $ node <- node

save (data, fit, file = 'driftFlat.RData')

rm (list = ls())
