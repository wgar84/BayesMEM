require (ape)
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
require (grid)
require (gridExtra)
require (gridBase)
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

attach ('TestPR.RData')

prior.test <- list()
prior.test $ ext <- llply (Test.Prior, function (L) L $ ext)
#prior.test $ ext <- prior.test $ ext [-7]

prior.test $ subtree <- extract.clade(Tree [[5]], 138)

### QUANTIS
prior.test $ diag.quant <- llply (prior.test $ ext, function (L)
                                DiagQuantilePop(138, L, OneDef, Tree [[5]],
                                                at = c(0, 0.25, 0.75, 1)),
                                .parallel = TRUE)

pdf (file = 'prior.diag.quant.pdf', width = 24, height = 12)
for (i in 1:8)
  print(prior.test $ diag.quant [[i]])
dev.off (dev.cur ())

### TRACE

prior.test $ combine <-
  extract (sflist2stanfit (llply (Test.Prior, function (L) L $ model)))

prior.test $ combine <-
  llply (prior.test $ combine,
         function (L)
         {
           change.k <- which (dim (L) == 39)
           if (length (change.k) > 0)
             for (i in change.k)
               {
                 dimnames (L) [[i]] <- dimnames (prior.test $ ext [[1]] [[1]]) $ trait
                 names (dimnames (L)) [i] <- 'trait'
               }
           change.m <- which (dim (L) == 22)
           if (length (change.m) > 0)
             for (i in change.m)
               {
                 dimnames (L) [[i]] <- dimnames (prior.test $ ext [[1]] [[1]]) $ node
                 names (dimnames (L)) [i] <- 'node'
               }
           L
         })
           
prior.test $ diag.trace <- DiagTrace(prior.test $ combine, OneDef [[30]] $ mean)
           
ggsave('prior.diag.trace.jpg',
       do.call (arrangeGrob, c(prior.test $ diag.trace, 'ncol' = 3)),
       width = 36, height = 12)

### BETAS
prior.test $ post.response <-
  llply (prior.test $ ext, PostDeltaZ, tree = prior.test $ subtree, beta = TRUE)

prior.test $ combine.response <-
  PostDeltaZ (prior.test $ combine, tree = prior.test $ subtree, beta = TRUE)

pdf (file = 'prior.diag.beta.pdf', width = 14, height = 7)
for (i in 1:8)
  print (DiagBeta(prior.test $ post.response [[i]] $ Beta, prior.test $ subtree))
dev.off (dev.cur ())

### GELMAN
prior.test $ convergence <- DiagGelman(prior.test $ ext)

prior.test $ convergence.composite <-
  do.call (arrangeGrob, c(prior.test $ convergence, 'ncol' = 3))

ggsave(filename = 'prior.gelmanDiag.jpg',
       plot = prior.test $ convergence.composite,
       width = 21, height = 7)

save (prior.test, file = 'priorsDiag.RData')
