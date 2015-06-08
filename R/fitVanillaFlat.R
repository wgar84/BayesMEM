require (ape)
require (evolqg)
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

attach('../../Tese/Data/allo.Results.RData')

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

data <- list ()
data <-
  within (data,
          {
            k <- 4
            m <- 109 
            C <- solve(vcvPhylo(Tree [[5]]))
            ni <- Aux $ sample.size
            ni_max <- max(ni)
            X <- array (0, c(m, ni_max, k))
            for (i in 1:m)
              X[i, 1:ni[i], ] <-
                as.matrix (subset (allo.Data $ onedef.df,
                                   animal == Tree [[1]] $ tip.label [i]) [, 3:6])
          })

data $ ni

fail.model <-
  failwith(NULL, function (i)
           stan('../Stan/oneSigma_Anc.stan', data = data,
                pars = c('terminal', 'ancestor', 'root',
                  'SigmaW', 'SigmaB'),
                warmup = 1000, iter = 2000, thin = 10, chains = 1,
                control = list ('chain_id' = i)))

fit <- alply (1:4, 1, fail.model, .parallel = TRUE)

save (data, fit, file = 'vanillaFlat.RData')


