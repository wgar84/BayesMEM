require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 12)
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

load ('pcaTest.RData')

pcaTest.ext <- llply (pcaTest, function (L) extract (L $ model))

str(pcaTest.ext [[1]])
     
pcaTest.ext <-
  llply (pcaTest.ext,
         function (M)
         llply (M, function (L)
                {
                  change.k <- which (dim (L) == 11)
                  if (length (change.k) > 0)
                    for (i in change.k)
                      {
                        dimnames (L) [[i]] <- paste ('PC', 1:11, sep = '')
                        names (dimnames (L)) [i] <- 'trait'
                      }
                  change.m <- which (dim (L) == 12)
                  if (length (change.m) > 0)
                    for (i in change.m)
                      {
                        dimnames (L) [[i]] <- extract.clade (Tree [[1]], 138) $ tip.label
                        names (dimnames (L)) [i] <- 'node'
                      }
                  change.m <- which (dim (L) == 10)
                  if (length (change.m) > 0)
                    for (i in change.m)
                      {
                        dimnames (L) [[i]] <- as.character (14:23)
                        names (dimnames (L)) [i] <- 'node'
                      }
                  L
                }))
              
pcaTest.gelman <- DiagGelman(pcaTest.ext)

pcaTest.gelman.comp <-
  do.call (arrangeGrob, c(pcaTest.gelman, 'ncol' = 3))

### WORKS!
