require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 3)
#registerDoMC (cores = 12)
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
#attach ('../../covTensor/Work/post.vcv.RData')

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

load ('pcaTest.RData')

pcaPost <- list ()

pcaPost $ subtree <- extract.clade (Tree [[1]], 138)
pcaPost $ ext <- llply (pcaTest, function (L) extract (L $ model))

pcaPost $ ext <- llply (pcaPost $ ext, nameStanExt,
                        names.list =
                        list ('trait' = paste ('PC', 1:11, sep = ''),
                              'node' = pcaPost $ subtree $ tip.label,
                              'node' = as.character (14:23)))

pcaPost $ gelman.plots <- DiagGelman(pcaPost $ ext)

pcaPost $ gelman.plot.comp <-
  do.call (arrangeGrob, c(pcaPost$gelman.plots, 'ncol' = 3))
ggsave('pca.gelman.pdf', pcaPost $ gelman.plot.comp, width = 24, height = 8)

pcaPost $ quantile <- llply (pcaPost $ ext, function (L)
                             DiagQuantilePop(138, L,
                                             or.pops = pcaTest [[1]] $ aux $ raw.rotated, 
                                             tree = Tree [[1]],
                                             at = c(0, 0.25, 0.75, 1)),
                             .parallel = TRUE)

pdf (file = 'pca.quantile.pdf', width = 24, height = 12)
for (i in 1:12)
  print(pcaPost $ quantile [[i]])
dev.off (dev.cur ())

pcaPost $ combine <- sflist2stanfit(llply (pcaTest, function (L) L $ model))

pcaPost $ combine.ext <- extract(pcaPost $ combine)

pcaPost $ combine.ext <-
  nameStanExt (pcaPost $ combine.ext, 
               names.list =
               list ('trait' = paste ('PC', 1:11, sep = ''),
                     'node' = pcaPost $ subtree $ tip.label,
                     'node' = as.character (14:23)))


pcaPost $ trace <- DiagTrace (pcaPost $ combine.ext, pcaTest [[1]] $ aux $ grand.mean)
pcaPost $ trace.comp <-
  do.call (arrangeGrob, c(pcaPost $ trace, 'ncol' = 3))
ggsave('pca.trace.pdf', pcaPost $ trace.comp, width = 24, height = 8)

###
pcaPost $ combine.response <- PostDeltaZ(pcaPost $ combine.ext, pcaPost $ subtree, TRUE)

DiagBeta(pcaPost $ combine.response, pcaPost $ subtree, 'PC5')

pcaPost $ response <- llply (pcaPost $ ext, PostDeltaZ,
                             tree = pcaPost $ subtree, beta = TRUE)

pdf (file = 'pca.beta.pdf', width = 16, height = 8)
for(i in 1:12) 
  print (DiagBeta (pcaPost $ response [[i]] $ Beta, pcaPost $ subtree, 'PC1'))
dev.off (dev.cur())

pcaPost $ ext.rotated <-
  llply (pcaPost $ ext, function (L)
         {
           llply (L, function (M)
                  {
                    rotation <- t(pcaTest [[1]] $ aux $ eigenW $ vectors [, 1:11])
                    to.use <- which(dim (M) != 11)
                    if (sum (dim (M) == 11) == 1)
                      M <- aaply (M, 1, function (X) X %*% rotation)
                    else if (sum (dim (M) == 11) == 2)
                      M <- aaply (M, 1, function (X)
                             ExtendMatrix (t (rotation) %*% X %*% rotation, cut.off = 11))
                    M
                  })
         })

pcaPost $ ext.rotated <-
  llply (pcaPost $ ext.rotated, nameStanExt,
         names.list =
         list ('trait' = names (OneDef [[1]] $ mean),
               'node' = pcaPost $ subtree $ tip.label,
               'node' = as.character (14:23)))

pcaPost $ rot.response <- llply (pcaPost $ ext.rotated, PostDeltaZ,
                             tree = pcaPost $ subtree, beta = TRUE)

pdf (file = 'pca.beta.pdf', width = 16, height = 8)
for(i in 1:12) 
  print (DiagBeta (pcaPost $ rot.response [[i]] $ Beta, pcaPost $ subtree, 'logCS'))
dev.off (dev.cur())

save (pcaPost, file = 'pcaPost.RData')
