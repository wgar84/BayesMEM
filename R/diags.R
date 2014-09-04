require (ape)
require (MCMCglmm)
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
require (mvtnorm)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')
attach  ('../../covTensor/Work/etd.all.RData')


source ('../../covTensor/Func/MeanMatrix.R')
source ('../../covTensor/Func/Pairwise.R')
source ('../../covTensor/Func/Frobenius.R')
source ('../../covTensor/Func/MatrixDist.R')
source ('../../covTensor/Func/BuildSigma.R')
source ('../../covTensor/Func/EigenTensorDecomposition.R')
source ('../../covTensor/Func/RebuildETD.R')
source ('../FuncR/DiagQuantile.R')

C.allSigma <- Callithrix
rm (list = ls(pat = 'Callithrix'))
attach('Callithrix.RData')

names (Callithrix)

CallPost <- list()
CallPost $ raw <- sflist2stanfit(Callithrix $ fit.parl)

CallPost $ ext <- extract(CallPost $ raw)

str (CallPost $ ext)

CallPost $ quantile <-
  foreach (i = 1:7) %dopar%
DiagQuantile (Callithrix $ ancSigma $ X [i, 1:Callithrix $ ancSigma $ ni[i], ],
              CallPost $ ext $ Xbar [, i, ],
              CallPost $ ext $ Sigma, at = seq(0, 1, 0.25))

pdf (file = 'diags_Call_1Sig.pdf', width = 15, height = 12)
for (i in 1:7)
  plot (CallPost $ quantile [[i]] +
        labs (title = colnames (Callithrix $ ancSigma $ invC) [i]) +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())

#### a posteriori pro modelo completo não bate com a posteriori dos modelos pra espécie
#### bom, isso até faz sentido se eu considerar a estrutura de cada modelo...
#### vou ter que esperar pra ver como ficam as posterioris associadas
#### a cada modelo de W anc.

CallTree <- extract.clade (Tree [[1]], 143)

paste (c(CallTree $ tip.label, as.character (8:13)) [CallTree $ edge [, 1]],
       c(CallTree $ tip.label, as.character (8:13)) [CallTree $ edge [, 2]], sep = '-')

CallPost $ response <- PostDeltaZ(CallPost $ ext, CallTree, TRUE)

names (dimnames (CallPost $ response $ Beta)) <- c('iter', 'edge', 'trait')

dimnames (CallPost $ response $ Beta) [[3]] <- colnames (OneDef [[1]] $ local)
dimnames (CallPost $ response $ Beta) [[2]] <-
  paste (c(CallTree $ tip.label, as.character (8:13)) [CallTree $ edge [, 1]],
         c(CallTree $ tip.label, as.character (8:13)) [CallTree $ edge [, 2]], sep = '-')

dimnames (CallPost $ response $ Beta) [[2]] <-
  gsub ('Callithrix_', '',
        dimnames (CallPost $ response $ Beta) [[2]])

names (dimnames (CallPost $ response $ DeltaZ)) <-
  names (dimnames (CallPost $ response $ Beta))

dimnames (CallPost $ response $ DeltaZ) <- dimnames (CallPost $ response $ Beta)

CallPost $ beta.df <- melt (CallPost $ response $ Beta)

head (CallPost $ beta.df)

pdf (file = 'beta_i_logCS.pdf', width = 10, height = 8)

layout (array (c(1, 2, 3, 3), c (2, 2)))

plot (CallTree)
nodelabels()

boxplot (OneDef [[32]] $ local [, 'logCS'], xlim = c (0, 8), ylim = c (4.2, 4.7), at = 1)
for (i in 33:38)
  boxplot(OneDef [[i]] $ local [, 'logCS'], at = i - 31, add = TRUE)
axis(side = 1, labels = names (OneDef) [32:38], las = 3, at = 1:7, cex.axis = 0.5)

boxplot (CallPost $ response $ Beta [, , 'logCS'], las = 3, cex.axis = 0.7,
         main = 'Beta_i (logCS)')

dev.off (dev.cur ())



