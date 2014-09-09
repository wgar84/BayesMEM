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
attach ('../../covTensor/Work/etd.all.RData')


source ('../../covTensor/Func/MeanMatrix.R')
source ('../../covTensor/Func/Pairwise.R')
source ('../../covTensor/Func/Frobenius.R')
source ('../../covTensor/Func/MatrixDist.R')
source ('../../covTensor/Func/BuildSigma.R')
source ('../../covTensor/Func/EigenTensorDecomposition.R')
source ('../../covTensor/Func/RebuildETD.R')
source ('../FuncR/DiagQuantile.R')
source ('../../Func/plot.R')

source ('../../Func/splines.R')


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


dim (CallPost $ response $ Beta)

CallPost $ beta.cred <-
  aaply (CallPost $ response $ Beta, c(2, 3), quantile, probs = c(0.025, 0.975))

CallPost $ beta.zero.out.df <-
  melt (aaply (CallPost $ beta.cred, c(1, 2), prod) > 0)

CallPost $ beta.zero.out.df $ value <-
  as.factor (c('', '*') [CallPost $ beta.zero.out.df $ value + 1])

par (mfrow = c(1, 2))

X11()
plot (CallTree)
nodelabels()

ggplot(CallPost $ beta.df, aes (x = edge, y = value)) +
  geom_boxplot (outlier.shape = NA) + facet_wrap(~ trait, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1))

### DRIFTSEL

par (mfrow = c(1, 2))

boxplot (aaply (CallPost $ ext $ Sigma, 1, function (x) eigen(x)$values), main = 'W')
abline (h = 0)

boxplot (aaply (CallPost $ ext $ Sigma_bm, 1, function (x) eigen(x)$values), main = 'B')
abline (h = 0)

CallPost $ driftsel <- list ()

CallPost $ driftsel $ evalW  <-
  aaply (CallPost $ ext $ Sigma, 1, function (x) eigen(x)$values)

CallPost $ driftsel $ evalB  <-
  aaply (CallPost $ ext $ Sigma_bm, 1, function (x) eigen(x)$values)

### cov (b, w) / var (w)

CallPost $ driftsel <-
  within (CallPost $ driftsel, 
          eval.slopes <- aaply (1:100, 1, function (i)
                                coef (lm (log (evalB[i, 1:6]) ~ log (evalW[i, 1:6]))) [2]))
                                

qplot (CallPost $ driftsel $ eval.slopes, geom = 'histogram')

par (mfrow = c(1, 2))

boxplot (aaply (CallPost $ ext $ Sigma [, -1, -1], 1,
                function (x) eigen(x)$values), main = 'W')
abline (h = 0)

boxplot (aaply (CallPost $ ext $ Sigma_bm [, -1, -1], 1,
                function (x) eigen(x)$values), main = 'B')
abline (h = 0)

## CallPost $ driftsel.ns <- list ()

## CallPost $ driftsel.ns $ evalW  <-
##   aaply (CallPost $ ext $ Sigma[, -1, -1], 1, function (x) eigen(x)$values)

## CallPost $ driftsel.ns $ evalB  <-
##   aaply (CallPost $ ext $ Sigma_bm[, -1, -1], 1, function (x) eigen(x)$values)

## ### cov (b, w) / var (w)

## CallPost $ driftsel.ns <-
##   within (CallPost $ driftsel.ns, 
##           eval.slopes <- aaply (1:100, 1, function (i)
##                                 coef (lm (log (evalB[i, 1:6]) ~ log (evalW[i, 1:6]))) [2]))
                                

## qplot (CallPost $ driftsel.ns $ eval.slopes, geom = 'histogram')

CallPost $ driftsel $ evecW <-
  aaply (CallPost $ ext $ Sigma, 1,
         function (x) eigen (x) $ vectors [, 1:6])

CallPost $ driftsel $ Gamma <-
  aaply (1:100, 1, function (i)
         {
           t (CallPost $ driftsel $ evecW [i, , ]) %*%
             CallPost $ ext $ Sigma_bm [i, , ] %*%
               CallPost $ driftsel $ evecW [i, , ]
         })

dim(CallPost $ driftsel $ Gamma)
