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

CallPost $ ext $ Sigma <- aperm (CallPost $ ext $ Sigma, c(2, 3, 1))

CallPost $ projSigma <- aaply (CallPost $ ext $ Sigma, 3, ProjectMatrix,
                               etd = etd.ss $ mean.etd,
                               mean.matrix = post.vcv $ ss.grand.mean, .parallel = TRUE)


CallPost $ projSigma_bm <- aaply (CallPost $ ext $ Sigma_bm, 1, ProjectMatrix,
                                  etd = etd.ss $ mean.etd,
                                  mean.matrix = post.vcv $ ss.grand.mean, .parallel = TRUE)

C.allSigma $ proj <-
  laply (C.allSigma $ post [4:10],
         function (L) aaply (L, 1, ProjectMatrix,
                             etd = etd.ss $ mean.etd,
                             mean.matrix = post.vcv $ ss.grand.mean), .parallel = TRUE)


Callithrix.who <- which (grepl ('Callithrix', names (OneDef)))

plot (etd.ss $ mean.etd $ projection [, 1:2], 
      pch = 20, cex = 1, xlim = c (-10, 12), ylim = c(-2.5, 2), xlab = 'PM1', ylab = 'PM2')

points (etd.ss $ mean.etd $ projection [Callithrix.who, 1:2], col = hsv (1:7/7, 1, 1, 1),
      pch = 20, cex = 2, xlim = c (-10, 12), ylim = c(-2.5, 2), xlab = 'PM1', ylab = 'PM2')


text(x = etd.ss $ mean.etd $ projection [Callithrix.who, 1:2],
     labels = names (OneDef) [Callithrix.who], pos = 2, cex = 0.5, )

for (i in Callithrix.who)
  {
    hull <- chull (t (etd.ss $ post.vcv.mean.proj [1:2, i, ]))
    polygon (t (etd.ss $ post.vcv.mean.proj [1:2, i, ]) [hull , ],
             col = hsv ((i-31)/7, 1, 1, 0.3), border = hsv ((i-31)/7, 1, 1, 0.5), lwd = 2)
  }

hull <- chull (CallPost $ projSigma [, 1:2])

polygon (CallPost $ projSigma [hull, 1:2],
         col = hsv (1, 0, 0.8, 0.3), border = hsv (1, 0, 0.8, 0.5), lwd = 2)

hull <- chull (t(etd.ss $ post.anc.proj [1:2, , '143']))

polygon (t(etd.ss $ post.anc.proj [1:2, hull, '143']),
         col = hsv (1, 0, 0.5, 0.1), border = hsv (1, 0, 0.5, 0.15), lwd = 2)


for (i in 1:7)
  {
    hull <- chull (C.allSigma $ proj[i, , 1:2])
    polygon (C.allSigma $ proj[i, hull, 1:2],
             col = hsv (i/7, 1, 1, 0.3), border = hsv (i/7, 1, 1, 0.5), lwd = 2)
  }

legend ('topleft', col = hsv ((1:7)/7, 1, 1, 0.5), pch = '', lty = 1, lwd = 3,
        legend = Tree [[1]] $ tip.label [Callithrix.who], bty = 'n', cex = 1)

X11()
plot (extract.clade(Tree [[1]], 143))

#### a posteriori pro modelo completo não bate com a posteriori dos modelos pra espécie
#### bom, isso até faz sentido se eu considerar a estrutura de cada modelo...
#### vou ter que esperar pra ver como ficam as posterioris associadas a cada modelo de W anc.
