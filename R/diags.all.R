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

attach ('All.RData')

AllPost <- list ()

AllPost $ ext <- extract (All $ fit.oneSigma)

dim (AllPost $ ext $ Sigma)

color2D.matplot (cov2cor (AllPost $ ext $ Sigma_bm [500, , ]))

AllPost $ quantile <-  foreach (i = 1:109) %dopar%
DiagQuantile (All $ oneSigma $ X [i, 1:All $ oneSigma $ ni[i], ],
              AllPost $ ext $ Xbar [, i, ],
              AllPost $ ext $ Sigma, at = seq(0, 1, 0.25))

pdf (file = 'diags_all_1Sig.pdf', width = 15, height = 12)
for (i in 1:109)
  plot (AllPost $ quantile [[i]] +
        labs (title = colnames (All $ oneSigma $ C) [i]) +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())
