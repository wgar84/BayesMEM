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

AllPost $ ext <- extract (All $ fit.SA)

dim (AllPost $ ext $ Sigma)

color2D.matplot (cov2cor (AllPost $ ext $ Sigma_bm [10, , ]))

AllPost $ quantile <-  foreach (i = 1:109) %dopar%
DiagQuantile (OneDef [[i]] $ local,
              AllPost $ ext $ Xbar [, i, ],
              AllPost $ ext $ Sigma, at = seq(0, 1, 0.25))

pdf (file = 'diags_all_1Sig.pdf', width = 15, height = 12)
for (i in 1:109)
  plot (AllPost $ quantile [[i]] +
        labs (title = names (OneDef) [i]) +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())

AllPost $ ml.means <- laply (OneDef, function (L) L $ mean)

dimnames (AllPost $ ml.means) <- list (names (OneDef), colnames (OneDef [[1]] $ local))

dim (AllPost $ ml.means)

AllPost $ quant.means <- DiagQuantile (AllPost $ ml.means,
                                       AllPost $ ext $ alpha,
                                       AllPost $ ext $ Sigma_bm, at = seq(0, 1, 0.25))

pdf (file = 'diags_mean_1Sig.pdf', width = 15, height = 12)
plot (AllPost $ quant.means) + theme (axis.text = element_text (size=6))
dev.off (dev.cur ())
