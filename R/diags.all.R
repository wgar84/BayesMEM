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

boxplot (t (apply (AllPost $ ext $ Sigma_bm, 1, function (x) eigen (x) $ values)))

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

AllPost $ response <- PostDeltaZ(AllPost $ ext, Tree [[1]], TRUE)

names (dimnames (AllPost $ response $ Beta)) <- c('iter', 'edge', 'trait')

dimnames (AllPost $ response $ Beta) [[3]] <- colnames (OneDef [[1]] $ local)
dimnames (AllPost $ response $ Beta) [[2]] <-
  paste (c(Tree [[1]] $ tip.label, as.character (110:217)) [Tree [[1]] $ edge [, 1]],
         c(Tree [[1]] $ tip.label, as.character (110:217)) [Tree [[1]] $ edge [, 2]],
         sep = '-')

names (dimnames (AllPost $ response $ DeltaZ)) <-
  names (dimnames (AllPost $ response $ Beta))

dimnames (AllPost $ response $ DeltaZ) <- dimnames (AllPost $ response $ Beta)

par (mar = c(10, 5, 5, 5))
boxplot (AllPost $ response $ Beta [, , 'logCS'], las = 3, cex.axis = 0.7,
         main = 'Beta_i (logCS)', cex = 0.5, cex.axis = 0.5)

abline (h = 0, lty = 3, col = 'red')

abs (colMeans (AllPost $ response $ Beta [, , 'logCS'])) /
     max(abs (colMeans (AllPost $ response $ Beta [, , 'logCS'])))



AllPost $ mean.beta.logcs <- colMeans (AllPost $ response $ Beta [, , 'logCS'])

ifelse (AllPost $ mean.beta.logcs > 0, rgb (1, 0, 0, 0.2), rgb (0, 0, 1, 0.2))

pdf ('tree.beta.logcs.pdf', width = 15, height = 15)
plot (Tree [[1]], cex = 0.5, use.edge.length = FALSE)
edgelabels (pch = 20,
            cex = 10 * abs (AllPost $ mean.beta.logcs) / max(AllPost $ mean.beta.logcs),
            col = ifelse (AllPost $ mean.beta.logcs > 0,
              rgb (1, 0, 0, 0.2), rgb (0, 0, 1, 0.2)))


edgelabels (pch = c('', '*') [aaply (apply (AllPost $ response $ Beta [, , 'logCS'], 2, range), 2, function (x) prod (x) > 0) + 1])

dev.off (dev.cur ())

aaply (apply (AllPost $ response $ Beta [, , 'logCS'], 2, range), 2,
       function (x) prod (x) > 0)

AllPost $ beta.df <- melt (AllPost $ response $ Beta)

pdf ('all.betas.pdf', width = 100, height = 100)
ggplot(AllPost $ beta.df, aes (x = trait, y = value)) +
  geom_boxplot (outlier.shape = NA) + facet_wrap(~ edge, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1))
dev.off (dev.cur ())



AllPost $ response $ Beta [, 1, ] %*%
  rbind (rep (0, 8), Aux $ ed.hyp [[1]] [-20, ])

fo
boxplot (t (apply (AllPost $ response $ Beta [, i, ], 1, Normalize)) %*%
  apply (rbind (rep (0, 8), Aux $ ed.hyp [[1]] [-20, ]), 2, Normalize))

boxplot (aaply (AllPost $ response $ DeltaZ, c(1, 2), Norm), las = 3, cex.axis = 0.3)

boxplot (aaply (AllPost $ ext $ Xbar [, 110:216, ], c(2, 3), mean))

points (1:39, OneDef [['Cebus_apella']] $ mean, col = 'red')
