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

### com priors

attach ('Test.prior.RData')

test.prior.post <- list()
test.prior.post $ ext <- extract(Test.prior)

test.prior.post $ quantile <-  foreach (i = 1:12) %dopar%
DiagQuantile (OneDef [[i+26]] $ local,
              test.prior.post $ ext $ Xbar [, i, ],
              test.prior.post $ ext $ Sigma, at = seq(0, 1, 0.25))

pdf (file = 'diags_test_prior.pdf', width = 15, height = 12)
for (i in 1:12)
  plot (test.prior.post $ quantile [[i]] +
        labs (title = names (OneDef) [i+26]) + theme_minimal() +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())

test.prior.post $ subtree <- extract.clade(Tree [[5]], 138)

test.prior.post $ response <- PostDeltaZ(test.prior.post $ ext,
                                         test.prior.post $ subtree, TRUE)

names (dimnames (test.prior.post $ response $ Beta)) <- c('iter', 'edge', 'trait')

dimnames (test.prior.post $ response $ Beta) [[3]] <- colnames (OneDef [[1]] $ local)
dimnames (test.prior.post $ response $ Beta) [[2]] <-
  paste (c(test.prior.post $ subtree $ tip.label,
           as.character (13:23)) [test.prior.post $ subtree $ edge [, 1]],
         c(test.prior.post $ subtree $ tip.label,
           as.character (13:23))  [test.prior.post $ subtree $ edge [, 2]],
         sep = '-')

names (dimnames (test.prior.post $ response $ DeltaZ)) <-
  names (dimnames (test.prior.post $ response $ Beta))

dimnames (test.prior.post $ response $ DeltaZ) <-
  dimnames (test.prior.post $ response $ Beta)

par (mar = c(10, 5, 5, 5))
boxplot (test.prior.post $ response $ Beta [, , 'logCS'], las = 3, cex.axis = 0.7,
         main = 'Beta_i (logCS)', cex = 0.5, cex.axis = 0.5)

boxplot (aaply (test.prior.post $ response $ DeltaZ, c(1, 2), Norm),
         las = 3, cex.axis = 0.3)

boxplot (aaply (test.prior.post $ ext $ Xbar [, 13:21, ], c(2, 3), mean))

points (1:39, OneDef [['Pithecia_pithecia']] $ mean, col = 'red')

test.prior.post $ beta.df <- melt (test.prior.post $ response $ Beta)

pdf ('test.prior.betas.pdf', width = 12, height = 12)
ggplot(test.prior.post $ beta.df, aes (x = trait, y = value)) +
  geom_boxplot (outlier.shape = NA) + facet_wrap(~ edge, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1))
dev.off (dev.cur ())

test.prior.post $ mean.beta.logcs <-
  colMeans (test.prior.post $ response $ Beta [, , 'logCS'])

pdf ('test.prior.beta.logcs.pdf', width = 15, height = 15)
plot (test.prior.post $ subtree, cex = 0.5, use.edge.length = FALSE)
edgelabels (pch = 20,
            cex = 10 * abs (test.prior.post $ mean.beta.logcs) /
            max(test.prior.post $ mean.beta.logcs),
            col = ifelse (test.prior.post $ mean.beta.logcs > 0,
              rgb (1, 0, 0, 0.2), rgb (0, 0, 1, 0.2)))

edgelabels (pch = c('', '*') [aaply (
              apply (test.prior.post $ response $ Beta [, , 'logCS'],
                     2, range), 2, function (x) prod (x) > 0) + 1])

### sem priors

attach ('Test.RData')

test.post <- list()
test.post $ ext <- extract(Test)

test.post $ quantile <-  foreach (i = 1:12) %dopar%
DiagQuantile (OneDef [[i+26]] $ local,
              test.post $ ext $ Xbar [, i, ],
              test.post $ ext $ Sigma, at = seq(0, 1, 0.25))

pdf (file = 'diags_test.pdf', width = 15, height = 12)
for (i in 1:12)
  plot (test.post $ quantile [[i]] +
        labs (title = names (OneDef) [i+26]) + theme_minimal() +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())

test.post $ subtree <- extract.clade(Tree [[5]], 138)

test.post $ response <- PostDeltaZ(test.post $ ext,
                                         test.post $ subtree, TRUE)

names (dimnames (test.post $ response $ Beta)) <- c('iter', 'edge', 'trait')

dimnames (test.post $ response $ Beta) [[3]] <- colnames (OneDef [[1]] $ local)
dimnames (test.post $ response $ Beta) [[2]] <-
  paste (c(test.post $ subtree $ tip.label,
           as.character (13:23)) [test.post $ subtree $ edge [, 1]],
         c(test.post $ subtree $ tip.label,
           as.character (13:23))  [test.post $ subtree $ edge [, 2]],
         sep = '-')

names (dimnames (test.post $ response $ DeltaZ)) <-
  names (dimnames (test.post $ response $ Beta))

dimnames (test.post $ response $ DeltaZ) <-
  dimnames (test.post $ response $ Beta)

par (mar = c(10, 5, 5, 5))
boxplot (test.post $ response $ Beta [, , 'logCS'], las = 3, cex.axis = 0.7,
         main = 'Beta_i (logCS)', cex = 0.5, cex.axis = 0.5)

boxplot (aaply (test.post $ response $ DeltaZ, c(1, 2), Norm),
         las = 3, cex.axis = 0.3)

boxplot (aaply (test.post $ ext $ Xbar [, 13:21, ], c(2, 3), mean))

points (1:39, OneDef [['Pithecia_pithecia']] $ mean, col = 'red')

test.post $ beta.df <- melt (test.post $ response $ Beta)

pdf ('test.betas.pdf', width = 12, height = 12)
ggplot(test.post $ beta.df, aes (x = trait, y = value)) +
  geom_boxplot (outlier.shape = NA) + facet_wrap(~ edge, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1))
dev.off (dev.cur ())

test.post $ mean.beta.logcs <- colMeans (test.post $ response $ Beta [, , 'logCS'])

pdf ('test.beta.logcs.pdf', width = 15, height = 15)
plot (test.post $ subtree, cex = 0.5, use.edge.length = FALSE)
edgelabels (pch = 20,
            cex = 10 * abs (test.post $ mean.beta.logcs) /
            max(test.post $ mean.beta.logcs),
            col = ifelse (test.post $ mean.beta.logcs > 0,
              rgb (1, 0, 0, 0.2), rgb (0, 0, 1, 0.2)))

edgelabels (pch = c('', '*') [aaply (apply (test.post $ response $ Beta [, , 'logCS'],
              2, range), 2, function (x) prod (x) > 0) + 1])

### diagnostico de convergencia

dimnames (test.prior.post $ ext $ alpha) <-
  list (c(1:100), colnames (OneDef [[1]] $ local))

names (dimnames (test.prior.post $ ext $ alpha)) <- c('iteration', 'trait')

test.prior.post $ alpha.df <- melt (test.prior.post $ ext $ alpha)
test.prior.post $ alpha.df [, 'iteration'] <-
  as.numeric (test.prior.post $ alpha.df [, 'iteration'])

test.prior.post $ alpha.prior.df <-
  melt (OneDef [['Pithecia_pithecia']] $ mean)

test.prior.post $ alpha.prior.df $ trait <- rownames (test.prior.post $ alpha.prior.df)

ggplot (test.prior.post $ alpha.df, aes (x = iteration, y = value)) +
  geom_point() + geom_line() +
  geom_hline(data = test.prior.post $ alpha.prior.df, aes (yintercept = value)) + 
  facet_wrap(~ trait)

dim (test.prior.post $ ext $ Xbar)

plot (test.prior.post $ subtree, direction = 'upwards', cex = 0.7)
nodelabels()

dimnames (test.prior.post $ ext $ Xbar) <-
  list (c(1:100),
        c(test.prior.post $ subtree $ tip.label, as.character (14:23)),
        colnames (OneDef [[1]] $ local))

names (dimnames (test.prior.post $ ext $ Xbar)) <- c ('iteration', 'node', 'trait')

test.prior.post $ xbar.df <- melt (test.prior.post $ ext $ Xbar)

test.prior.post $ xbar.df [, 'iteration'] <-
  as.numeric (test.prior.post $ xbar.df [, 'iteration'])

X11()
ggplot (subset (test.prior.post $ xbar.df,
                node == '21'),
        aes (x = iteration, y = value)) +
  geom_point() + geom_line() +
  geom_hline(data = test.prior.post $ alpha.prior.df, aes (yintercept = value)) + 
  facet_wrap(~ trait)

dimnames (test.post $ ext $ alpha) <-
  list (c(1:100), colnames (OneDef [[1]] $ local))

names (dimnames (test.post $ ext $ alpha)) <- c('iteration', 'trait')

test.post $ alpha.df <- melt (test.post $ ext $ alpha)
test.post $ alpha.df [, 'iteration'] <-
  as.numeric (test.post $ alpha.df [, 'iteration'])

test.post $ alpha.start.df <-
  melt (OneDef [['Pithecia_pithecia']] $ mean)

test.post $ alpha.start.df $ trait <- rownames (test.post $ alpha.start.df)

ggplot (test.post $ alpha.df, aes (x = iteration, y = value)) +
  geom_point() + geom_line() +
  geom_hline(data = test.post $ alpha.start.df, aes (yintercept = value)) + 
  facet_wrap(~ trait)

dim (test.post $ ext $ Xbar)

plot (test.post $ subtree, direction = 'upwards', cex = 0.7)
nodelabels()

dimnames (test.post $ ext $ Xbar) <-
  list (c(1:100),
        c(test.post $ subtree $ tip.label, as.character (14:23)),
        colnames (OneDef [[1]] $ local))

names (dimnames (test.post $ ext $ Xbar)) <- c ('iteration', 'node', 'trait')

test.post $ xbar.df <- melt (test.post $ ext $ Xbar)

test.post $ xbar.df [, 'iteration'] <-
  as.numeric (test.post $ xbar.df [, 'iteration'])

X11()
ggplot (subset (test.post $ xbar.df,
                node == '21'),
        aes (x = iteration, y = value)) +
  geom_point() + geom_line() +
  geom_hline(data = test.post $ alpha.start.df, aes (yintercept = value)) + 
  facet_wrap(~ trait)
