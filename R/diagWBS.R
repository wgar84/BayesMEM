require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 4)
#registerDoMC (cores = 8)
#registerDoMC (cores = 60)
require (reshape2)
require (ggplot2)
require (gridBase)
require (gridExtra)
require (grid)
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

load ('testRuncieWBS.RData')

folder <- 'RuncieWBS/'

runciePost <- list()
runciePost $ ext <- llply (runcie.test, extract)

runciePost $ subtree <- extract.clade (Tree [[1]], 138)

runciePost $ ext <-
  llply (runciePost $ ext, nameStanExt,
         names.list = list ('trait' = colnames (OneDef [[1]] $ local),
           'node' = runciePost $ subtree $ tip.label,
           'node' = as.character (14:23),
           'factor' = paste ('lambda', 1:6, sep = ''),
           'iterations' = 1:100))

runciePost $ combine <- sflist2stanfit(runcie.test)

runciePost $ combine.ext <- extract(runciePost $ combine)

runciePost $ combine.ext <-
  nameStanExt (runciePost $ combine.ext, 
               names.list = list ('trait' = colnames (OneDef [[1]] $ local),
                 'node' = runciePost $ subtree $ tip.label,
                 'node' = as.character (14:23),
                 'factor' = paste ('lambda', 1:6, sep = ''),
                 'iterations' = 1:400))


runciePost $ gelman.plot <- DiagGelman(runciePost $ ext)

runciePost $ gelman.plot.comp <-
  do.call (arrangeGrob, c(runciePost$gelman.plot, 'nrow' = 3))
ggsave(paste (folder, 'gelman.pdf', sep = ''),
       runciePost $ gelman.plot.comp, width = 24, height = 8)
              
runciePost $ trace <-
  DiagTrace (runciePost $ combine.ext)
runciePost $ trace.comp <-
  do.call (arrangeGrob, c(runciePost $ trace, 'nrow' = 3))
ggsave(paste (folder, 'trace.pdf', sep = ''),
       runciePost $ trace.comp, width = 24, height = 24)

runciePost $ response <- llply (runciePost $ ext, PostDeltaZ,
                                tree = runciePost $ subtree,
                                beta = TRUE)

runciePost $ combine.response <-
  PostDeltaZ(runciePost $ combine.ext, runciePost $ subtree, TRUE)

pdf(file = paste (folder, 'beta.pdf', sep = ''), width = 16, height = 8)
llply(runciePost $ response, DiagBeta, runciePost $ subtree, slice = 'logCS')
dev.off (dev.cur ())

pdf(file = paste (folder, 'combine.beta.pdf', sep = ''), width = 16, height = 8)
DiagBeta(runciePost $ combine.response, runciePost $ subtree, slice = 'logCS')
dev.off (dev.cur ())

runciePost $ diagW <-
  llply (runciePost $ ext, DiagW,
         vcv.terminal.list = llply (OneDef[27:38],
           function (L) L $ ml.vcv),
         sample.size = Aux $ sample.size [27:38], 
         parallel = FALSE, .parallel = TRUE)

runciePost $ diagW.comp <- do.call(arrangeGrob, c(runciePost $ diagW, ncol = 2, nrow = 2))
ggsave(paste (folder, 'W.pdf', sep = ''),
       runciePost $ diagW.comp, width = 16, height = 16)

runciePost $ quantile <-
  foreach(i = 1:4) %dopar%
DiagQuantilePop(138, runciePost $ ext [[i]], OneDef, tree = Tree [[5]])
runciePost $ quantile.comp <-
  do.call(arrangeGrob, c(runciePost $ quantile, nrow = 4))
ggsave(paste (folder, 'quantile.pdf', sep = ''),
       runciePost $ quantile.comp, width = 16, height = 32)

runciePost $ Drift <- DiagDrift (runciePost $ combine.ext)
runciePost $ Drift.comp <-
  do.call(arrangeGrob, c(runciePost $ Drift, 'nrow' = 3))
ggsave(paste (folder, 'drift.pdf', sep = ''),
       runciePost $ Drift.comp, width = 12, height = 36)

runciePost $ Drift.mat <- DiagDrift (runciePost $ combine.ext, mode = 'matrix')
runciePost $ Drift.mat.comp <-
  do.call(arrangeGrob, c(runciePost $ Drift.mat, 'nrow' = 3))
ggsave(paste (folder, 'drift.mat.pdf', sep = ''),
       runciePost $ Drift.mat.comp, width = 12, height = 36)

### vamos olhar para os lambdas

runciePost $ LambdaW.df <- melt (llply (runciePost $ ext, function (L) L $ LambdaW))

runciePost $ LambdaW.plot <- 
  ggplot (runciePost $ LambdaW.df) +
  geom_violin(aes (x = trait, y = value, color = L1),
              alpha = 0.5, scale = 'width') +
  facet_wrap(~ Var2, scales = 'free', ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste (folder, 'LambdaW.pdf', sep = ''),
       runciePost $ LambdaW.plot, width = 24, height = 48)

runciePost $ LambdaS.df <- melt (llply (runciePost $ ext, function (L) L $ LambdaS))

runciePost $ LambdaS.plot <- 
  ggplot (runciePost $ LambdaS.df) +
  geom_violin(aes (x = trait, y = value, color = L1),
              alpha = 0.5, scale = 'width') +
  facet_wrap(~ Var2, scales = 'free', ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste (folder, 'LambdaS.pdf', sep = ''),
       runciePost $ LambdaS.plot, width = 24, height = 48)

save(runciePost, file = paste (folder, 'post.RData', sep = ''))

### escrever DiagEvol

Norm.hyp <- t (aaply (Aux $ def.hyp [, 1:6], 2, Normalize))

runciePost $ combine.evol.W <- 
  aaply (runciePost $ combine.ext $ SigmaW, 1,
         function (C) aaply (Norm.hyp, 2, Evolvability, cov.matrix = C))

runciePost $ combine.evol.S <- 
  aaply (runciePost $ combine.ext $ SigmaS, 1,
         function (C) aaply (Norm.hyp, 2, Evolvability, cov.matrix = C))

runciePost $ evol.df <-
  cbind (melt(runciePost $ combine.evol.S, value.name = 'S'),
         melt(runciePost $ combine.evol.W, value.name = 'W'))

runciePost $ evol.df <- runciePost $ evol.df [, -(1:2)]

ggplot (runciePost $ evol.df) +
  geom_point (aes (x = S, y = W, color = Var2), alpha = 1, size = 5, shape = '+') +
  scale_x_log10() + scale_y_log10() +
  theme_minimal() + geom_smooth(aes(x = S, y = W), formula = y ~ x, method = 'lm')

runciePost $ evol.df $ ratio <- runciePost $ evol.df $ B / runciePost $ evol.df $ W

ggplot(runciePost $ evol.df) +
  geom_violin(aes (y = log(B/W), x = Var2))
