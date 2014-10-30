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

load ('testRuncieDrift.RData')

folder <- 'RuncieDrift/'

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
ggsave(paste (folder, 'runcie.gelman.pdf', sep = ''),
       runciePost $ gelman.plot.comp, width = 8, height = 24)
              

runciePost $ trace <-
  DiagTrace (runciePost $ combine.ext)
runciePost $ trace.comp <-
  do.call (arrangeGrob, c(runciePost $ trace, 'nrow' = 3))
ggsave(paste (folder, 'runcie.trace.pdf', sep = ''),
       runciePost $ trace.comp, width = 24, height = 24)

runciePost $ response <- llply (runciePost $ ext, PostDeltaZ,
                                tree = runciePost $ subtree,
                                beta = TRUE)

runciePost $ combine.response <-
  PostDeltaZ(runciePost $ combine.ext, runciePost $ subtree, TRUE)

pdf(file = paste (folder, 'runcie.beta.pdf', sep = ''), width = 16, height = 8)
llply(runciePost $ response, DiagBeta, runciePost $ subtree, slice = 'logCS')
dev.off (dev.cur ())

pdf(file = paste (folder, 'combine.beta.pdf', sep = ''), width = 16, height = 8)
DiagBeta(runciePost $ combine.response, runciePost $ subtree, 'logCS')
dev.off (dev.cur ())


runciePost $ diagW <-
  llply (runciePost $ ext, DiagW,
         vcv.terminal.list = llply (OneDef[27:38],
           function (L) L $ ml.vcv),
         sample.size = Aux $ sample.size [27:38], 
         parallel = FALSE, .parallel = TRUE)

runciePost $ diagW.comp <- do.call(arrangeGrob, c(runciePost $ diagW, ncol = 2, nrow = 2))
ggsave(paste (folder, 'runcie.W.pdf', sep = ''),
       runciePost $ diagW.comp, width = 16, height = 16)

runciePost $ quantile <-
  foreach(i = 1:4) %dopar%
DiagQuantilePop(138, runciePost $ ext [[i]], OneDef, tree = Tree [[5]])
runciePost $ quantile.comp <-
  do.call(arrangeGrob, c(runciePost $ quantile, nrow = 4))
ggsave(paste (folder, 'runcie.quantile.pdf', sep = ''),
       runciePost $ quantile.comp, width = 16, height = 32)

runciePost $ Drift <- DiagDrift (runciePost $ combine.ext)
runciePost $ Drift.comp <-
  do.call(arrangeGrob, c(runciePost $ Drift, ncol = 3))
ggsave(paste (folder, 'drift.pdf', sep = ''),
       runciePost $ Drift.comp, width = 36, height = 12)

### vamos olhar para os lambdas

runciePost $ LambdaW.df <- melt (llply (runciePost $ ext, function (L) L $ LambdaW))

runciePost $ LambdaW.plot <- 
  ggplot (runciePost $ LambdaW.df) +
  geom_violin(aes (x = interaction (trait, L1), y = value, color = L1),
              alpha = 0.5, scale = 'width') +
  facet_wrap(~ Var2, scales = 'free', ncol = 1) +
  scale_x_discrete (breaks = NULL) + theme_minimal()

ggsave(paste (folder, 'LambdaW.pdf', sep = ''),
       runciePost $ LambdaW.plot, width = 24, height = 48)

save(runciePost, file = paste (folder, 'post.RData', sep = ''))

### escrever DiagEvol

mean (runciePost $ combine.ext $ lp__) +
   0.5 * var (runciePost $ combine.ext $ lp__)

