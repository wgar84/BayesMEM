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

load ('testRuncie.RData')

runciePost <- list()
runciePost $ ext <- llply (runcie.test, extract)

runciePost $ subtree <- extract.clade (Tree [[5]], 138)

runciePost $ ext <-
  llply (runciePost $ ext, nameStanExt,
         names.list = list ('trait' = colnames (OneDef [[1]] $ local),
           'node' = runciePost $ subtree $ tip.label,
           'node' = as.character (14:23),
           'factor' = paste ('lambda', 1:6, sep = ''),
           'iterations' = 1:100))

runciePost $ gelman.plot <- DiagGelman(runciePost $ ext)

runciePost $ gelman.plot.comp <-
  do.call (arrangeGrob, c(runciePost$gelman.plot, 'ncol' = 3))
ggsave('runcie.gelman.pdf', runciePost $ gelman.plot.comp, width = 24, height = 8)

runciePost $ combine <- sflist2stanfit(runcie.test)

runciePost $ combine.ext <- extract(runciePost $ combine)

runciePost $ combine.ext <-
  nameStanExt (runciePost $ combine.ext, 
               names.list = list ('trait' = colnames (OneDef [[1]] $ local),
                 'node' = runciePost $ subtree $ tip.label,
                 'node' = as.character (14:23),
                 'factor' = paste ('lambda', 1:6, sep = ''),
                 'iterations' = 1:400))
               

runciePost $ trace <-
  DiagTrace (runciePost $ combine.ext)
runciePost $ trace.comp <-
  do.call (arrangeGrob, c(runciePost $ trace, 'nrow' = 3))
ggsave('runcie.trace.pdf', runciePost $ trace.comp, width = 24, height = 24)

runciePost $ response <- llply (runciePost $ ext, PostDeltaZ,
                                tree = runciePost $ subtree,
                                beta = TRUE)


pdf(file = 'runcie.beta.pdf', width = 16, height = 8)
llply(runciePost $ response, DiagBeta, runciePost $ subtree, slice = 'logCS')
dev.off (dev.cur ())

runciePost $ diagW <-
  llply (runciePost $ ext, DiagW,
         vcv.terminal.list = llply (OneDef[27:38],
           function (L) L $ ml.vcv),
         sample.size = Aux $ sample.size [27:38], 
         parallel = FALSE, .parallel = TRUE)

runciePost $ diagW.comp <- do.call(arrangeGrob, c(runciePost $ diagW, ncol = 2, nrow = 2))

### vamos olhar para os lambdas

par (mfrow = c(2, 2))
for(i in 1:4)
  {
    for (j in 1:4)
      {
        boxplot (runciePost $ ext [[j]] $ LambdaW [, i, ], las = 3, cex.axis = 0.75,
                 main = i, border = hsv(h = j/4, v = 0.5), add = ifelse (j != 1, T, F))
      }
    abline (h = 0, lty = 3)
  }

par (mfrow = c(2, 2))
for(i in 1:4)
  {
    for (j in 1:4)
      {
        boxplot (runciePost $ ext [[j]] $ LambdaB [, i, ], las = 3, cex.axis = 0.75,
                 main = i, border = hsv(h = j/4, v = 0.5), add = ifelse (j != 1, T, F))
      }
    abline (h = 0, lty = 3)
  }

par (mfrow = c(2, 2))
for(i in 1:4)
  {
    boxplot (runciePost $ ext [[i]] $ PsiW, las = 3, cex.axis = 0.75,
             main = i, border = hsv(h = j/4, v = 0.7))
  }

par (mfrow = c(2, 2))
for(i in 1:4)
  {
    boxplot (runciePost $ ext [[i]] $ PsiB, las = 3, cex.axis = 0.75,
             main = i, border = hsv(h = j/4, v = 0.7))
  }
