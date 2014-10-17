require (ape)
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
  do.call (arrangeGrob, c(runciePost $ trace, 'ncol' = 3))
ggsave('runcie.trace.pdf', runciePost $ trace.comp, width = 24, height = 8)

runciePost $ trace $ root
