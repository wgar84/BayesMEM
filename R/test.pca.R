require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 12)
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

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

failMainModel <- failwith(NULL, mainModel)

pcaTest <-
  alply(1:12, 1, function (i)
        failMainModel(138, OneDef, Tree [[1]], 'local',
                      list ('mean' = OneDef [[i+26]] $ mean,
                            'vcv' = 10*OneDef[['Callithrix_kuhlii']] $ml.vcv),
                      model = 'pcaSigma_Anc',
                      pars = c('terminal', 'root', 'ancestor',
                        'SigmaW', 'SigmaB_W', 'SigmaB'),
                      control = list ('chain_id' = i),
                      warmup = 500, iter = 1000, thin = 5),
        .parallel = TRUE)


save (pcaTest, file = 'pcaTest.RData')

