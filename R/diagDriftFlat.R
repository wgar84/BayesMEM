require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
registerDoMC (cores = 4)
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

## Aux $ def.hyp <- rbind (rep (0, 8), Aux $ ed.hyp [[1]] [-20, ])
## rownames (Aux $ def.hyp) [1] <- 'logCS'
## save (Aux, file = '../../Databases/Aux.RData')
## rm (Aux)
## detach(file:../../Databases/Aux.RData)

load ('driftFlat.RData')

system('mkdir DriftFlat')
folder <- 'DriftFlat/'

post <- list()
post $ ext <- llply (fit, extract)

data $ node <- 138
data $ subtree <- extract.clade (Tree [[5]], data $ node)

post <-
  within(post,
         {
           ### ARRANGE
           ext <-
             llply (ext, nameStanExt,
                    names.list = list ('trait' = colnames (OneDef [[1]] $ local),
                      'node' = data $ subtree $ tip.label,
                      'node' = as.character ((data $ m) + 2:(data $ m - 1)),
                      'factor' = paste ('lambda', 1:data $ n_fac, sep = ''),
                      'iterations' = length (ext [[1]] $ lp__)))
           ## ext <-
           ##   llply (ext, function (L)
           ##          {
           ##            dimnames (L $ ancestor) [[2]] <-
           ##              as.character ((data $ m) + 2:(data $ m - 1))
           ##            names (dimnames (L $ ancestor)) [2] <- 'node'
           ##            L
           ##          })
                    
           combine.fit <- sflist2stanfit(fit)

           combine.ext <- extract(combine.fit)

           combine.ext <-
             nameStanExt (combine.ext, 
                          names.list = list ('trait' = colnames (OneDef [[1]] $ local),
                            'node' = data $ subtree $ tip.label,
                            'node' = as.character ((data $ m) + 2:(data $ m - 1)),
                            'factor' = paste ('lambda', 1:data $ n_fac, sep = ''),
                            'iterations' = length (post $ combine.ext $ lp__)))

           ## names (dimnames (combine.ext $ ancestor)) [2] <- 'node'
           ## dimnames (combine.ext $ ancestor) [[2]] <-
           ##   as.character ((data $ m) + 2:(data $ m - 1))
         })




### DIAGNOSTICS
post <-
  within(post,
         {
           ### GELMAN
           gelman.plot <- DiagGelman(ext)
           gelman.plot.comp <-
             do.call(arrangeGrob, c(gelman.plot, 'nrow' = 3))
           ggsave(paste (folder, 'gelman.pdf', sep = ''),
                  gelman.plot.comp, width = 8, height = 24)
              
           ### TRACE
           trace <- DiagTrace (combine.ext)
           trace.comp <- do.call (arrangeGrob, c(trace, 'nrow' = 3))
           ggsave(paste (folder, 'trace.pdf', sep = ''),
                  trace.comp, width = 24, height = 24)

           ### BETA
           response <- llply (ext, PostDeltaZ, tree = data $ subtree, beta = TRUE)
           combine.response <- PostDeltaZ(combine.ext, data $ subtree, TRUE)

           pdf(file = paste (folder, 'beta.pdf', sep = ''),
               width = 16, height = 8)
           llply(response, DiagBeta, data $ subtree, slice = 'logCS')
           dev.off (dev.cur ())

           pdf(file = paste (folder, 'combine.beta.pdf', sep = ''),
               width = 16, height = 8)
           DiagBeta(combine.response, data $ subtree, 'logCS')
           dev.off (dev.cur ())

           ### W
           diagW <- llply (ext, DiagW,
                           vcv.terminal.list = llply (OneDef[data $ subtree $ tip.label],
                             function (L) L $ ml.vcv),
                           sample.size = Aux $ sample.size [data $ subtree $ tip.label], 
                           parallel = FALSE, .parallel = TRUE)

           diagW.comp <- do.call(arrangeGrob, c(diagW, ncol = 2, nrow = 2))
           ggsave(paste (folder, 'compW.pdf', sep = ''),
                  diagW.comp, width = 16, height = 16)

           ### QUANTILE
           quantile <- foreach(i = 1:4) %dopar%
           DiagQuantilePop(data $ node, ext [[i]], OneDef, tree = Tree [[5]])
           quantile.comp <- do.call(arrangeGrob, c(quantile, nrow = 4))
           ggsave(paste (folder, 'quantile.pdf', sep = ''),
                  quantile.comp, width = 16, height = 32)

           ### DRIFT
           drift <- DiagDrift (combine.ext)
           drift.comp <-
             do.call(arrangeGrob, c(drift, ncol = 2, nrow = 2))
           ggsave(paste (folder, 'drift.pdf', sep = ''),
                  drift.comp, width = 24, height = 24)
         })

post <-
  within(post,
         {
           ### LAMBDA
           LambdaW.df <- melt (llply (ext, function (L) L $ LambdaW))

           LambdaW.plot <- 
             ggplot (LambdaW.df) +
               geom_violin(aes (y = value, x = L1, color = L1, fill = L1),
                           alpha = 0.5, scale = 'width') +
                             scale_x_discrete(breaks = NULL) +
                             facet_grid(factor ~ trait, scales = 'free') +
                               theme_bw() +
                                 theme(axis.text.x = element_text(angle = 90))

           ggsave(paste (folder, 'LambdaW.pdf', sep = ''),
                  LambdaW.plot, width = 24, height = 48)
         })

save(post, file = paste (folder, 'post.RData', sep = ''))

rm (list = ls())

### escrever DiagEvol
