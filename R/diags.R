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
attach ('../../covTensor/Work/post.vcv.RData')
attach  ('../../covTensor/Work/etd.all.RData')

source ('../FuncR/DiagQuantile.R')

Diagnostics <-
  foreach (i = 1:7) %dopar% DiagQuantile (Callithrix $ data [[i]],
             Callithrix $ post $ B [, i, ],
             Callithrix $ post [[i + 3]], at = seq(0, 1, 0.25))

pdf (file = 'diags_Callithrix.pdf', width = 15, height = 12)
for (i in 1:7)
  plot (Diagnostics [[i]] + labs (title = names (Callithrix $ data) [i]) +
        theme (axis.text = element_text (size=6)))
dev.off (dev.cur ())

Diagnostics $ B <- DiagQuantile (laply (Callithrix $ data, colMeans),
                                 Callithrix $ post $ alpha, Callithrix $ post $ V)


plot (Diagnostics $ B)
