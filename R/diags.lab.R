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
require (grid)
require (gridBase)
require (gridExtra)
require (rstan)
require (phytools)
require (geiger)
require (mvtnorm)
require (boa)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')

.source.files <- dir('../FuncR/', pattern = '.R', full.names = TRUE)
.source.files <- .source.files [!grepl ('~', .source.files)]
for (i in 1:length (.source.files))
  source (.source.files [i])

### oito (ou sete posterioris)

attach ('sevenTest.RData')

new.test <- list()
new.test $ ext <- llply (Test.new, function (L) L $ ext)
new.test $ ext <- new.test $ ext [-7]

new.test $ subtree <- extract.clade(Tree [[5]], 138)

### QUANTIS
new.test $ diag.quant <- llply (new.test $ ext, function (L)
                                DiagQuantilePop(138, L, OneDef, Tree [[5]],
                                                at = c(0, 0.25, 0.75, 1)),
                                .parallel = TRUE)

pdf (file = 'multi.diag.quant.pdf', width = 24, height = 12)
for (i in 1:7)
  print(new.test $ diag.quant [[i]])
dev.off (dev.cur ())

### TRACE
new.test $ diag.trace <- alply (1:7, 1, function (i)
                                DiagTrace(new.test $ ext [[i]],
                                          OneDef [[c(26:31, 33) [i]]] $ mean),
                                .parallel = TRUE)

pdf (file = 'multi.diag.trace.pdf', width = 36, height = 12)
for (i in 1:7)
  {
    test.plot <- do.call (arrangeGrob, c (new.test $ diag.trace [[i]], ncol = 3))
    print(test.plot)
  }
dev.off (dev.cur ())

### BETAS
new.test $ post.response <-
  llply (new.test $ ext, PostDeltaZ, tree = new.test $ subtree, beta = TRUE)

pdf (file = 'multi.diag.beta.pdf', width = 24, height = 12)
for (i in 1:7)
  {
    print (DiagBeta(new.test $ post.response [[i]], new.test $ subtree))
  }
dev.off (dev.cur ())

### GELMAN
new.test $ convergence <- DiagGelman(new.test $ ext)

new.test $ convergence.composite <-
  do.call (arrangeGrob, c(new.test $ convergence, 'ncol' = 3))

ggsave(filename = 'gelmanDiag.jpg',
       plot = new.test $ convergence.composite,
       width = 21, height = 7)

new.test $ convergence $ lp.hist

save (new.test, file = 'new.test.RData')
