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

attach ('ancTest.RData')

anc.post <- list()
anc.post $ ext <- llply (Test.Anc, function (L) L $ ext)
#anc.post $ ext <- anc.post $ ext [-7]

anc.post $ subtree <- extract.clade(Tree [[5]], 138)

### QUANTIS
anc.post $ diag.quant <- llply (anc.post $ ext, function (L)
                                DiagQuantilePop(138, L, OneDef, Tree [[5]],
                                                at = c(0, 0.25, 0.75, 1)),
                                .parallel = TRUE)

pdf (file = 'multi.diag.quant.pdf', width = 24, height = 12)
for (i in 1:7)
  print(anc.post $ diag.quant [[i]])
dev.off (dev.cur ())

### TRACE
anc.post $ diag.trace <- alply (1:7, 1, function (i)
                                DiagTrace(anc.post $ ext [[i]],
                                          OneDef [[c(26:31, 33) [i]]] $ mean),
                                .parallel = TRUE)

pdf (file = 'multi.diag.trace.pdf', width = 36, height = 12)
for (i in 1:7)
  {
    test.plot <- do.call (arrangeGrob, c (anc.post $ diag.trace [[i]], ncol = 3))
    print(test.plot)
  }
dev.off (dev.cur ())

### BETAS
anc.post $ post.response <-
  llply (anc.post $ ext, PostDeltaZ, tree = anc.post $ subtree, beta = TRUE)

pdf (file = 'multi.diag.beta.pdf', width = 24, height = 12)
for (i in 1:7)
  {
    print (DiagBeta(anc.post $ post.response [[i]], anc.post $ subtree))
  }
dev.off (dev.cur ())

### GELMAN
anc.post $ convergence <- DiagGelman(anc.post $ ext)

anc.post $ convergence.composite <-
  do.call (arrangeGrob, c(anc.post $ convergence, 'ncol' = 3))

ggsave(filename = 'gelmanDiag.jpg',
       plot = anc.post $ convergence.composite,
       width = 21, height = 7)

anc.post $ convergence $ lp.hist

save (anc.post, file = 'anc.post.RData')
