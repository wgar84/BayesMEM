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

anc.post $ ext <- llply(anc.post $ ext, function (L)
                        { names (dimnames (L $ ancestor)) [2] <- 'node'
                          dimnames (L $ ancestor) [[2]] <- as.character(14:23)
                          L })

anc.post $ subtree <- extract.clade(Tree [[5]], 138)

### TRACE
anc.post $ diag.trace <- alply (1:8, 1, function (i)
                                DiagTraceAlt(anc.post $ ext [[i]],
                                             OneDef [[c(26:33) [i]]] $ mean),
                                .parallel = TRUE)

pdf (file = 'multi.diag.trace.anc.pdf', width = 24, height = 12)
for (i in 1:8)
  {
    test.plot <- do.call (arrangeGrob, c (anc.post $ diag.trace [[i]], ncol = 2))
    print(test.plot)
  }
dev.off (dev.cur ())

anc.post $ one.model <- sflist2stanfit(llply (Test.Anc, function (L) L $ model))
anc.post $ one.model.ext <- extract(anc.post $ one.model)

for (i in 1:3)
names (dimnames (anc.post $ one.model.ext [[i]])) <-
  names (dimnames (anc.post $ ext [[1]] [[i]]))

anc.post $ trace.one.model <- DiagTraceAlt(anc.post $ one.model.ext)

anc.post $ trace.one.model [[1]]
anc.post $ trace.one.model [[2]]


### GELMAN
anc.post $ convergence <- DiagGelman(anc.post $ ext)

anc.post $ convergence.composite <-
  do.call (arrangeGrob, c(anc.post $ convergence, 'ncol' = 3))

ggsave(filename = 'gelmanDiag.anc.jpg',
       plot = anc.post $ convergence.composite,
       width = 21, height = 7)

anc.post $ convergence $ lp.hist

save (anc.post, file = 'anc.post.RData')
