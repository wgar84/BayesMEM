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
require (gridExtra)
require (gridBase)
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

attach ('TestNO.RData')

no.test <- list()
no.test $ ext <- llply (Test.NoPrior, function (L) L $ ext)
#no.test $ ext <- no.test $ ext [-7]

no.test $ subtree <- extract.clade(Tree [[5]], 138)

### QUANTIS
no.test $ diag.quant <- llply (no.test $ ext, function (L)
                                DiagQuantilePop(138, L, OneDef, Tree [[5]],
                                                at = c(0, 0.25, 0.75, 1)),
                                .parallel = TRUE)

pdf (file = 'no.diag.quant.pdf', width = 24, height = 12)
for (i in 1:8)
  print(no.test $ diag.quant [[i]])
dev.off (dev.cur ())

### TRACE

no.test $ combine <-
  extract (sflist2stanfit (llply (Test.NoPrior, function (L) L $ model)))

no.test $ combine <-
  llply (no.test $ combine,
         function (L)
         {
           change.k <- which (dim (L) == 39)
           if (length (change.k) > 0)
             for (i in change.k)
               {
                 dimnames (L) [[i]] <- dimnames (no.test $ ext [[1]] [[1]]) $ trait
                 names (dimnames (L)) [i] <- 'trait'
               }
           change.m <- which (dim (L) == 22)
           if (length (change.m) > 0)
             for (i in change.m)
               {
                 dimnames (L) [[i]] <- dimnames (no.test $ ext [[1]] [[1]]) $ node
                 names (dimnames (L)) [i] <- 'node'
               }
           L
         })
           
no.test $ diag.trace <- DiagTrace(no.test $ combine, OneDef [[30]] $ mean)
           
ggsave('no.diag.trace.jpg',
       do.call (arrangeGrob, c(no.test $ diag.trace, 'ncol' = 3)),
       width = 36, height = 12)

### BETAS
no.test $ post.response <-
  llply (no.test $ ext, PostDeltaZ, tree = no.test $ subtree, beta = TRUE)

no.test $ combine.response <-
  PostDeltaZ (no.test $ combine, tree = no.test $ subtree, beta = TRUE)

DiagBeta(no.test $ combine.response $ Beta, no.test $ subtree)

pdf (file = 'no.diag.beta.pdf', width = 14, height = 7)
for (i in 1:8)
  print (DiagBeta(no.test $ post.response [[i]] $ Beta, no.test $ subtree))
dev.off (dev.cur ())

### GELMAN
no.test $ convergence <- DiagGelman(no.test $ ext)

no.test $ convergence.composite <-
  do.call (arrangeGrob, c(no.test $ convergence, 'ncol' = 3))

ggsave(filename = 'no.gelmanDiag.jpg',
       plot = no.test $ convergence.composite,
       width = 21, height = 7)

save (no.test, file = 'noDiag.RData')

load ('noDiag.RData')

print (DiagW(no.test $ ext [[1]], llply (OneDef [27:38], function (L) L $ ml.vcv),
             Aux $ sample.size [27:38], TRUE))


#### W vs P
no.test $ RS <- 
  llply (no.test $ ext, DiagW, vcv.terminal.list =
         llply (OneDef [27:38], function (L) L $ ml.vcv),
         sample.size = Aux $ sample.size [27:38], 
         parallel = FALSE,
         .parallel = TRUE) ### paralelo entre cadeias, não dentro

no.test $ RS.panel <- do.call (arrangeGrob, c(no.test $ RS, 'nrow' = 2, 'ncol' = 4))

ggsave('WvsP.panel.jpg', no.test $ RS.panel, width = 32, height = 16)

### B (under construction)

no.test $ B.ml <- var (laply (OneDef [27:38], function (L) L $ mean))

no.test $ B.term <- llply (no.test $ ext,
                           function (L) aaply (L $ Xbar [, 1:12, ], 1, var))

no.test $ RS.B.term <-
  laply (no.test $ B.term, DiagB, no.test $ B.ml, parallel = FALSE, .parallel = TRUE)

boxplot (t (no.test $ RS.B.term))

### ML.B, for real
options(contrasts = c('contr.sum', 'contr.poly'))

ml.BW.calc <- list ()
ml.BW.calc $ data <- ldply (OneDef [27:38], function (L) L $ local)

ml.BW.calc $ formula <- paste (colnames (ml.BW.calc $ data) [-1], collapse = ', ')
ml.BW.calc $ formula <- paste ('cbind(', ml.BW.calc $ formula, ') ~ .id', sep = '')
ml.BW.calc $ formula <- as.formula(ml.BW.calc $ formula)

ml.BW.calc $ formula

ml.BW.calc $ model <- lm (ml.BW.calc $ formula, data = ml.BW.calc $ data)

str (ml.BW.calc $ model)

ml.BW.calc $ sigT <- var (as.matrix (ml.BW.calc $ data [-1]))
ml.BW.calc $ sigT <- ml.BW.calc $ sigT * (nrow (ml.BW.calc $ data) - 1)

ml.BW.calc $ sigW <- CalculateMatrix(ml.BW.calc $ model)
ml.BW.calc $ sigW <- ml.BW.calc $ sigW *
  (nrow (ml.BW.calc $ data) - length (unique (ml.BW.calc $ data $ .id)) - 1)

ml.BW.calc $ sigB <- ml.BW.calc $ sigT - ml.BW.calc $ sigW
ml.BW.calc $ sigB <- ml.BW.calc $ sigB / (length (unique (ml.BW.calc $ data $ .id)) - 1)

KrzCor (ml.BW.calc $ sigB, no.test $ B.ml)

plot (eigen (ml.BW.calc $ sigB) $ values)

no.test $ Krz.B <-
  laply (no.test $ ext, DiagB, B.ml = ml.BW.calc $ sigB,
         parallel = FALSE, .parallel = TRUE, ret.dim = 12)

boxplot (t (no.test $ Krz.B [, , 1]))

boxplot (t (laply (no.test $ ext, function (L)
                   aaply(L $ Sigma, 1, RandomSkewers, cov.y = ml.BW.calc $ sigW) [, 1],
                   .parallel = TRUE)))
