require (ape)
require (Morphometrics)
require (expm)
require (plyr)
require (plotrix)
require (doMC)
#registerDoMC (cores = 3)
registerDoMC (cores = 6)
#registerDoMC (cores = 60)
require (reshape2)
require (ggplot2)
require (rstan)
require (phytools)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

source ('../FuncR/mainModel.R')
attach ('Nodes.RData')

safe.mainModel <- failwith (NULL, mainModel, FALSE)

plot (Tree [[1]], direction = 'upwards', cex = 0.5)

### Nodes falhou nestes nós... (Sigma_bm cagada)
nodelabels(node = Nodes $ intermediate [laply (Nodes $ inter.model, is.null)],
           cex = 0.8, pch = 22, col = 'black', bg = 'blue')
nodelabels (node = Nodes $ intermediate [59],
            cex = 0.8, pch = 22, col = 'black', bg = 'red')

## Nodes $ inter.model [[59]] <-
##   safe.mainModel(Nodes $ intermediate [59],
##                  main.data = OneDef, tree = Tree [[1]], what = 'local',
##                  prior.list = list ('mean' = OneDef [['Cebus_apella']] $ mean,
##                    'vcv' = post.vcv $ ss.grand.mean),
##                  model = 'oneSigma_Anc',
##                  pars = c('Xbar', 'alpha', 'Sigma', 'Sigma_bm'),
##                  warmup = 500, iter = 1000, thin = 5)

NodesPost <- list ()

NodesPost $ inter.ok <- which (!laply (Nodes $ inter.model, is.null))
NodesPost $ inter.ok <- NodesPost $ inter.ok [-length (NodesPost $ inter.ok)]

NodesPost $ inter.ok.nodes <- Nodes $ intermediate [NodesPost $ inter.ok]
nodelabels(node = NodesPost $ inter.ok.nodes,
           cex = 0.8, pch = 22, col = 'red', bg = 'yellow')
nodelabels (node = 111)

NodesPost $ inter.ext <-
  llply (Nodes $ inter.model [NodesPost $ inter.ok],
         extract, .parallel = TRUE)



NodesPost $ inter.response <-
  foreach (i = 1:length (NodesPost $ inter.ext)) %dopar% {
    PostDeltaZ(NodesPost $ inter.ext [[i]],
               extract.clade(Tree [[1]], NodesPost $ inter.ok.nodes [i]), beta = TRUE)
  }

for (i in 1:length (NodesPost $ inter.response))
  {
    names (dimnames(NodesPost $ inter.response [[i]] $ DeltaZ)) <-
      c('iter', 'edge', 'trait')
    dimnames (NodesPost $ inter.response [[i]] $ DeltaZ) [[3]] <-
      colnames (OneDef [[1]] $ local)
    subtree <- extract.clade (Tree [[1]], NodesPost $ inter.ok.nodes [i])
    n.taxa <- length (subtree $ tip.label)
    labels <- c (subtree $ tip.label, as.character ((1:(n.taxa - 1)) + n.taxa))
    dimnames (NodesPost $ inter.response [[i]] $ DeltaZ) [[2]] <-
      paste (labels [subtree $ edge [, 1]], labels [subtree $ edge [, 2]], sep = '-')
    names (dimnames (NodesPost $ inter.response [[i]] $ Beta)) <-
      names (dimnames (NodesPost $ inter.response [[i]] $ DeltaZ))
    dimnames (NodesPost $ inter.response [[i]] $ Beta) <-
      dimnames (NodesPost $ inter.response [[i]] $ DeltaZ)
    }

boxplot (aaply (NodesPost $ inter.response [[1]] $ DeltaZ, c(1, 2), Norm),
         las = 3, cex.axis = 0.3)

pdf (file = 'betarecon.pdf', width = 12, height = 8)

for (i in 1:length (NodesPost $ inter.ok.nodes))
  {
    par (mfrow = c(1, 2))
    plot (extract.clade (Tree [[1]], NodesPost $ inter.ok.nodes [i]),
          cex = 0.5, use.edge.length = FALSE)
    edgelabels (pch = 20,
                cex = 10 *
                abs (colMeans (NodesPost $ inter.response [[i]] $ Beta [, , 'logCS'])) /
                max (colMeans (NodesPost $ inter.response [[i]] $ Beta [, , 'logCS'])),
                col =
                ifelse (colMeans (NodesPost $ inter.response [[i]] $ Beta [, , 'logCS']) > 0,
                        rgb (1, 0, 0, 0.2), rgb (0, 0, 1, 0.2)))
    edgelabels (pch =
                c('', '*') [aaply (apply (
                  NodesPost $ inter.response [[i]] $ Beta [, , 'logCS'], 2, range), 2,
                                   function (x) prod (x) > 0) + 1])
    
    boxplot (NodesPost $ inter.response [[i]] $ Beta [, , 'logCS'], las = 3, cex.axis = 0.5,
             main = 'Beta_i (logCS)', cex = 0.5)
  }

dev.off (dev.cur())

boxplot (aaply (NodesPost $ inter.ext [[1]] $ Sigma_bm, 1, function (x)  eigen (x) $ values))

dim (NodesPost $ inter.ext [[1]] $ Xbar)


boxplot (aaply (NodesPost $ inter.ext [[1]] $ Xbar [, 41:78 , ], c(2, 3), mean))

aaply (NodesPost $ inter.ext [[1]] $ Xbar [, 41:78 , ], c(2, 3), mean) [, 1] < 0
