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
require (rstan)
require (phytools)
require (geiger)
require (mvtnorm)

##require (devtools)
##install_github('Statisticat/LaplacesDemonCpp')

require (LaplacesDemonCpp)
require (MCMCpack)

attach ('../../Databases/Reference.RData')
attach ('../../Databases/OneDef/ED.RData')
attach ('../../Databases/OneDef/OneDef.RData')
attach ('../../Databases/Tree.RData')
attach ('../../Databases/Aux.RData')
attach ('../../covTensor/Work/post.vcv.RData')

subtree <- extract.clade(Tree [[5]], 138)

Data.example <- list()
Data.example <-
  within(Data.example,
         {
           
           C <- solve(vcvPhylo(subtree))
           k <- 39
           m <- length(subtree $ tip.label)
           n.taxa <- nrow(C)
           X <- llply (OneDef [subtree $ tip.label], function (L) L $ local [, 1:k])
           N <- sum (laply (X, nrow))
           priorX <- c(4.5, rep (0, times = k-1))
           priorS <- post.vcv $ ss.grand.mean [1:k, 1:k]
           parm.names <-
             as.parm.names (
               list ('terminal' = array (0, c(m, k)),
                     'root' = rep (0, times = k),
                     'ancestor' = array (0, c(m - 2, k)),
                     'sigmaW' = rep (0, times = k),
                     'sigmaB' = rep (0, times = k),
                     'GammaW' = array (0, c(k, k)),
                     'GammaB' = array (0, c(k, k))),
               uppertri = c(0, 0, 0, 0, 0, 1, 1))
           mon.names = c ('LP', as.parm.names (
               list ('terminal' = array (0, c(m, k)),
                     'root' = rep (0, times = k),
                     'ancestor' = array (0, c(m - 2, k)),
                     'SigmaW' = array (0, c(k, k)),
                     'SigmaB' = array (0, c(k, k))),
               uppertri = c(0, 0, 0, 1, 1)))
         })

oneSigmaAnc <- function (parm, Data)
  {
    with (Data,
          {
### Parameter Reconstruction
            terminal <- array (parm [grepl ('terminal', parm.names)], c(m, k))
            root <- parm [grepl ('root', parm.names)]
            ancestor <- array (parm [grepl ('ancestor', parm.names)], c(m - 2, k))
            sigmaW <- parm [grepl ('sigmaW', parm.names)]
            sigmaB <- parm [grepl ('sigmaB', parm.names)]
            GammaW <- array (0, c(k, k))
            GammaW[upper.tri(GammaW, diag=TRUE)] <-
              parm [grepl ('GammaW', parm.names)]
            GammaB <- array (0, c(k, k))
            GammaB[upper.tri(GammaB, diag=TRUE)] <-
              parm [grepl ('GammaB', parm.names)]
            
            ### exponential (lower bounds for variances)
            sigmaB <- exp (sigmaW)
            sigmaW <- exp (sigmaB)
            diag (GammaB) <- exp (diag (GammaB))
            diag (GammaW) <- exp (diag (GammaW))
            
            ### Prior
            PriorLogDensity <-
              c (dmvnorm(terminal, priorX, priorS, log = TRUE),
                 dmvnorm(root, priorX, priorS, log = TRUE),
                 dmvnorm(ancestor, priorX, priorS, log = TRUE),
                 ###log(diwish(t (GammaB) %*% GammaB, k+1, diag(k))), 
                 ###log(diwish(t (GammaW) %*% GammaW, k+1, diag(k))),
                 dchisq(sigmaB, k+1, log = TRUE), +
                 dchisq(sigmaW, k+1, log = TRUE))

            ### VCV
            SigmaB <- diag (sigmaB) %*% t(GammaB) %*% GammaB %*% diag (sigmaB)
            SigmaW <- diag (sigmaW) %*% t(GammaW) %*% GammaW %*% diag (sigmaW)

            invPriorS <- solve(priorS)

            
            ### Likelihood
            ldC <- - sum (log (eigen (C) $ values))
            ldSigB <- sum (log (eigen (SigmaB) $ values))

            LogLik <-
              c (aaply (1:length (X), 1,
                     function (i) sum (dmvnorm (X [[i]], terminal [i, ],
                                                SigmaW, log = TRUE))))

            ### Means
            term.center <- aaply (terminal, 1, function (x) x - root)
            anc.center <- aaply (ancestor, 1, function (x) x - root)
            Xbar <- rbind (term.center, anc.center)
            LogLik <-
              c (aaply (1:length (X), 1,
                     function (i) sum (dmvnorm (X [[i]], terminal [i, ],
                                                SigmaW, log = TRUE))), 

                 - 0.5 * (sum (diag (C %*% Xbar %*% solve(SigmaB) %*% t (Xbar))) +
                          k * ldC + n.taxa * ldSigB))
            
            ### LogPosterior
            
            LP <- sum (PriorLogDensity) + sum (LogLik)

            #print (LP)
            
            mon.parm <- as.parm.names (
              list ('terminal' = array (0, c(m, k)),
                    'root' = rep (0, times = k),
                    'ancestor' = array (0, c(m - 2, k)),
                    'SigmaW' = array (0, c(k, k)),
                    'SigmaB' = array (0, c(k, k))),
              uppertri = c(0, 0, 0, 1, 1))
            
            mon.parm <- paste ('c(', paste (mon.parm, collapse = ', '), ')', sep = '')
            Monitor.out <- c (LP, eval (parse(text = mon.parm)))
            
            Modelout <- list (LP = LP, Dev = -2 * sum (LogLik),
                              Monitor = Monitor.out,
                              yhat = terminal [1, 1], parm = parm)
          })
  }

Init.example <- function ()
  {
    c(t(rmvnorm(mean = Data.example $ priorX, sigma = Data.example $ priorS,
                n = 2 * Data.example $ m - 1)),
      rchisq(2 * Data.example $k, Data.example $k + 1),
      chol(riwish(Data.example $ k + 1, Data.example $ priorS)) [upper.tri(
        diag (Data.example $ k), TRUE)],
      chol(riwish(Data.example $k + 1, Data.example $ priorS)) [upper.tri(
        diag (Data.example $ k), TRUE)])
  }



LD.test <- LaplacesDemon.hpc(oneSigmaAnc, Data.example, Init.example,
                             Iterations = 1000,
                             Thinning = 5, Status = 100,
                             Algorithm = 'HARM', CPUs = 3)

plot (LD.test)
LD.test $ Posterior1

save(Data.example, Init.example, LD.test, oneSigmaAnc, file = 'DemonTest.RData')
