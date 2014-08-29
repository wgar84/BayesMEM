All $ post.oneSigma <- extract(All $ fit.oneSigma)

names (All $ post.oneSigma)

attach ('../../covTensor/Work/etd.all.RData')

source ('../../covTensor/Func/MeanMatrix.R')
source ('../../covTensor/Func/Pairwise.R')
source ('../../covTensor/Func/Frobenius.R')
source ('../../covTensor/Func/MatrixDist.R')
source ('../../covTensor/Func/BuildSigma.R')
source ('../../covTensor/Func/EigenTensorDecomposition.R')
source ('../../covTensor/Func/RebuildETD.R')

dim (All $ post.oneSigma $ Sigma)

All $ projectedSigma <- aaply (All $ post.oneSigma $ Sigma, 1, ProjectMatrix,
                               etd = etd.ss $ mean.etd,
                               mean.matrix = post.vcv $ ss.grand.mean, .parallel = TRUE)

color2D.matplot (All $ post.oneSigma $ Sigma [1, , ])

plot (etd.ss $ mean.etd $ projection [, 1:2], col = hsv ((1:109)/109, 1, 1, 0.75),
      pch = 20, cex = 2, xlim = c (-10, 12), ylim = c(-2.5, 2), xlab = 'PM1', ylab = 'PM2')
for (i in 1:109)
  {
    hull <- chull (t (etd.ss $ post.vcv.mean.proj [1:2, i, ]))
    polygon (t (etd.ss $ post.vcv.mean.proj [1:2, i, ]) [hull , ],
             col = hsv (i/109, 1, 1, 0.3), border = hsv (i/109, 1, 1, 0.5), lwd = 2)
  }
legend ('topright', col = hsv ((1:109)/109, 1, 1, 0.5), pch = '', lty = 1, lwd = 3,
        legend = Tree [[1]] $ tip.label, bty = 'n', cex = 0.385)

hull <- chull (All $ projectedSigma [,1:2 ])

polygon (All $ projectedSigma [hull , ],
         col = hsv (1, 1, 0, 0.3), border = hsv (1, 1, 0, 0.5), lwd = 2)
