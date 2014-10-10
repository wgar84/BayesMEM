DiagGelman <- function (extract.list, parallel = TRUE)
  {
    nchain <- length (extract.list)
    if(nchain < 2)
      stop('supply at least two chains')
    theta.names <- names (extract.list [[1]])
    
    theta.loa <- alply (theta.names, 1, function (theta)
                        laply (extract.list, function (L) L [[theta]]),
                        .parallel = parallel)
    names (theta.loa) <- theta.names
    theta.loa <- llply (theta.loa,
                        function (L) {names (dimnames (L)) [1] <- 'chain'; L},
                        .parallel = parallel)
    theta.lod <- llply (theta.loa, function (L) melt (L), .parallel = parallel)
    theta.lod <-
      llply (theta.lod, function (L)
             {
               to.replace <- which (colnames (L) == 'trait')
               if (length (to.replace) > 1)
                 colnames (L) [to.replace [2]] <- 'trait2'
               with.node <- which (colnames (L) == 'node')
               if (length (with.node) > 0)
                 L $ anc <- 
                   c('ancestor', 'terminal') [is.na (
                     as.numeric (as.character(L $ node))) + 1]
               L
             }, .parallel = parallel)
    LogP <- theta.lod $ lp__
    theta.lod <- theta.lod [-length (theta.lod)]
    theta.varW.form <-
      llply (theta.lod, function (L)
             {
               colL <- colnames (L)
               formL <- colL [-c(which (colL == 'iterations'),
                                 which (colL == 'value'))]
                       formL
             }, .parallel = parallel)
    theta.varT.form <-
      llply (theta.lod, function (L)
             {
               colL <- colnames (L)
               formL <- colL [-c(which (colL == 'chain'),
                                 which (colL == 'iterations'),
                                 which (colL == 'value'))]
               formL
             }, .parallel = parallel)
    theta.varW <-
      alply (1:length (theta.lod), 1, function (i)
             ddply (theta.lod [[i]], theta.varW.form [[i]],
                    summarise,
                    var.w = var (value)), .parallel = parallel)
    theta.varW <-
      alply (1:length (theta.varW), 1, function (i)
             ddply (theta.varW [[i]], theta.varT.form [[i]],
                    summarise,
                    mean.var.w = mean (var.w)), .parallel = parallel)
    theta.varT <-
      alply (1:length (theta.lod), 1, function (i)
             ddply (theta.lod [[i]], theta.varT.form [[i]], summarise,
                    var.t = var (value)), .parallel = parallel)
    names (theta.varT) <- names (theta.varW) <- names (theta.lod)
    theta.varT <- melt (theta.varT)
    theta.varW <- melt (theta.varW)
    theta.varW $ var.t <- theta.varT $ value
    theta.var <- theta.varW
    theta.var <- theta.var [-which (colnames (theta.var) == 'variable')]
    colnames (theta.var) [which (colnames (theta.var) == 'value')] <- 'var.w'
    
    Plots <- list ()
    Plots <-
      within (Plots,
              {
                ### scatter
                wt.scatter <-
                  ggplot (theta.var) +
                    geom_point (aes (x = var.w, y = var.t, colour = L1),
                                alpha = 0.4) +
                                  scale_x_log10() +
                                    scale_y_log10() +
                                      geom_abline(intercept = 0, slope = 1)
                wt.scatter <-
                  wt.scatter +
                    theme_minimal() +
                      xlab (expression (plain (log) ~~ sigma[w])) +
                        ylab (expression (plain (log) ~~ sigma[t])) +
                          scale_colour_discrete(name = 'Parameter')

                ### var.w/var.t
                ratio.jitter <-
                  ggplot (theta.var) +
                    geom_point(aes (x = L1, y = var.w/var.t, colour = anc),
                               alpha = 0.4,
                               position = position_jitter (width = 0.25)) 
                ratio.jitter <-
                  ratio.jitter +
                    theme_minimal() +
                      xlab('Parameter') +
                        ylab (expression (frac (sigma[w], sigma[t]))) +
                          scale_colour_discrete(name = 'Parameter')


                ### LP
                lp.hist <-
                  ggplot(LogP) +
                    geom_histogram(aes(x = value, fill = factor (chain)),
                                   alpha = 0.7,
                                   position = 'identity')
                lp.hist <- lp.hist +
                  theme_minimal() +
                  xlab ('Log Probability') +
                    ylab('Count') +
                      scale_fill_discrete(name = 'Chain ID')
              })
    
    
    Plots
  }
