DiagB <- function (extraction, B.ml, sample.size = NULL, parallel = TRUE, ...)
  {
    require (Morphometrics)
    B.post <- extraction $ Sigma_bm
    RS.mat <- 
      aaply(B.post, 1, function (A)
            RandomSkewers(B.ml, A, ...),
            .parallel = parallel)
    RS.mat
  }
    
  ##   names (dimnames (RS.mat)) <- c('iterations', 'terminal')
  ##   RS.df <- melt (RS.mat)
  ##   if (!is.null (sample.size))
  ##     {
  ##       ### color breaks function
  ##       colorBreaks <- function (x)
  ##         {
  ##           breaks <- c (floor (exp (quantile(log (x), probs = 0.20))),
  ##                        floor (exp (quantile(log (x), probs = 0.70))),
  ##                        floor (exp (quantile(log (x), probs = 0.95))))
  ##           names (breaks) <- NULL
  ##           breaks
  ##         }
        
  ##       ### errorbar
  ##       RS.summary <- ddply(RS.df, .(terminal), summarize,
  ##                           'mean' = mean (value), 
  ##                           'ci.m' = quantile (value, probs = c(0.025)),
  ##                           'ci.p' = quantile (value, probs = c(0.975)))
  ##       RS.summary <- cbind (RS.summary, sample.size)
  ##       RS.plot <-
  ##         ggplot (RS.summary) +
  ##           geom_point (aes(x = terminal, y = mean, colour = sample.size)) +
  ##             geom_errorbar(aes(x = terminal, ymin = ci.m, ymax = ci.p,
  ##                               colour = sample.size), width = 0.25) +
  ##                                 scale_colour_continuous(
  ##                                   low = 'red', high = 'green',
  ##                                   trans = 'log',
  ##                                   breaks = colorBreaks (sample.size),
  ##                                   name = 'Sample Size')
  ##       RS.plot <-
  ##         RS.plot +
  ##           theme_minimal() +
  ##             theme(axis.text.x = element_text(angle = 90, size = 6)) +
  ##               xlab ('OTU') + ylab ('Random Skewers OTU vs W')
  ##     }
  ##   else
  ##     {
  ##       RS.summary <- ddply(RS.df, .(terminal), summarize,
  ##                           'mean' = mean (value), 
  ##                           'ci.m' = quantile (value, probs = c(0.025)),
  ##                           'ci.p' = quantile (value, probs = c(0.975)))
  ##       RS.plot <-
  ##         ggplot (RS.summary) +
  ##           geom_point (aes(x = terminal, y = mean)) +
  ##             geom_errorbar(aes(x = terminal, ymin = ci.m, ymax = ci.p), width = 0.25)
                                
  ##       RS.plot <-
  ##         RS.plot +
  ##           theme_minimal() +
  ##             theme(axis.text.x = element_text(angle = 90, size = 6)) +
  ##               xlab ('OTU') + ylab ('Random Skewers OTU vs W')
  ##     }
  ##   RS.plot
  ## }
