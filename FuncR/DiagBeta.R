DiagBeta <- function (PostBeta, Tree, slice = 'logCS')
  {
    
    Beta.df <- melt (PostBeta)
    Beta.sub <- subset (Beta.df, trait == slice)

    Beta.summary <- ddply (Beta.sub, .(edge), summarize,
                           mean = mean (value),
                           ci.m = quantile(value, probs = 0.025),
                           ci.p = quantile(value, probs = 0.975))

    max.upper <- max (Beta.summary [2:4])
    
    Beta.summary $ zero.out <-
      Beta.summary $ ci.m * Beta.summary $ ci.p > 0

    plot.new()
    gr.lout <- grid.layout(nrow = 1, ncol = 2) 
    vp.phylo <- viewport(layout.pos.col = 1)
    vp.gg <- viewport(layout.pos.col = 2)
    pushViewport(viewport(layout = gr.lout))

    pushViewport(vp.phylo)
    par(new = TRUE, fig = gridFIG())
    plot.phylo(Tree, type = 'cladogram', label.offset = 0.1,
               direction = 'rightwards', use.edge.length = FALSE, cex = 0.75,
               edge.width = c(1, 3)[Beta.summary $ zero.out + 1])
    nodelabels(frame = 'none', adj = -1, cex = 0.75)
    popViewport()

    pushViewport(vp.gg)
    Beta.plot <-
      ggplot(Beta.summary) +
        geom_point(aes (x = edge, y = mean), size = 3) +
          geom_errorbar(aes(x = edge, ymin = ci.m, ymax = ci.p), width = 0.4) +
            theme_minimal() + labs(title = paste ('Selection for ', slice, sep = '')) +
              theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1))
    
    ## Beta.plot <- 
    ##   Beta.plot + annotate('text', x = (1:nrow (Beta.summary)) [Beta.summary $ zero.out],
    ##                        y = max.upper + (max.upper / 10),
    ##                        label = rep ('*', times = sum (Beta.summary $ zero.out)),
    ##                        size = 7)

    print(Beta.plot, newpage = FALSE)
    popViewport()
  }
