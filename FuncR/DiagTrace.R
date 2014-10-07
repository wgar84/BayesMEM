DiagTrace <- function (extraction, reference.value = NULL)
  {
    alpha.df <- melt (extraction $ alpha)
    alpha.df [, 'iterations'] <-
      as.numeric (alpha.df [, 'iterations'])

    Plots <- list()
    Plots $ alpha <-
      ggplot (alpha.df, aes (x = iterations, y = value)) +
        geom_line() +
          facet_wrap(~ trait, scales = 'free_y')

    xbar.df <- melt (extraction $ Xbar)
    xbar.df [, 'ancestor'] <- !is.na (as.numeric (as.character (xbar.df [, 'node'])))

    Plots $ extant <-
      ggplot (subset (xbar.df, !ancestor), 
              aes (x = iterations, y = value, color = node)) +
                geom_line() +
                  facet_wrap(~ trait, scales = 'free_y')

    Plots $ ancestor <-
      ggplot (subset (xbar.df, ancestor), 
              aes (x = iterations, y = value, color = node)) +
                geom_line() +
                  facet_wrap(~ trait, scales = 'free_y')

    if (!is.null(reference.value))
      {
        ref.df <- melt (reference.value)
        ref.df $ trait <- names (reference.value)

        Plots <-
          llply (Plots, function (Gr)
                 Gr + geom_hline (aes (yintercept = value), ref.df))
      }                             
    Plots
  }
