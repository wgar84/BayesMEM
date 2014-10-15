DiagTrace <- function (extraction, reference.value = NULL)
  {
    root.df <- melt (extraction $ root)
    root.df [, 'iterations'] <-
      as.numeric (root.df [, 'iterations'])

    Plots <- list()
    Plots $ root <-
      ggplot (root.df, aes (x = iterations, y = value)) +
        geom_line() +
          theme_minimal() +
            facet_wrap(~ trait, scales = 'free_y')

    terminal.df <- melt (extraction $ terminal)
    
    Plots $ terminal <-
      ggplot (terminal.df, 
              aes (x = iterations, y = value, color = node)) +
                geom_line() +
                  theme_minimal() +
                    facet_wrap(~ trait, scales = 'free_y')

    ancestor.df <- melt (extraction $ ancestor)
    
    Plots $ ancestor <-
      ggplot (ancestor.df, 
              aes (x = iterations, y = value, color = node)) +
                geom_line() +
                  theme_minimal() +
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

