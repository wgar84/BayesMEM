nameStanExt <- function (stan.ext, names.list)
  {
    stan.ext <- llply (stan.ext,
                       function (L)
                       {
                         for (i in 1:length (names.list))
                           {
                             change <- which (dim (L) == length (names.list [[i]]))
                             if (length (change) > 0)
                               for (j in change)
                                 {
                                   dimnames (L) [[j]] <- names.list [[i]]
                                   names (dimnames (L)) [j] <- names (names.list) [i]
                                 }
                           }
                         L
                       })
    stan.ext
  }
