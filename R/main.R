#' Main function
#' @import dplyr tidyr
#' @export
#' @examples
#' data(staurosporineTPP)
#' nparc(x = staurosporineTPP$temperature,
#'       y = staurosporineTPP$relAbundance,
#'       group = staurosporineTPP$compoundConcentration)
#'
nparc <- function(x, y, id,
                  groupsNull,
                  groupsAlt,
                  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                  seed,
                  return_models,
                  maxAttempts,
                  alwaysPermute,
                  df_type = c("theoretical", "estimate")){

  fitStatsNull <- invokeParallelFits(x = x,
                                     y = y,
                                     id = id,
                                     groups = groupsNull,
                                     BPPARAM = BPPARAM,
                                     seed = seed,
                                     maxAttempts = maxAttempts,
                                     alwaysPermute = alwaysPermute,
                                     return_models = return_models)

  fitStatsAlt <-  invokeParallelFits(x = x,
                                     y = y,
                                     id = id,
                                     groups = groupsAlt,
                                     BPPARAM = BPPARAM,
                                     seed = seed,
                                     maxAttempts = maxAttempts,
                                     alwaysPermute = alwaysPermute,
                                     return_models = return_models)


  # to do: add next steps:
  # plt <-

  # pAdj <-  nparFtest(rss0 = rssDiff$rss0,
  #                       rss1 = rssDiff$rss1,
  #                       df_type = df_type,
  #                       n0 = rssDiff$n0,
  #                       n1 = rssDiff$n1,
  #                       pars0 = rssDiff$pars0,
  #                       pars1 = rssDiff$pars1) # tb filled
  # # return p-values, and later data.frame() with all stats

  return(fitStatsNull)
}
