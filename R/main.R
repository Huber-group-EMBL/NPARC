#' Main function
#'
#' @import dplyr tidyr
#' @export
#' @examples
#' data(staurosporineTPP)
#' nparc(x = staurosporineTPP$temperature,
#'       y = staurosporineTPP$relAbundance,
#'       id = staurosporineTPP$uniqueID,
#'       groupsAlt = staurosporineTPP$compoundConcentration)
#'
nparc <- function(x, y, id,
                  groupsNull = NULL,
                  groupsAlt,
                  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                  seed = 123,
                  return_models = FALSE,
                  maxAttempts = 100,
                  alwaysPermute = TRUE,
                  df_type = c("theoretical", "estimate")){

  fitResNull <- invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsNull,
                                   BPPARAM = BPPARAM,
                                   seed = seed,
                                   maxAttempts = maxAttempts,
                                   alwaysPermute = alwaysPermute,
                                   return_models = return_models)

  fitResAlt <-  invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsAlt,
                                   BPPARAM = BPPARAM,
                                   seed = seed,
                                   maxAttempts = maxAttempts,
                                   alwaysPermute = alwaysPermute,
                                   return_models = return_models)

  fitStatsNull <- fitResNull$modelMetrics
  fitStatsAlt <- fitResAlt$modelMetrics

  rss <- combineRSS(idNull = fitStatsNull$id,
                    rssNull = fitStatsNull$rss,
                    nCoeffsNull = fitStatsNull$nCoeffs,
                    nFittedNull = fitStatsNull$nFitted,
                    idAlt = fitStatsAlt$id,
                    rssAlt = fitStatsAlt$rss,
                    nCoeffsAlt = fitStatsAlt$nCoeffs,
                    nFittedAlt = fitStatsAlt$nFitted)

  testRes <-  nparcFtest(id = rss$id,
                      rss0 = rss$rssNull,
                      rss1 = rss$rssAlt,
                      df_type = df_type,
                      n0 = rss$nFittedNull,
                      n1 = rss$nFittedAlt,
                      pars0 = rss$nCoeffsNull,
                      pars1 = rss$nCoeffsAlt)

  # to do:
  # plt <-


  return(testRes)
}
