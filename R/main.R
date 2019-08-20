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
                  df_type = c("theoretical", "empirical")){

  fits <- nparcFit <- function(x = x,
                               y = y,
                               id = id,
                               groupsNull = groupsNull,
                               groupsAlt = groupsAlt,
                               BPPARAM = BPPARAM,
                               seed = seed,
                               return_models = return_models,
                               verbose = verbose,
                               maxAttempts = maxAttempts,
                               alwaysPermute = alwaysPermute)

    modelMetrics <- fits$metrics

    testRes <-  nparcFtest(modelMetrics, df_type = df_type)

  # to do:
  # plt <-


  return(testRes)
}
