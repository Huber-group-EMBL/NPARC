#' Non-parametric analysis of response curves
#'
#' Wrapper function for melting curve fitting and hypothesis testing.
#'
#' @param x numeric vector of the independent variables (typically temperature)
#' @param y numeric vector of the dependent variables (typically relative abundance measurements)
#' @param control list of parameters used to control specific parts of the analyse
#' @param BPPARAM BiocParallel parameter object to invoke curve fitting in parallel. Default: BiocParallel::SerialParam()
#' @param seed Random seed to control resampling in case of unsuccessful model fits
#' @param maxAttemps number of resampling steps in case of unsuccessful model fits
#' @param df_type character value indicating the method for degrees of freedom computation for the F-test. Theoretical yields the text-book solution. Empirical yields estimates derived from the distribution moments of the RSS.
#'
#' @export
#'
#' @examples
#' data(stauro_TPP_data_tidy)
#' df <- head(stauro_TPP_data_tidy, 100)
#' testResults <- runNPARC(x = df$temperature,
#'                      y = df$relAbundance,
#'                      id = df$uniqueID,
#'                      groupsAlt = df$compoundConcentration,
#'                      df_type = "empirical")
runNPARC <- function(x, y, id,
                  groupsNull = NULL,
                  groupsAlt,
                  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                  df_type = c("theoretical", "empirical"),
                  control = getParams()){

  fits <- nparcFit(x = x,
                   y = y,
                   id = id,
                   groupsNull = groupsNull,
                   groupsAlt = groupsAlt,
                   BPPARAM = BPPARAM,
                   return_models = FALSE,
                   control = control)

    modelMetrics <- fits$metrics

    testRes <-  nparcFtest(modelMetrics, df_type = df_type)

  return(testRes)
}
