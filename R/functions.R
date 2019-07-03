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

  return(fitStatsNull$rss)
}

nparFtest <- function(rss0, rss1,
                      df_type = c("estimate", "theoretical"),
                      n0 = NULL, n1 = NULL, pars0 = NULL, pars1 = NULL){

  pars0 = 3
  pars1 = ifelse(n1 > 40, yes = 9, no = 6)
  rssDiff <- rss0 - rss1

  if (df_type == "theoretical"){

    d1 = pars1 - pars0
    d2 = n1 - pars1
  }

  f = (rssDiff/d1) / (rss1/d2)
  pVal = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj = p.adjust(pVal, "BH")

  return(pAdj)
}


# invokeRSSdiff <- function(x, y, group, id=NULL, BPPARAM = BiocParallel::SerialParam()){
#
#   tab <- tibble(x, y, group, id)
#   ct <- 1
#
#   if (!is.null(id)) {
#
#     allRSS <- BiocParallel::bplapply(unique(id),
#                                      function(i){
#                                        idx <- which(id == i)
#                                        message(ct, ", ", i)
#                                        ct <<- ct + 1
#                                        rssDiff <- computeRSSdiff(x = x[idx],
#                                                                  y = y[idx],
#                                                                  treatment = group[idx]) %>%
#                                          mutate(id = i)
#                                      }, BPPARAM = BPPARAM)
#     allRSS <-  bind_rows(allRSS)
#   } else {
#
#     allRSS <- computeRSSdiff(x = x, y = y, treatment = group)
#
#   }
#
#   rm(ct)
#
#   message("\nFor ", sum(allRSS$repeats == 0) , " proteins, the models successfully converged and produced nonnegative RSS-differences in the first iteration.",
#           "\nFor ", sum(allRSS$repeats > 0 & allRSS$repeats < 10), " proteins, nonnegative RSS-differences could be obtained after repeating the fit with different start parameters.")
#
#   return(allRSS)
# }





# computeRSS <- function(x, y, start = c(Pl = 0, a = 550, b = 10), seed = NULL,
#                        alwaysPermute = FALSE, maxAttempts = 100){
#
#   # Start model fitting
#   fit <- repeatFits(x = x, y = y, start = start, seed = seed,
#                     alwaysPermute = alwaysPermute, maxAttempts = maxAttempts)
#
#   if (!inherits(fit, "try-error")){
#     # If model fit converged, obtain RSS and parameters
#     resid <- residuals(fit)
#     rss <- sum(resid^2, na.rm = TRUE)
#     fittedValues <- sum(!is.na(resid))
#     params <- coefficients(fit)
#     # Compute Tm according to the Eq. given in the Methods Section of the paper.
#     # It will not be used for hypothesis testing, but for Paper Fig. 4.
#     tm <- params[["a"]]/(params[["b"]] - log((1-params[["Pl"]])/(1/2 - params[["Pl"]]) - 1))
#   } else {
#     # If model fit did not converge, return default values
#     rss <- NA
#     fittedValues <- 0
#     params <- c(Pl = NA, a = NA, b = NA)
#     tm <- NA
#   }
#
#   out <- tibble(rss = rss, fittedValues = fittedValues, tm = tm,
#                 a = params[["a"]], b = params[["b"]], Pl = params[["Pl"]])
#   return(out)
# }


# computeRSSdiff <- function(x, y, treatment, maxAttempts = 100, repeatsIfNeg = 10){
#
#   rssDiff <- -1
#   repeats <- 0
#   alwaysPermute <- FALSE
#   start0 = start1 <- c(Pl = 0, a = 550, b = 10)
#
#   while((is.na(rssDiff) | rssDiff < 0) & repeats <= repeatsIfNeg){
#
#     nullResults <- computeRSS(x = x, y = y, start = start0, seed = repeats,
#                               maxAttempts = maxAttempts,
#                               alwaysPermute = alwaysPermute)
#
#     altResults <- tibble(x, y, treatment) %>%
#       group_by(treatment) %>%
#       do({
#         fit = computeRSS(x = .$x, y = .$y, start = start1, seed = repeats,
#                          maxAttempts = maxAttempts,
#                          alwaysPermute = alwaysPermute)
#       }) %>% ungroup
#
#     rss0 <- nullResults$rss
#     rss1 <- sum(altResults$rss)
#     rssDiff <- rss0 - rss1
#
#     if (is.na(rssDiff) | rssDiff < 0){
#       repeats <- repeats + 1
#       alwaysPermute <- TRUE
#       start1 <- c(Pl = nullResults[["Pl"]], a = nullResults[["a"]], b = nullResults[["b"]])
#     }
#   }
#
#   n0 <- nullResults$fittedValues
#   n1 <- sum(altResults$fittedValues)
#
#   tm <- altResults %>%
#     mutate(key = paste0("tm_", treatment)) %>%
#     dplyr::select(key, tm) %>% spread(key, tm)
#
#   out <- tibble(rss0, rss1, rssDiff, n0, n1, repeats) %>%
#     cbind(tm)
#
#   return(out)
# }
