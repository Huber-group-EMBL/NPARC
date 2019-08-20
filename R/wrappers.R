#' Fit null and alternative models
#'
#' @export
nparcFit <- function(x, y, id,
                     groupsNull = NULL,
                     groupsAlt,
                     BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                     seed = 123,
                     return_models = FALSE,
                     verbose = FALSE,
                     maxAttempts = 100,
                     alwaysPermute = TRUE){

  fitResNull <- invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsNull,
                                   BPPARAM = BPPARAM,
                                   seed = seed,
                                   maxAttempts = maxAttempts,
                                   alwaysPermute = alwaysPermute,
                                   return_models = return_models,
                                   verbose = verbose)

  fitResAlt <-  invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsAlt,
                                   BPPARAM = BPPARAM,
                                   seed = seed,
                                   maxAttempts = maxAttempts,
                                   alwaysPermute = alwaysPermute,
                                   return_models = return_models,
                                   verbose = verbose)

  predictions <- bind_rows(null = fitResNull$modelPredictions,
                           alternative = fitResAlt$modelPredictions,
                           .id = "modelType")

  metrics <- bind_rows(null = fitResNull$modelMetrics,
                       alternative = fitResAlt$modelMetrics,
                       .id = "modelType")

  return(list(predictions = predictions,
              metrics = metrics))

}


invokeParallelFits <- function(x, y,
                               id,
                               groups,
                               BPPARAM,
                               seed,
                               maxAttempts,
                               alwaysPermute,
                               return_models,
                               verbose){

  if (is.null(groups)){
    groups <- as.data.frame(matrix(nrow = length(id), ncol = 0))
  } else if (is.vector(groups)){
    groups <- data.frame(group = groups, stringsAsFactors = FALSE)
  }
  group_names <- colnames(groups)

  groups <- groups %>%
    mutate(id = id) %>%
    unite("iter", !!c("id", group_names), remove = FALSE)

  # ---- Fit models ----
  t1 <- Sys.time()
  message("Starting model fitting...")

  models <- performParallelFits(x = x,
                                y = y,
                                iter = groups$iter,
                                BPPARAM = BPPARAM,
                                seed = seed,
                                maxAttempts = maxAttempts,
                                alwaysPermute = alwaysPermute,
                                verbose = verbose)


  message("... complete")
  timeDiff <- Sys.time() - t1
  message("Elapsed time: ", round(timeDiff, 2), " ", units(timeDiff), "\n")

  # ---- Flag successful model fits ----
  message("Flagging successful model fits...")

  fits <- groups %>%
    group_by_(.dots = c("iter", "id", group_names)) %>%
    distinct_(.dots = groups(.)) %>%
    do(fittedModel = models[[.$iter]]) %>%
    ungroup %>%
    ## Mark proteins where model fit was not successful
    mutate(successfulFit = (sapply(fittedModel, class) != "try-error"))

  message("... complete.\n")

  # ---- Evaluate successful model fits ----
  message("Evaluating model fits...")

  fitResults <- evalModels(fits = fits,
                           x = x,
                           groups = group_names)

  if (return_models) {
    fitResults$fittedModels <- fits
  } else {
    fitResults$fittedModels <- NULL
  }

  return(fitResults)
}




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
# # functionality to permute negative differences to be added later
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
