#' Fit null and alternative models for Non-parametric analysis of response curves
#'
#' Fit melting curve and return model metrics as well as predictions for the null and alternative models.
#'
#' @param x numeric vector of the independent variables (typically temperature)
#' @param y numeric vector of the dependent variables (typically relative abundance measurements)
#' @param id character vector with the protein ID to which each each data point belongs.
#' @param control list of parameters used to control specific parts of the analyse
#' @param BPPARAM BiocParallel parameter object to invoke curve fitting in parallel. Default: BiocParallel::SerialParam()
#' @param returnModels boolean value. If true, the fitted models are returned together with the test results
#' @param groupsNull one or more vectors with grouping variables for the null models. See details.
#' @param groupsAlt one or more vectors with grouping variables for the alternative models. See details.
#' @return data frame with fitted model parameters and additional columns listing e.g. residuals sum of squares 
#' @details
#' \code{groupsNull} or \code{groupsAlt} can either be a single vector each, or data.frames of the same length as \code{x} and \code{y} with one column per factor
#'
#' @export
#' @examples
#' data(stauro_TPP_data_tidy)
#' df <- dplyr::filter(stauro_TPP_data_tidy, grepl("CDK|GTP|CRK", uniqueID))
#' testResults <- NPARCfit(x = df$temperature,
#'                      y = df$relAbundance,
#'                      id = df$uniqueID,
#'                      groupsAlt = df$compoundConcentration)
#'
#'
#' @export
NPARCfit <- function(x, y,
                     id,
                     control = getParams(),
                     groupsNull = NULL,
                     groupsAlt,
                     BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                     returnModels = FALSE){

  fitResNull <- invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsNull,
                                   BPPARAM = BPPARAM,
                                   maxAttempts = control$maxAttempts,
                                   returnModels = returnModels,
                                   start = control$start)

  fitResAlt <-  invokeParallelFits(x = x,
                                   y = y,
                                   id = id,
                                   groups = groupsAlt,
                                   BPPARAM = BPPARAM,
                                   maxAttempts = control$maxAttempts,
                                   returnModels = returnModels,
                                   start = control$start)

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
                               returnModels,
                               start){

  if (is.null(groups)){
    groups <- as.data.frame(matrix(nrow = length(id), ncol = 0))
  } else if (is.vector(groups)){
    groups <- data.frame(group = groups, stringsAsFactors = FALSE)
  }
  group_names <- colnames(groups)

  groups <- groups %>%
    as_tibble() %>% 
    mutate(id = id) %>%
    unite("iter", !!c("id", group_names), remove = FALSE)

  # ---- Fit models ----
  t1 <- Sys.time()
  message("Starting model fitting...")

  models <- fitAllModels(x = x,
                         y = y,
                         iter = groups$iter,
                         BPPARAM = BPPARAM,
                         maxAttempts = maxAttempts,
                         start = start)


  message("... complete")
  timeDiff <- Sys.time() - t1
  message("Elapsed time: ", round(timeDiff, 2), " ", units(timeDiff), "\n")

  # ---- Flag successful model fits ----
  message("Flagging successful model fits...")

  fits <- groups %>%
    group_by_at(c("iter", "id", group_names)) %>%
    distinct_at(c("iter", "id", group_names)) %>%
    do(fittedModel = models[[.data$iter]]) %>%
    ungroup %>%
    ## Mark proteins where model fit was not successful
    mutate(successfulFit = 
               (vapply(.data$fittedModel, class, FUN.VALUE = "nls") != 
                    "try-error"))

  message("... complete.\n")

  # ---- Evaluate successful model fits ----
  message("Evaluating model fits...")

  fitResults <- evalModels(fits = fits,
                           x = x,
                           groups = group_names)

  if (returnModels) {
    fitResults$fittedModels <- fits
  } else {
    fitResults$fittedModels <- NULL
  }

  return(fitResults)
}

fitAllModels <- function(x,
                         y,
                         iter,
                         BPPARAM,
                         seed,
                         maxAttempts,
                         start){

  # ---- Fit sigmoid models ----
  models <- BiocParallel::bplapply(unique(iter),
                                   function(i) fitToSubset(subset = i,
                                                           x = x,
                                                           y = y,
                                                           iter = iter,
                                                           maxAttempts = maxAttempts,
                                                           start = start),
                                   BPPARAM = BPPARAM)
  gc()

  # ---- Return model fits ----
  model_names <- vapply(models, function(m) attributes(m)$iter, 
                        FUN.VALUE = "geneName_ipiId")

  names(models) <- model_names

  #   message("\nFor ", sum(allRSS$repeats == 0) , " proteins, the models successfully converged and produced nonnegative RSS-differences in the first iteration.",
  #           "\nFor ", sum(allRSS$repeats > 0 & allRSS$repeats < 10), " proteins, nonnegative RSS-differences could be obtained after repeating the fit with different start parameters.")


  return(models)
}

fitToSubset <- function(subset, x, y, iter, seed, maxAttempts, start){

  idx <- which(iter == subset)

  model <- repeatSingleFit(x = x[idx],
                           y = y[idx],
                           maxAttempts = maxAttempts,
                           start = start)

  attr(model, "iter") <- subset

  return(model)

}


repeatSingleFit <- function(x, y,
                            start,
                            maxAttempts = 100){
  # Wrapper to repeat the fit until model has converged

  i <- 0
  doFit <- TRUE
  doVaryPars <- FALSE

  while (doFit){
    startTmp <- start * (1 + doVaryPars*runif(1, -0.5, 0.5))
    m <- fitSingleSigmoid(x = x, y = y, start = startTmp)
    i <- i + 1
    doFit <- inherits(m, "try-error") & i < maxAttempts
    doVaryPars <- TRUE
  }

  return(m)
}

#' Fit sigmoid model
#'
#' @param x numeric vector of the independent variables (typically temperature)
#' @param y numeric vector of the dependent variables (typically relative abundance measurements)
#' @param start numeric vector of start parameters for the melting curve equation
#' @return model summary of type "nls"
#' @export
#' @details
#' Fits the following function to the data:
#' \eqn{y = (1 - Pl)  / (1+exp((b - a/x))) + Pl}
#' @examples
#' data(stauro_TPP_data_tidy)
#' stk4 <- dplyr::filter(stauro_TPP_data_tidy, grepl("STK4", uniqueID))
#' fitSingleSigmoid(stk4$temperature, stk4$relAbundance)
fitSingleSigmoid <- function(x, y, start=c(Pl = 0, a = 550, b = 10)){
  try(nls(formula = y ~ (1 - Pl)  / (1+exp((b - a/x))) + Pl,
          start = start,
          data = list(x=x, y=y),
          na.action = na.omit,
          algorithm = "port",
          lower = c(0.0, 1e-5, 1e-5),
          upper = c(1.5, 15000, 250),
          control = nls.control(maxiter=50)),
      silent = TRUE)
}

#' Control parameters for model fitting
#'
#' @param maxAttempts Number of resampling steps in case of unsuccessful model fits
#' @param start Numeric vector of start parameters for the melting curve equation
#' @return list of two elements: 1) "start" listing the starting parameters for
#'  melting curve fitting, 2) "maxAttempts" listing the maximal number of attempts
#'  the fit should be allowed
#' @examples
#' data(stauro_TPP_data_tidy)
#' df <- dplyr::filter(stauro_TPP_data_tidy, grepl("MAPK|ATP|CDK|GTP|CRK", uniqueID))
#' testResults <- runNPARC(x = df$temperature,
#'                      y = df$relAbundance,
#'                      id = df$uniqueID,
#'                      groupsAlt = df$compoundConcentration,
#'                      dfType = "empirical",
#'                      control = getParams(maxAttempts = 50))
#'
#' @export
#'
getParams <- function(start = c(Pl = 0, a = 550, b = 10),
                      maxAttempts = 100){

  list(start = start,
       maxAttempts = maxAttempts)
}

