#' Invoke model evaluation
#'
evalModels <- function(fits, groups, x){

  modelMetrics <- assessModelMetrics(fits = fits, x = x, groups = groups)
  modelPredictions <- augmentModels(fits = fits, groups = groups)

  return(list(modelMetrics = modelMetrics, modelPredictions = modelPredictions))
}

#' Augment model to obtain data frame with measurements, predictions and residuals.
#'
#'@importFrom broom augment
augmentModels <- function(fits, groups){

  # ---- Computing model predictions and residuals ----
  message("Computing model predictions and residuals ...")

  modelPredictions <- fits %>%
    filter(successfulFit) %>%
    group_by_(.dots = c("iter", "id", groups)) %>%
    do(augment(.$fittedModel[[1]])) %>%
    ungroup %>%
    right_join(fits %>% select(!!c("iter", "id", groups))) %>% # fill for unsuccessful fits
    select(-iter)

  message("... complete.\n")

  return(modelPredictions)

}

#' Invoke calculation of performance metrics per model
#'
assessModelMetrics <- function(fits, x, groups){
  message("Evaluating models ...")

  metrics <- fits %>%
    #filter(successfulFit) %>%
    group_by_(.dots = c("iter", "id", groups)) %>%
    do(assessSingleModel(nls_obj = .$fittedModel[[1]],
                         xVec = unique(x))) %>%
    ungroup %>%
    # right_join(fits %>% select(!!c("iter", "id", groups))) %>% # fill for unsuccessful fits
    select(-iter)

  message("... complete.\n")

  return(metrics)
}

#' Retrieve fitted parameters from model
#'
assessSingleModel <- function(nls_obj, xVec){

  if (any(is.na(xVec))){
    stop("Temperature vector contains missing values. Cannot integrate over these values to compute the curve area.")
  }

  tm = a = b = pl = aumc = tm_sd = rss = logL = nCoeffs = nFitted = resid_sd <- NA
  conv <- FALSE
  if (class(nls_obj) != "try-error"){
    #meltCurve <- Reduce(paste, deparse(formula(nls_obj))) %>% gsub("y ~ ", "", .)
    pars <- coefficients(nls_obj)
    nCoeffs <- length(pars)
    a <- pars[["a"]]
    b <- pars[["b"]]
    pl <- pars[["Pl"]]
    yVec <- fitted(nls_obj) + resid(nls_obj)
    tm_fct <- parse(text = "a / (b - log((1-pl)/(1/2 - pl) - 1))")
    suppressWarnings(tm <- eval(tm_fct))
    dtm_da <-D(tm_fct, "a")
    dtm_db <-D(tm_fct, "b")
    dtm_dpl <-D(tm_fct, "pl")
    suppressWarnings(grad_tm <- c(a = eval(dtm_da), b = eval(dtm_db), pl = eval(dtm_dpl))) # derivatives of 'tm' .w.r.t. 'a' and 'b'
    covmat <- try(vcov(nls_obj), silent = TRUE) # var-covar matrix of estimators 'a' and 'b'
    if (!inherits(covmat, "try-error")){
      var_tm <- grad_tm %*% covmat[c("a", "b", "Pl"), c("a", "b", "Pl")] %*% grad_tm
      tm_sd <- as.numeric(sqrt(var_tm)) # squash matrix to vector
    }
    ## Compute area under the melting curve:
    int <- try(integrate(function(x,m) predict(m, newdata = list(x=x)),
                         min(xVec), max(xVec), nls_obj), silent = TRUE)
    if (!inherits(int, "try-error")){
      aumc <- int$value
    }

    nFitted <- nobs(nls_obj)
    rss <- sum(resid(nls_obj)^2, na.rm = TRUE)
    resid_sd  <- sqrt(rss/nFitted)
    logL <- -nFitted/2 * log(2*pi*resid_sd^2) - rss/(2*resid_sd^2) #loglik <- logLik(m)
    nonNAs <- sum(!is.na(resid(nls_obj)))

    conv <- nls_obj$convInfo$isConv
  }
  return(data.frame(tm = tm, a = a,  b = b,  pl = pl, aumc = aumc,
                    resid_sd = resid_sd, rss = rss, loglik = logL,
                    tm_sd = tm_sd, nCoeffs = nCoeffs, nFitted = nFitted, conv = conv))
}


# augmentSingleModel <- function(nls_obj){
#   #' Compute predictions and residuals for a single model
#   #'
#   #' @importFrom broom augment
#
#   modelPredictions <- broom::augment(nls_obj) #%>%
#     # mutate(type = "measurement_available")
#
#   # xGrid <- seq(min(modelPredMeasurements$x), max(modelPredMeasurements$x), length.out = 100)
#   #
#   # modelPredGrid <- data.frame(x = xGrid, y = predict(nls_obj, newdata = list(x = xGrid)))  %>%
#   #   mutate(type = "measurement_not_available")
#   #
#   # out <- full_join(modelPredMeasurements, modelPredGrid, by = c("x", "y", "type")) %>%
#   #   arrange(x)
#'
#   return(modelPredictions)
# }
