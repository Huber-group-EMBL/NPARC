
invokeParallelFits <- function(x, y,
                               id,
                               groups,
                               BPPARAM,
                               seed,
                               maxAttempts,
                               alwaysPermute,
                               return_models){

  if (is.null(groups)){
    groups <- as.data.frame(matrix(nrow = length(id), ncol = 0))
  } else if (is.vector(groups)){
    groups <- data.frame(group = groups, stringsAsFactors = FALSE)
  }

  # Determing factor names for individual model fitting
  group_names <- colnames(groups)

  # Prepare long data for fitting (df and matrix format)
  dat_long <- data.frame(id, x, y, groups, stringsAsFactors = FALSE)

  # Assign unique iterator for parallelization
  dat_long <- dat_long %>% unite_("iter", c("id", group_names), remove = FALSE)

  # ---- Fit models ----
  t1 <- Sys.time()
  message("Starting model fitting...")

  models <- performParallelFits(x = dat_long$x,
                                y = dat_long$y,
                                iter = dat_long$iter,
                                BPPARAM = BPPARAM,
                                seed = seed,
                                maxAttempts = maxAttempts,
                                alwaysPermute = alwaysPermute)


  message("... complete")
  timeDiff <- Sys.time() - t1
  message("Elapsed time: ", round(timeDiff, 2), " ", units(timeDiff), "\n")

  # ---- Flag successful model fits ----
  message("Flagging successful model fits...")

  fits <- dat_long %>%
    group_by_(.dots = c("iter", "id", group_names)) %>%
    distinct_(.dots = groups(.)) %>%
    do(fitted_model = models[[.$iter]]) %>%
    ungroup %>%
    ## Mark proteins where model fit was not successful
    mutate(successfulFit = (sapply(fitted_model, class) != "try-error"))

  message("... complete.\n")

  # ---- Evaluate successful model fits ----
  message("Evaluating model fits...")

  fitResults <- evalModels(fits = fits,
                           x = x,
                           groups = group_names,
                           return_models = return_models)


  return(fitResults)
}

fitToSubset <- function(subset, x, y, iter, seed, alwaysPermute, maxAttempts){

  message(subset)

  idx <- which(iter == subset)

  model <- repeatSingleFit(x = x[idx],
                           y = y[idx],
                           seed = seed,
                           alwaysPermute = alwaysPermute,
                           maxAttempts = maxAttempts)

  attr(model, "iter") <- subset

  return(model)

}

performParallelFits <- function(x,
                                y,
                                iter,
                                BPPARAM,
                                seed,
                                maxAttempts,
                                alwaysPermute){

  # ---- Fit sigmoid models ----
  models <- BiocParallel::bplapply(unique(iter),
                                   function(i) fitToSubset(subset = i,
                                                           x = x,
                                                           y = y,
                                                           iter = iter,
                                                           seed = seed,
                                                           alwaysPermute = alwaysPermute,
                                                           maxAttempts = maxAttempts),
                                   BPPARAM = BPPARAM)
  gc()

  # ---- Return model fits ----
  model_names <- sapply(models, function(m) attributes(m)$iter)

  names(models) <- model_names

  return(models)
}

repeatSingleFit <- function(x, y,
                            start,
                            seed = NULL,
                            alwaysPermute = FALSE,
                            maxAttempts = 100){

  i <- 0
  doFit <- TRUE
  doVaryPars <- alwaysPermute
  start <- c(Pl = 0, a = 550, b = 10)

  if (!is.null(seed)){
    set.seed(seed)
  }

  while (doFit){
    startTmp <- start * (1 + doVaryPars*runif(1, -0.5, 0.5))
    m <- fitSingleSigmoid(x = x, y = y, start = startTmp)
    i <- i + 1
    doFit <- inherits(m, "try-error") & i < maxAttempts
    doVaryPars <- TRUE
  }

  return(m)
}

fitSingleSigmoid <- function(x, y, start){
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

evalModels <- function(fits, groups, x, return_models){

  modelStats       <- assessModels(fits = fits, x = x, groups = groups, return_models = return_models)
  modelPredictions <- augmentModels(fits = fits, groups = groups)

  return(list(stats = modelStats, predictions = modelPredictions))
}

augmentModels <- function(fits, groups){
  # ---- Computing model predictions and residuals ----
  message("Computing model predictions and residuals ...")

  model_table <- fits %>%
    filter(successfulFit) %>%
    group_by_(.dots = c("iter", "id", groups)) %>%
    do(augmentSingleModel(nls_obj = .$fitted_model[[1]])) %>%
    ungroup %>%
    arrange(id)
  message("... complete.\n")

  return(model_table)

}

assessModels <- function(fits, x, groups, return_models){
  # ---- Evaluate models ----
  message("Evaluating models ...")

  fit_stats <- fits %>%
    filter(successfulFit) %>%
    group_by_(.dots = c("iter", "id", groups)) %>%
    do(assessSingleModel(nls_obj = .$fitted_model[[1]],
                         xVec = unique(x))) %>%
    ungroup %>%
    arrange(id)
  message("... complete.\n")

  # ---- Create return value ----
  # Return model summaries and (optional) model objects:
  out <- fits %>%
    left_join(fit_stats) %>%
    arrange(id) %>%
    ungroup %>%
    dplyr::select(-iter)

  if (!return_models) {
    out [["fitted_model"]] <- NULL
  }

  return(out)
}

augmentSingleModel <- function(nls_obj){
  #' Compute predictions and residuals for a single model
  #'
  #' @importFrom broom augment

  modelPredMeasurements <- broom::augment(nls_obj) %>%
    mutate(type = "measurement_available")

  xGrid <- seq(min(modelPredMeasurements$x), max(modelPredMeasurements$x), length.out = 100)

  modelPredGrid <- data.frame(x = xGrid, y = predict(nls_obj, newdata = list(x = xGrid)))  %>%
    mutate(type = "measurement_not_available")

  out <- full_join(modelPredMeasurements, modelPredGrid, by = c("x", "y", "type")) %>%
    arrange(x)

  return(out)
}

assessSingleModel <- function(nls_obj, xVec){
  #' Retrieve fitted parameters from model
  #'
  #' @importFrom sme AICc

  if (any(is.na(xVec))){
    stop("Temperature vector contains missing values. Cannot integrate over these values to compute the curve area.")
  }

  meltCurve <- Reduce(paste, deparse(formula(nls_obj))) %>% gsub("y ~ ", "", .)
  deriv1 <- D(parse(text = meltCurve), "x")
  deriv2 <- D(deriv1, "x")

  tm = a = b = pl = aumc = tm_sd = rss = logL = nCoeffs = nObs = resid_sd = aicc <- NA
  conv <- FALSE
  if (class(nls_obj) != "try-error"){
    pars <- coefficients(nls_obj)
    nCoeffs <- length(pars)
    a <- pars[["a"]]
    b <- pars[["b"]]
    pl <- pars[["Pl"]]
    yVec <- fitted(nls_obj) + resid(nls_obj)
    tm_fct <- parse(text = "a / (b - log((1-pl)/(1/2 - pl) - 1))")
    suppressWarnings(tm <- eval(tm_fct))
    # suppressWarnings(slope <- TPP:::meltingCurveSlope(model = nls_obj, xInfl = tm))
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

    nObs <- nobs(nls_obj)
    rss <- sum(resid(nls_obj)^2, na.rm = TRUE)
    resid_sd  <- sqrt(rss/nObs)
    logL <- -nObs/2 * log(2*pi*resid_sd^2) - rss/(2*resid_sd^2) #loglik <- logLik(m)
    aicc <- AICc(nls_obj)
    nonNAs <- sum(!is.na(resid(nls_obj)))

    conv <- nls_obj$convInfo$isConv
  }
  return(data.frame(tm = tm, a = a,  b = b,  pl = pl, aumc = aumc,
                    resid_sd = resid_sd, rss = rss, loglik = logL,  AICc = aicc,
                    tm_sd = tm_sd, nCoeffs = nCoeffs, nObs = nObs, conv = conv))
}


