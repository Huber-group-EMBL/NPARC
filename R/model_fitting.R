performParallelFits <- function(x,
                                y,
                                iter,
                                BPPARAM,
                                seed,
                                maxAttempts,
                                alwaysPermute,
                                verbose){
  # rename to fitModels

  # ---- Fit sigmoid models ----
  models <- BiocParallel::bplapply(unique(iter),
                                   function(i) fitToSubset(subset = i,
                                                           x = x,
                                                           y = y,
                                                           iter = iter,
                                                           seed = seed,
                                                           alwaysPermute = alwaysPermute,
                                                           maxAttempts = maxAttempts,
                                                           verbose = verbose),
                                   BPPARAM = BPPARAM)
  gc()

  # ---- Return model fits ----
  model_names <- sapply(models, function(m) attributes(m)$iter)

  names(models) <- model_names

  #   message("\nFor ", sum(allRSS$repeats == 0) , " proteins, the models successfully converged and produced nonnegative RSS-differences in the first iteration.",
  #           "\nFor ", sum(allRSS$repeats > 0 & allRSS$repeats < 10), " proteins, nonnegative RSS-differences could be obtained after repeating the fit with different start parameters.")


  return(models)
}

fitToSubset <- function(subset, x, y, iter, seed, alwaysPermute, maxAttempts, verbose){

  if (verbose) message(subset)

  idx <- which(iter == subset)

  model <- repeatSingleFit(x = x[idx],
                           y = y[idx],
                           seed = seed,
                           alwaysPermute = alwaysPermute,
                           maxAttempts = maxAttempts)

  attr(model, "iter") <- subset

  return(model)

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

#' Fit sigmoid model
#'
#' @export
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



