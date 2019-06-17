fitSigmoidsParallel <- function(x, y, group, BPPARAM, seed, maxAttempts, alwaysPermute){

  # ---- Fit sigmoid models ----
  models <- BiocParallel::bplapply(X = unique(dat$iter),
                                   function(i){
                                     current_data <- dat[which(dat$iter == i), ]
                                     model <- repeatFits(x = current_data[["x"]],
                                                         y = current_data[["y"]],
                                                         seed = seed,
                                                         alwaysPermute = alwaysPermute,
                                                         maxAttempts = maxAttempts)
                                     attr(model, "iter") <- i
                                     return(model)
                                   }, BPPARAM = BPPARAM)
  gc()

  # ---- Return model fits ----
  model_names <- sapply(models, function(m) attributes(m)$iter)

  names(models) <- model_names

  return(models)
}

repeatFits <- function(x, y, seed = NULL, alwaysPermute = FALSE, maxAttempts = 100){

  i <- 0
  doFit <- TRUE
  doVaryPars <- alwaysPermute

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

fitSingleSigmoid <- function(x, y){
  try(nls(formula = y ~ (1 - Pl)  / (1+exp((b - a/x))) + Pl,
          start = c(Pl = 0, a = 550, b = 10),
          data = list(x=x, y=y),
          na.action = na.exclude,
          algorithm = "port",
          lower = c(0.0, 1e-5, 1e-5),
          upper = c(1.5, 15000, 250),
          control = nls.control(maxiter=50)),
      silent = TRUE)
}

