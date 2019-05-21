#' Main function
#' @import dplyr tidyr
#' @export
#' @examples
#' data(staurosporineTPP)
#' nparc(x = staurosporineTPP$temperature,
#'       y = staurosporineTPP$relAbundance,
#'       group = staurosporineTPP$compoundConcentration)
#'
nparc <- function(x, y, group, id=NULL, params){

  rssDiff <- invokeRSSdiff(x = x, y = y, group = group, id = id)

  return(rssDiff)
}

invokeRSSdiff <- function(x, y, group, id=NULL){

  if (!is.null(id)) {

    tibble(x, y, group,id) %>%
      group_by(id) %>%
      do(computeRSSdiff(x = .$x, y = .$y, treatment = .$group)) %>%
      ungroup

  } else {

    rssDiff <- computeRSSdiff(x = x, y = y, treatment = group)

  }

  return(rssDiff)
}

fitSingleSigmoid <- function(x, y, start=c(Pl = 0, a = 550, b = 10)){
  try(nls(formula = y ~ (1 - Pl)  / (1+exp((b - a/x))) + Pl,
          start = start,
          data = list(x=x, y=y),
          na.action = na.exclude,
          algorithm = "port",
          lower = c(0.0, 1e-5, 1e-5),
          upper = c(1.5, 15000, 250),
          control = nls.control(maxiter=50)),
      silent = TRUE)
  }

repeatFits <- function(x, y, start = c(Pl = 0, a = 550, b = 10),
                       seed = NULL, alwaysPermute = FALSE, maxAttempts = 100){

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


computeRSS <- function(x, y, start = c(Pl = 0, a = 550, b = 10), seed = NULL,
                       alwaysPermute = FALSE, maxAttempts = 100){

  # Start model fitting
  fit <- repeatFits(x = x, y = y, start = start, seed = seed,
                    alwaysPermute = alwaysPermute, maxAttempts = maxAttempts)

  if (!inherits(fit, "try-error")){
    # If model fit converged, obtain RSS and parameters
    resid <- residuals(fit)
    rss <- sum(resid^2, na.rm = TRUE)
    fittedValues <- sum(!is.na(resid))
    params <- coefficients(fit)
    # Compute Tm according to the Eq. given in the Methods Section of the paper.
    # It will not be used for hypothesis testing, but for Paper Fig. 4.
    tm <- params[["a"]]/(params[["b"]] - log((1-params[["Pl"]])/(1/2 - params[["Pl"]]) - 1))
  } else {
    # If model fit did not converge, return default values
    rss <- NA
    fittedValues <- 0
    params <- c(Pl = NA, a = NA, b = NA)
    tm <- NA
  }

  out <- tibble(rss = rss, fittedValues = fittedValues, tm = tm,
                a = params[["a"]], b = params[["b"]], Pl = params[["Pl"]])
  return(out)
}


computeRSSdiff <- function(x, y, treatment, maxAttempts = 100, repeatsIfNeg = 10){

  rssDiff <- -1
  repeats <- 0
  alwaysPermute <- FALSE
  start0 = start1 <- c(Pl = 0, a = 550, b = 10)

  while((is.na(rssDiff) | rssDiff < 0) & repeats <= repeatsIfNeg){

    nullResults <- computeRSS(x = x, y = y, start = start0, seed = repeats,
                              maxAttempts = maxAttempts,
                              alwaysPermute = alwaysPermute)

    altResults <- tibble(x, y, treatment) %>%
      group_by(treatment) %>%
      do({
        fit = computeRSS(x = .$x, y = .$y, start = start1, seed = repeats,
                         maxAttempts = maxAttempts,
                         alwaysPermute = alwaysPermute)
      }) %>% ungroup

    rss0 <- nullResults$rss
    rss1 <- sum(altResults$rss)
    rssDiff <- rss0 - rss1

    if (is.na(rssDiff) | rssDiff < 0){
      repeats <- repeats + 1
      alwaysPermute <- TRUE
      start1 <- c(Pl = nullResults[["Pl"]], a = nullResults[["a"]], b = nullResults[["b"]])
    }
  }

  n0 <- nullResults$fittedValues
  n1 <- sum(altResults$fittedValues)

  tm <- altResults %>%
    mutate(key = paste0("tm_", treatment)) %>%
    dplyr::select(key, tm) %>% spread(key, tm)

  out <- tibble(rss0, rss1, rssDiff, n0, n1, repeats) %>%
    cbind(tm)

  return(out)
}
