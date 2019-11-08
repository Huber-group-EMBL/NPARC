combineRSS <- function(idNull, rssNull, nCoeffsNull, nFittedNull,
                       idAlt, rssAlt, nCoeffsAlt, nFittedAlt){

  rssAggrNull <- aggregateRSS(id = idNull, rss = rssNull, nCoeffs = nCoeffsNull, nFitted = nFittedNull)
  rssAggrAlt  <- aggregateRSS(id = idAlt, rss = rssAlt,  nCoeffs = nCoeffsAlt,  nFitted = nFittedAlt)

  colnames(rssAggrNull)[-1] %<>% paste0(.data, "Null")
  colnames(rssAggrAlt)[-1] %<>% paste0(.data, "Alt")

  rss <- full_join(rssAggrNull, rssAggrAlt, by = "id") %>%
    mutate(rssDiff = rssNull - rssAlt)

  return(rss)
}

aggregateRSS <- function(id, rss, nCoeffs, nFitted){

  out <- tibble(id, rss, nCoeffs, nFitted) %>%
    group_by(id) %>%
    summarise(rss = sum(rss, na.rm = FALSE),
              nCoeffs = sum(nCoeffs, na.rm = FALSE),
              nFitted = sum(nFitted, na.rm = FALSE)) %>%
    ungroup

  return(out)
}


#' Perform F-test
#'
#' @param modelMetrics data.frame with results of the model fit in long format.
#' @param df_type character value indicating the method for degrees of freedom computation for the F-test. Theoretical yields the text-book solution. Empirical yields estimates derived from the distribution moments of the RSS.
#' @return data frame with fitted model parameters and additional columns listing e.g. residuals sum of squares of
#'  null and alterantive model and raw and adjusted p values retrieved from testing
#' @export
#' @importFrom stats pf p.adjust
#' @examples
#' data(stauro_TPP_data_tidy)
#' df <- dplyr::filter(stauro_TPP_data_tidy, grepl("CDK|GTP|CRK", uniqueID))
#' fits <- NPARCfit(x = df$temperature, 
#'                  y = df$relAbundance, 
#'                  id = df$uniqueID, 
#'                  groupsNull = NULL, 
#'                  groupsAlt = df$compoundConcentration, 
#'                  return_models = FALSE)
#' modelMetrics <- fits$metrics
#' testRes <-  NPARCtest(modelMetrics, df_type = "theoretical")                     
NPARCtest <- function(modelMetrics, df_type = c("empirical", "theoretical")){

  metricsNull <- filter(modelMetrics, .data$modelType == "null")
  metricsAlt <- filter(modelMetrics, .data$modelType == "alternative")

  metrics <- combineRSS(idNull = metricsNull$id,
                        rssNull = metricsNull$rss,
                        nCoeffsNull = metricsNull$nCoeffs,
                        nFittedNull = metricsNull$nFitted,
                        idAlt = metricsAlt$id,
                        rssAlt = metricsAlt$rss,
                        nCoeffsAlt = metricsAlt$nCoeffs,
                        nFittedAlt = metricsAlt$nFitted)

  id <- metrics$id
  rss0 <- metrics$rssNull
  rss1 <- metrics$rssAlt
  pars0 <- metrics$nCoeffsNull
  pars1 <- metrics$nCoeffsAlt
  n0 <- metrics$nFittedNull
  n1 <- metrics$nFittedAlt
  rssDiff <- metrics$rssDiff

  if (df_type == "theoretical"){

    d1 = pars1 - pars0
    d2 = n1 - pars1

  } else if (df_type == "empirical"){

    distr_pars <- estimate_df(rss1 = rss1, rssDiff = rssDiff)
    d1 <- distr_pars$d1
    d2 <- distr_pars$d2
    s0_sq <- distr_pars$s0_sq
    rssDiff = rssDiff/s0_sq
    rss1 = rss1/s0_sq

  }

  f = (rssDiff/d1) / (rss1/d2)
  pVal = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj = p.adjust(pVal, "BH")

  out <- tibble(id, rssDiff, fStat = f, pVal, pAdj, df1 = d1, df2 = d2, rssNull = rss0, rssAlt = rss1,
                nFittedNull = n0, nFittedAlt = n1, nCoeffsNull = pars0, nCoeffsAlt = pars1)

  return(out)
}

estimate_df <- function(rss1, rssDiff){!is.na(rssDiff)
  # Estimate degrees of freedom from distribution moments of the RSS

  rm_idx <- is.na(rssDiff) | (rssDiff <= 0)
  rss1 <- rss1[!rm_idx]
  rssDiff <- rssDiff[!rm_idx]

  M = median(rssDiff, na.rm = TRUE)
  V = mad(rssDiff, na.rm = TRUE)^2
  s0_sq = 1/2 * V/M
  rssDiff = rssDiff/s0_sq
  rss1 = rss1/s0_sq
  d1 = MASS::fitdistr(x = rssDiff, densfun = "chi-squared", start = list(df = 1), method = "Brent", lower = 0, upper = length(rssDiff))[["estimate"]]
  d2 = MASS::fitdistr(x = rss1, densfun = "chi-squared", start = list(df = 1),  method = "Brent", lower = 0, upper = length(rssDiff))[["estimate"]]

  out <- list(d1 = d1, d2 = d2, s0_sq = s0_sq)
  return(out)
}
