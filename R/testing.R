combineRSS <- function(idNull, rssNull, nCoeffsNull, nFittedNull,
                       idAlt, rssAlt, nCoeffsAlt, nFittedAlt){

  rssAggrNull <- aggregateRSS(id = idNull, rss = rssNull, nCoeffs = nCoeffsNull, nFitted = nFittedNull)
  rssAggrAlt  <- aggregateRSS(id = idAlt, rss = rssAlt,  nCoeffs = nCoeffsAlt,  nFitted = nFittedAlt)

  colnames(rssAggrNull)[-1] %<>% paste0(., "Null")
  colnames(rssAggrAlt)[-1] %<>% paste0(., "Alt")

  rss <- full_join(rssAggrNull, rssAggrAlt, by = "id") %>%
    mutate(rssDiff = rssNull - rssAlt)

  return(rss)
}

aggregateRSS <- function(id, rss, nCoeffs, nFitted){

  out <- tibble(id, rss, nCoeffs, nFitted) %>%
    group_by(id) %>%
    summarise(rss = sum(rss, na.rm = TRUE),
              nCoeffs = sum(nCoeffs, na.rm = TRUE),
              nFitted = sum(nFitted, na.rm = TRUE)) %>%
    ungroup

  return(out)
}


nparcFtest <- function(id, rss0, rss1,
                      df_type = c("estimate", "theoretical"),
                      n0 = NULL, n1 = NULL, pars0 = NULL, pars1 = NULL){

  if (df_type == "theoretical"){
    d1 = pars1 - pars0
    d2 = n1 - pars1
  }

  rssDiff <- rss0 - rss1

  f = ((rssDiff)/d1) / (rss1/d2)
  pVal = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj = p.adjust(pVal, "BH")

  out <- tibble(id, rssDiff, F = f, pVal, pAdj)

  return(out)
}
