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
