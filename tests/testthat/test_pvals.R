test_that("nparFtest_realdata", {
  pars0 = 3
  pars1 = 6

  pAdj <- nparFtest(rss0=rss0_ref,
                    rss1=rss1_ref,
                    df_type = "theoretical",
                    n0=n0_ref,
                    n1=n1_ref,
                    pars0=pars0,
                    pars1=pars1)

  # Compute expected values:
  d1 = pars1 - pars0
  d2 = n1_ref - pars1

  f = ((rss0_ref-rss1_ref)/d1) / (rss1_ref/d2)

  p = 1 - pf(f, df1 = d1, df2 = d2)

  pAdj_exp = p.adjust(p, "BH")

  expect_equal(pAdj, pAdj_exp)

})

test_that("nparFtest_single", {

  rss0 <- 2
  rss1 <- 1
  rssDiff <- rss0 - rss1

  n0 <- 40
  n1 <- 30 # in 1 condition, only 1 replicate could be fitted

  pars0 = 3
  pars1 = 6

  # Compute expected values:
  d1 = pars1 - pars0
  d2 = n1 - pars1

  f = (rssDiff/d1) / (rss1/d2)
  p = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj_exp = p.adjust(p, "BH")

  pAdj <- nparFtest(rss0=rss0, rss1=rss1, df_type = "theoretical", n0=n0, n1=n1, pars0=pars0, pars1=pars1)

  expect_equal(pAdj, pAdj_exp)
})
