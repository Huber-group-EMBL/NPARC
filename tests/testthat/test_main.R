source("tests/init_tests.R")

test_that("invokeNPARC_allok", {
  results <- nparc(x = ATP_targets_stauro$temperature,
                   y = ATP_targets_stauro$relAbundance,
                   id = ATP_targets_stauro$uniqueID,
                   groupsNull = NULL,
                   groupsAlt = ATP_targets_stauro$compoundConcentration,
                   BPPARAM = BPPARAM,
                   seed = 123,
                   maxAttempts = 100,
                   alwaysPermute = FALSE,
                   return_models = TRUE)

})

test_that("invokeRSSdiff_allok", {

  rssDiffs <- invokeRSSdiff(x = ATP_targets_stauro$temperature,
                            y = ATP_targets_stauro$relAbundance,
                            group = ATP_targets_stauro$compoundConcentration,
                            id = ATP_targets_stauro$uniqueID,
                            BPPARAM = BPPARAM)

  expect_equal(rssDiffs$rss0, rss0_ref)
  expect_equal(rssDiffs$rss1, rss1_ref)
  expect_equal(rssDiffs$rssDiff, rssDiff_ref)
  expect_equal(rssDiffs$n0, n0_ref)
  expect_equal(rssDiffs$n1, n1_ref)
  expect_equal(rssDiffs$repeats, repeats_ref)
  expect_equal(rssDiffs$tm_0, tm_0_ref)
  expect_equal(rssDiffs$tm_20, tm_20_ref)
})

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


# Rprof(line.profiling = TRUE)

test_that("nparc_allok_smalldata", {

system.time({
  pAdj <- nparc(x = ATP_targets_stauro$temperature,
                y = ATP_targets_stauro$relAbundance,
                group = ATP_targets_stauro$compoundConcentration,
                id = ATP_targets_stauro$uniqueID,
                BPPARAM = BPPARAM)
})
  expect_equal(pAdj, pAdj_ref, tolerance = 1e-6)
})


test_that("nparc_allok_largedata", {

  system.time({
    pAdj <- nparc(x = staurosporineTPP$temperature,
                  y = staurosporineTPP$relAbundance,
                  group = staurosporineTPP$compoundConcentration,
                  id = staurosporineTPP$uniqueID,
                  BPPARAM = BPPARAM)

  })

  # expect_equal(pAdj, pAdj_ref, tolerance = 1e-6)
})



# Rprof(NULL)

# summaryRprof(lines = "show")


# test_that("benchmark_loops", {
#
#   x <- rnorm(1e3)
#   # try ddply:
#   system.time(
#     out1 <- expand.grid(id = 1:1e4, x = x) %>% plyr::ddply("id", function(dat) t.test(dat$x)$p.val, .progress = "text")
#   )
#
#   # try adply:
#   system.time(
#     out2 <- expand.grid(id = 1:1e4, x = x) %>% as.matrix() %>% plyr::adply(., .margins = 1, function(i) t.test(x)$p.val, .progress = "text")
#   )
#
#   # try dplyr::do:
#   system.time(
#     out3 <- expand.grid(id = 1:1e4, x = x) %>% group_by(id) %>% do(data.frame(p = t.test(x)$p.val))
#   )
#
#   system.time(
#     res_ddply <- invokeRSSdiff_do(x = staurosporineTPP$temperature,
#                                   y = staurosporineTPP$relAbundance,
#                                   group = staurosporineTPP$compoundConcentration,
#                                   id = staurosporineTPP$uniqueID)
#   )
# })
