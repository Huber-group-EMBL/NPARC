source(file.path(rprojroot::find_package_root_file(), "tests/init_tests.R"))

test_that("nparcFTest_realdata", {

  metrics <- data.frame(modelType = c(rep("null", each = length(rss0_ref)), rep("alternative", each = length(rss1_ref))),
                        id = c(unique(atp_ids), unique(atp_ids)),
                        rss = c(rss0_ref, rss1_ref),
                        nFitted = c(n0_ref, n1_ref),
                        nCoeffs = c(rep(3, each = length(rss0_ref)), rep(6, each = length(rss1_ref))))

  testRes <- NPARCtest(metrics, dfType = "theoretical")

  # Compute expected values:
  d1 = 6 - 3
  d2 = n1_ref - 6
  f = ((rss0_ref-rss1_ref)/d1) / (rss1_ref/d2)
  p = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj_exp = p.adjust(p, "BH")

  expect_equal(testRes$pAdj, pAdj_exp)

})

test_that("nparcFTest_single", {

  metrics <- data.frame(modelType = c("null", "alternative"),
                        id = c(1,1),
                        rss = c(2, 1),
                        nFitted = c(40, 30),
                        nCoeffs = c(3,6))

  # Compute expected values:
  d1 = 6 - 3
  d2 = 30 - 6
  f = ((2-1)/d1) / (1/d2)
  p = 1 - pf(f, df1 = d1, df2 = d2)
  pAdj_exp = p.adjust(p, "BH")

  pAdj <- NPARCtest(modelMetrics = metrics, dfType = "theoretical")$pAdj

  expect_equal(pAdj, pAdj_exp)
})
