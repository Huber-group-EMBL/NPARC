library(testthat)
library(NPARC)

test_check("NPARC")

test_that("invoke_fits", {

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
