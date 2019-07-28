source("tests/init_tests.R")

test_that("invokeNPARC_allok", {
  results <- nparc(x = ATP_targets_stauro$temperature,
                   y = ATP_targets_stauro$relAbundance,
                   id = ATP_targets_stauro$uniqueID,
                   groupsNull = NULL,
                   groupsAlt = ATP_targets_stauro$compoundConcentration,
                   BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                   seed = 123,
                   maxAttempts = 100,
                   alwaysPermute = FALSE,
                   return_models = TRUE,
                   df_type = "theoretical")

  expect_equal(results$pAdj, pAdj_ref, tolerance = 1e-6)

})
