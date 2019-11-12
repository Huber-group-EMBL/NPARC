source(file.path(rprojroot::find_package_root_file(), "tests/init_tests.R"))

test_that("invokeNPARC_allok", {
  results <- runNPARC(x = ATP_targets_stauro$temperature,
                      y = ATP_targets_stauro$relAbundance,
                      id = ATP_targets_stauro$uniqueID,
                      groupsAlt = ATP_targets_stauro$compoundConcentration,
                      dfType = "theoretical")

  expect_equal(results$pAdj, pAdj_ref, tolerance = 1e-6)

})
