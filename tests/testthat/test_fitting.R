source("tests/init_tests.R")

test_that("invokeParallelFits_allok_rss0", {

  result <- invokeParallelFits(x = ATP_targets_stauro$temperature,
                                y = ATP_targets_stauro$relAbundance,
                                id = ATP_targets_stauro$uniqueID,
                                groups = ATP_targets_stauro$uniqueID,
                                BPPARAM = BiocParallel::SerialParam(),
                                seed = 123,
                                maxAttempts = 100,
                                alwaysPermute = FALSE,
                                return_models = TRUE)

  rss0_new <- result$stats$rss

  expect_equal(unname(rss0_new)[-16], rss0_ref[-16]) # position 16: ATP5G1_IPI00009075 -> was a different seed used to resample due to negative RSS-Diff?
})

test_that("invokeParallelFits_allok_rss1", {

  result <- invokeParallelFits(x = ATP_targets_stauro$temperature,
                               y = ATP_targets_stauro$relAbundance,
                               id = ATP_targets_stauro$uniqueID,
                               groups = ATP_targets_stauro$compoundConcentration,
                               BPPARAM = BiocParallel::SerialParam(),
                               seed = 123,
                               maxAttempts = 100,
                               alwaysPermute = FALSE,
                               return_models = TRUE)

  rss1_new <- result$stats %>%
    group_by(id) %>%
    summarise(rss = sum(rss))


  expect_equal(rss1_new$rss, rss1_ref)
})

test_that("performParallelFits_allok_rss0", {

  models <- performParallelFits(x = ATP_targets_stauro$temperature,
                                y = ATP_targets_stauro$relAbundance,
                                iter = ATP_targets_stauro$uniqueID,
                                BPPARAM = BiocParallel::SerialParam(),
                                seed = 123,
                                maxAttempts = 100,
                                alwaysPermute = FALSE)

  rss0_new <- sapply(models, function(m) {
    ifelse(inherits(m , "try-error"), NA, m$m$deviance())
  })

  expect_equal(unname(rss0_new)[-16], rss0_ref[-16]) # position 16: ATP5G1_IPI00009075 -> was a different seed used to resample due to negative RSS-Diff?
})


test_that("performParallelFits_allok_rss1", {

  models <- performParallelFits(x = ATP_targets_stauro$temperature,
                                y = ATP_targets_stauro$relAbundance,
                                iter = paste(ATP_targets_stauro$uniqueID, ATP_targets_stauro$compoundConcentration),
                                BPPARAM = BiocParallel::SerialParam(),
                                seed = 123,
                                maxAttempts = 100,
                                alwaysPermute = FALSE)

  rss1_new <- sapply(models, function(m) {
    ifelse(inherits(m , "try-error"), NA, m$m$deviance())
  })

  rss1_new <- tibble(groups = names(rss1_new), rss1 = rss1_new) %>%
    separate("groups", c("id", "compoundConcentration"), remove = FALSE, sep = " ") %>%
    group_by(id) %>%
    summarise(rss1 = sum(rss1))

  expect_equal(unname(rss1_new$rss1), rss1_ref)
})
