source("tests/init_tests.R")

test_that("fitSigmoidsParallel_allok_rss0", {

  models <- fitWrapper(x = ATP_targets_stauro$temperature,
                       y = ATP_targets_stauro$relAbundance,
                       groups = ATP_targets_stauro$uniqueID,
                       BPPARAM = BPPARAM,
                       seed = 123,
                       maxAttempts = 100,
                       alwaysPermute = FALSE)

  rss0_new <- sapply(models, function(m) {
    ifelse(inherits(m , "try-error"), NA, m$m$deviance())
  })

  expect_equal(unname(rss0_new)[-16], rss0_ref[-16]) # position 16: ATP5G1_IPI00009075 -> case for negative RSS-Diff before resampling?
})


test_that("fitSigmoidsParallel_allok_rss1", {

  models <- fitSigmoidsParallel(x = ATP_targets_stauro$temperature,
                                y = ATP_targets_stauro$relAbundance,
                                group = paste(ATP_targets_stauro$uniqueID, ATP_targets_stauro$compoundConcentration),
                                BPPARAM = BPPARAM,
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
