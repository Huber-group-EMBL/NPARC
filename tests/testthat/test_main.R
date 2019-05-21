data(staurosporineTPP)

test_that("nparc_allok", {
  res <- nparc(x = staurosporineTPP$temperature,
               y = staurosporineTPP$relAbundance,
               group = staurosporineTPP$compoundConcentration,
               id = staurosporineTPP$uniqueID)
})
