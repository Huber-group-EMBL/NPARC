data(staurosporineTPP)
library(dplyr)
library(magrittr)
library(NPARC)
library(testthat)

ATP_targets_stauro <- filter(staurosporineTPP, grepl("ATP|STK4_", uniqueID))
rss0_ref <- c( 0.042328688, 0.574210782, 1.105606894, 66.141814549, 1.279980245, 1.445979905, NA, 0.085035397, NA, 0.086908456, 0.040746912, 0.229219646, 0.185040763, 0.069901302, 3.219166488, 0.890621716, 1.482869417, 0.111385691, 0.890474276, 0.096161260, 0.005317976, 0.165361410, 0.261380998, 2.004720774, NA, 1.017639038, 0.037805608, 0.079748130, 0.080865091, 0.023917017, 0.302070404, 0.438928745, 0.021234921, 0.207783991, 0.511561035, NA, 0.047061699, 0.022204837, 0.018557014, 82.669418181, 0.171696257, 0.246755638, 0.027504901, 0.042566445, 0.685432518, 1.2181315672231198821)
rss1_ref <- c( NA, 0.52452012, 1.02728630, NA, 1.16253040, NA, NA, NA, NA, 0.08085710, 0.03391284, 0.19178616, 0.11715310, 0.05738689, 2.85347619, NA, 0.88566800, NA, 0.86233078, 0.05913776, NA, NA, 0.24112396, NA, NA, NA, 0.03726856, 0.07672426, 0.07933303, 0.02046256, NA, 0.42376129, 0.01929703, 0.19330061, 0.49355960, NA, 0.03926680, 0.01932926, NA, 82.39894648, 0.14698662, 0.22814391, 0.02189395, 0.03691358, 0.57373856, 0.083147452021042442)
rssDiff_ref <- c(NA, 0.0496906604, 0.0783205890, NA, 0.1174498426, NA, NA, NA, NA, 0.0060513570, 0.0068340678, 0.0374334892, 0.0678876612, 0.0125144109, 0.3656902966, NA, 0.5972014163, NA, 0.0281434961, 0.0370234978, NA, NA, 0.0202570335, NA, NA, NA, 0.0005370464, 0.0030238719, 0.0015320641, 0.0034544609, NA, 0.0151674541, 0.0019378859, 0.0144833862, 0.0180014361, NA, 0.0077948980, 0.0028755747, NA, 0.2704717055, 0.0247096421, 0.0186117273, 0.0056109477, 0.0056528626, 0.1116939612, 1.13498411520207742598)
n0_ref <- c(10, 40, 40, 20, 20, 20, 0, 20, 0, 40, 40, 40, 40, 40, 20, 10, 20, 10, 40, 30, 10, 10, 20, 20, 0, 10, 30, 40, 40, 40, 10, 40, 40, 30, 40, 0, 40, 30, 10, 20, 40, 40, 40, 40, 40, 40)
n1_ref <- c(10, 40, 40, 20, 20, 20, 0, 20, 0, 40, 40, 40, 40, 40, 20, 10, 20, 10, 40, 30, 10, 10, 20, 20, 0, 10, 30, 40, 40, 40, 10, 40, 40, 30, 40, 0, 40, 30, 10, 20, 40, 40, 40, 40, 40, 40)
repeats_ref <- c(11, 0, 0, 11, 0, 11, 11, 11, 11, 0, 0, 0, 0, 0, 0, 11, 0, 11, 0, 0, 11, 11, 0, 11, 11, 11, 0, 0, 0, 0, 11, 0, 0, 0, 0, 11, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0)
tm_0_ref <- c(NA, 61.733018503993208, NaN, NA, NaN, NaN, NA, NA, NA, 45.507090381927227, 48.370656179268586, 46.177862386428778, NaN, 50.747974622093842, 51.108100718477552, NA, 53.716668890572350, NA, 49.944382649054042, 53.268268043773986, NA, 64.083619744519751, 52.630782071773126, 54.403632566800773, NA, NA, 49.973318102061327, 54.987190448048963, 55.650774547079479, 53.776217585739715, NA, 63.931766238670086, 53.474065831051554, 59.660513721432004, 53.651458551339147, NA, 53.400317830126923, 56.166341210006678, NA, 49.957133804531374, 58.379308939431617, 61.739135919243360, 49.739613672120782, 51.424194381078600, NaN, 51.075584610286953)
tm_20_ref <- c(51.826211134785211, 66.548397028205514, 77.195149424582979, 55.304171825879109, 54.826147979730159, NA, NA, 56.362667583506500, NA, 45.626792065465864, 48.448244898573385, 47.460807063834899, 60.921156578436033, 50.709466817179830, NaN,                NaN,                NaN, 53.705720351822045, 50.199251315345073, 51.808623346193485, 54.052321351744041, NA, 51.186531862923864, NA, NA,                NaN, 50.030558335114812, 54.702758421926561, 55.813688161147894, 53.533856781624372, 49.760889873698368, 63.528237092905130, 53.782382421569636, 60.003138034648828, 53.162997375385423, NA, 53.292278053288818, 56.235749170195618, 53.673086380853192, 55.280037025265500, 59.676701562973001, 60.694252715919731, 49.732243787451694, 51.158615360263845, 70.641541915911802, 59.740907732406420)
pAdj_ref <- c(NA, 0.699795397624508708, 0.753523752550566694, NA, 0.861339846144721699, NA, NA, NA, NA, 0.753523752550566694, 0.315762492309564990, 0.315762492309564990, 0.019064553849503785, 0.315762492309564990, 0.861339846144721699, NA, 0.315762492309564990, NA, 0.861339846144721699, 0.077448008418172787, NA, NA, 0.861339846144721699, NA, NA, NA, 0.983055252694920556, 0.861339846144721699, 0.945666731555014795, 0.368321920466761210, NA, 0.861339846144721699, 0.695132899776564450, 0.861339846144721699, 0.861339846144721699, NA, 0.315762492309564990, 0.695132899776564450, NA, 0.997288089402401878, 0.368321920466761210, 0.753523752550566694, 0.315762492309564990, 0.411224424237822794, 0.315762492309564990, 0.000000000000000000)


test_that("invokeRSSdiff_allok", {

  rssDiffs <- invokeRSSdiff(x = ATP_targets_stauro$temperature,
                            y = ATP_targets_stauro$relAbundance,
                            group = ATP_targets_stauro$compoundConcentration,
                            id = ATP_targets_stauro$uniqueID,
                            BPPARAM = BiocParallel::SnowParam(4, progressbar = TRUE))

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

  pAdj <- nparFtest(rss0=rss0_ref, rss1=rss1_ref, df_type = "theoretical", n0=n0_ref, n1=n1_ref, pars0=pars0, pars1=pars1)

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

  pAdj <- nparc(x = ATP_targets_stauro$temperature,
                y = ATP_targets_stauro$relAbundance,
                group = ATP_targets_stauro$compoundConcentration,
                id = ATP_targets_stauro$uniqueID,
                BPPARAM = BiocParallel::SnowParam(4, progressbar = TRUE))

  expect_equal(pAdj, pAdj_ref, tolerance = 1e-6)
})


test_that("nparc_allok_largedata", {

  system.time({
    pAdj <- nparc(x = staurosporineTPP$temperature,
                  y = staurosporineTPP$relAbundance,
                  group = staurosporineTPP$compoundConcentration,
                  id = staurosporineTPP$uniqueID,
                  BPPARAM = BiocParallel::SnowParam(4, progressbar = TRUE))

  })

  expect_equal(pAdj, pAdj_ref, tolerance = 1e-6)
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
