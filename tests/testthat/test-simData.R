library(PLSDAbatch)

context('simData')
test_that("data simulation works", {

  simu_data <- simData(bat.mean = 7,
                       sd.bat = 8,
                       trt.mean = 3,
                       sd.trt = 2,
                       noise = 0.2,
                       N = 40,
                       p_total = 300,
                       p_trt_relevant = 60,
                       p_bat_relevant = 150,
                       percentage_samples = 0.5,
                       percentage_overlap_variables = 0.5)

  X <- simu_data$data
  trt <- simu_data$Y.trt
  batch <- simu_data$Y.bat
  true.trt <- simu_data$true.trt
  true.batch <- simu_data$true.batch

  expect_is(X, 'matrix')
  expect_is(trt, 'character')
  expect_is(batch, 'character')
  expect_is(true.trt, 'integer')
  expect_is(true.batch, 'integer')
})
