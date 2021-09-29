test_that("simulation returns a data.frame", {
  result.data <- driftSim:::runSimulation(
    N=20
    ,v=1.0
    ,c=2.0
    ,is_pure=T
    ,mu=0.01
    ,max_time=100
    ,pHawk_init=0.5
    ,output_nth_generation=5
    ,sd_pHawkMixed = 0.0
  )
  print(paste("nu dan: ", result.data[1,"freq_Hawk"]))
  expect_type(result.data,"list")
})

test_that("there are nonzero hawks",{
  result.data <- driftSim:::runSimulation(
    N=500
    ,v=1.0
    ,c=2.0
    ,is_pure=T
    ,mu=0.01
    ,max_time=10
    ,pHawk_init=1.0
    ,output_nth_generation=1
    ,sd_pHawkMixed = 0.0
  )

  # check whether values are larger or equal than 0
  expect_true(result.data[1,"freq_Hawk"] >= 0.5)

  # check whether everything is below 1.0
  expect_false(F %in% (result.data[,"freq_Hawk"] <= 1.0))
})

test_that("continous evolution results in nonzero standard dev",{
  result.data <- driftSim:::runSimulation(
    N=500
    ,v=1.0
    ,c=2.0
    ,is_pure=F
    ,mu=0.01
    ,max_time=100
    ,pHawk_init=1.0
    ,output_nth_generation=1
    ,sd_pHawkMixed = 0.1
  )

  # check whether everything is below 1.0
  expect_false(F %in% (result.data[,"sd_pHawkMixed"] > 0.0))
})

test_that("lack of continuous evolution: 0 dev",{
  result.data <- driftSim:::runSimulation(
    N=500
    ,v=1.0
    ,c=2.0
    ,is_pure=T
    ,mu=0.01
    ,max_time=100
    ,pHawk_init=1.0
    ,output_nth_generation=1
    ,sd_pHawkMixed = 0.1
  )

  # check whether everything is below 1.0
  expect_false(F %in% (result.data[,"sd_pHawkMixed"] == 0.0))
})
