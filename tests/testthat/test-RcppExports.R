test_that("simulation returns a data.frame", {
  result.data <- driftSim:::runSimulation(
    N=20
    ,v=1.0
    ,c=0.5
    ,mu=0.01
    ,max_time=100
    ,pHawk_init=0.5
    ,output_nth_generation=5
  )
  print(paste("nu dan: ", result.data[1,"pHawk"]))
  expect_type(result.data,"list")
})

test_that("there are nonzero hawks",{
  result.data <- driftSim:::runSimulation(
    N=500
    ,v=1.0
    ,c=0.5
    ,mu=0.01
    ,max_time=10
    ,pHawk_init=1.0
    ,output_nth_generation=1
  )

  # check whether values are larger or equal than 0
  expect_true(result.data[1,"pHawk"] >= 0.5)

  # check whether everything is below 1.0
  expect_false(F %in% (result.data[,"pHawk"] <= 1.0))
})
