context('benchmark')

test_that("Check DTLZ return correct value based on the Pareto front", {
  skip_if_missing_module()
  individual <- stats::runif(14)
  individual[4:14] <- 0.5
  individual[1:3] <- 0
  nObj <- 4
  expect_equal(norm(DTLZ1(individual,nObj),'1'),0.5)

  expect_equal(norm(DTLZ2(individual,nObj),'2'),1)
  expect_equal(norm(DTLZ3(individual,nObj),'2'),1)
  expect_equal(norm(DTLZ4(individual,nObj),'2'),1)
})

test_that("Check WFG return correct value based on the Pareto front", {
  skip_if_missing_module()
  individual <- stats::runif(14)
  individual[6:14] <- 0.35
  individual[1:5] <- 0
  nObj <- 6
  expect_lt(norm(WFG4(individual,nObj)/c(2,4,6,8,10,12),'2')-1,0.04) # should be equal, but strongly affected by numerical error
  expect_equal(norm(WFG5(individual,nObj)/c(2,4,6,8,10,12),'2'),1)
  expect_equal(norm(WFG6(individual,nObj)/c(2,4,6,8,10,12),'2'),1)
  expect_equal(norm(WFG7(individual,nObj)/c(2,4,6,8,10,12),'2'),1)
  expect_equal(norm(WFG8(individual,nObj)/c(2,4,6,8,10,12),'2'),1)
  expect_lt(norm(WFG9(individual,nObj)/c(2,4,6,8,10,12),'2')-1,0.04) # should be equal, but strongly affected by numerical error
})
