context("hypervolume")

test_that("Check hypervolume function working properly", {
  skip_if_missing_module()
  point <- matrix(c(0,0,0),nrow=3)
  reference <- (c(1,1,1))
  expect_equal( GetHypervolume(point,reference), 1)

  point <- matrix(c(0,0,1,1,0,0,0,1,0),nrow=3)
  reference <- (c(1.1,1.1,1.1))
  expect_equal( GetHypervolume(point,reference), 0.331)
})
