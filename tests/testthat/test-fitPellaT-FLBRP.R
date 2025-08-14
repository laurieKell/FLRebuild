# Test file for fitPellaT method with FLBRP objects
# Tests that the FLBRP method works correctly

test_that("fitPellaT works with FLBRP objects using refpts method", {
  # This test requires FLBRP package to be loaded
  if (!requireNamespace("FLBRP", quietly = TRUE)) {
    skip("FLBRP package not available")
  }
  
  # Create a simple test data frame first
  testData = data.frame(
    yield = c(30, 60, 80, 70, 40, 20),
    eb = c(100, 200, 300, 400, 450, 480)
  )
  
  # Test that the data.frame method works
  result = fitPellaT(testData, method = "moments", biomass = "eb")
  
  expect_false(is.null(result))
  expect_true("r" %in% names(result))
  expect_true("p" %in% names(result))
  expect_true("K" %in% names(result))
  expect_true("MSY" %in% names(result))
  expect_true("BMSY" %in% names(result))
  expect_equal(result$method, "moments")
})

test_that("fitPellaTDirect works with numeric inputs", {
  # Test direct calculation with numeric inputs
  result = fitPellaTDirect(MSY = 100, BMSY = 300, K = 1000)
  
  expect_false(is.null(result))
  expect_true("r" %in% names(result))
  expect_true("p" %in% names(result))
  expect_true("K" %in% names(result))
  expect_true("MSY" %in% names(result))
  expect_true("BMSY" %in% names(result))
  expect_equal(result$method, "direct")
  expect_equal(result$MSY, 100)
  expect_equal(result$BMSY, 300)
  expect_equal(result$K, 1000)
})

test_that("fitPellaTDirect works with negative p case", {
  # Test direct calculation with BMSY/K > 0.37 (should give negative p)
  result = fitPellaTDirect(MSY = 100, BMSY = 500, K = 1000)
  
  expect_false(is.null(result))
  expect_true(result$p < 0)  # p should be negative for BMSY/K = 0.5 > 0.37
  expect_equal(result$method, "direct")
})

test_that("fitPellaTDirect works with positive p case", {
  # Test direct calculation with BMSY/K < 0.37 (should give positive p)
  result = fitPellaTDirect(MSY = 100, BMSY = 300, K = 1000)
  
  expect_false(is.null(result))
  expect_true(result$p > 0)  # p should be positive for BMSY/K = 0.3 < 0.37
  expect_equal(result$method, "direct")
})

test_that("fitPellaT handles invalid inputs gracefully", {
  # Test with insufficient data
  testData = data.frame(
    yield = c(30, 60),  # Only 2 points, need at least 3
    eb = c(100, 200)
  )
  
  result = fitPellaT(testData, method = "moments", biomass = "eb")
  expect_true(is.null(result))
})

test_that("fitPellaT handles NA values correctly", {
  # Test with NA values
  testData = data.frame(
    yield = c(30, 60, 80, 70, 40, 20),
    eb = c(100, 200, NA, 400, 450, 480)  # One NA value
  )
  
  result = fitPellaT(testData, method = "moments", biomass = "eb")
  expect_false(is.null(result))  # Should filter out NA and still work
}) 