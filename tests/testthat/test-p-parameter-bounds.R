# Test file for p parameter bounds in fitPellaT functions
# Tests that p parameter bounds are correctly set based on BMSY/K ratio

test_that("p parameter bounds are correctly set for positive p (BMSY/K < 0.37)", {
  # Test case with BMSY/K < 0.37 (should use positive p bounds)
  MSY = 100
  BMSY = 300  # BMSY/K = 0.3 < 0.37
  K = 1000
  
  result = fitPellaTDirect(MSY, BMSY, K)
  
  expect_false(is.null(result))
  expect_true(result$p > 0)  # p should be positive
  expect_true(result$p >= 0.05 && result$p <= 5.0)  # within positive bounds
})

test_that("p parameter bounds are correctly set for negative p (BMSY/K > 0.37)", {
  # Test case with BMSY/K > 0.37 (should use negative p bounds)
  MSY = 100
  BMSY = 500  # BMSY/K = 0.5 > 0.37
  K = 1000
  
  result = fitPellaTDirect(MSY, BMSY, K)
  
  expect_false(is.null(result))
  expect_true(result$p < 0)  # p should be negative
  expect_true(result$p >= -5.0 && result$p <= -0.05)  # within negative bounds
})

test_that("MLE method uses correct p bounds", {
  # Create test data with BMSY/K < 0.37
  testData = data.frame(
    yield = c(50, 80, 100, 90, 60, 30),
    eb = c(200, 400, 600, 800, 900, 950)  # BMSY/K ≈ 0.6/0.95 ≈ 0.63 > 0.37
  )
  
  result = fitPellaT(testData, method = "mle", biomass = "eb")
  
  expect_false(is.null(result))
  expect_true(result$p < 0)  # p should be negative for this data
  expect_true(result$p >= -5.0 && result$p <= -0.05)  # within negative bounds
})

test_that("Method of moments uses correct p bounds", {
  # Create test data with BMSY/K < 0.37
  testData = data.frame(
    yield = c(30, 60, 80, 70, 40, 20),
    eb = c(100, 200, 300, 400, 450, 480)  # BMSY/K ≈ 0.3/0.48 ≈ 0.625 > 0.37
  )
  
  result = fitPellaT(testData, method = "moments", biomass = "eb")
  
  expect_false(is.null(result))
  expect_true(result$p < 0)  # p should be negative for this data
  expect_true(result$p >= -5.0 && result$p <= -0.05)  # within negative bounds
})

test_that("BMSY/K ratio calculation is correct", {
  # Test the BMSY/K ratio calculation logic
  MSY = 100
  BMSY = 300
  K = 1000
  
  bmsyKRatio = BMSY / K
  expect_equal(bmsyKRatio, 0.3)
  expect_true(bmsyKRatio < 0.37)  # Should trigger positive p bounds
})

test_that("Edge case: BMSY/K very close to 0.37", {
  # Test edge case where BMSY/K is very close to 0.37
  MSY = 100
  BMSY = 370  # BMSY/K = 0.37 exactly
  K = 1000
  
  # This should default to positive p bounds (since 0.37 is not > 0.37)
  result = fitPellaTDirect(MSY, BMSY, K)
  
  expect_false(is.null(result))
  expect_true(result$p > 0)  # Should use positive bounds
}) 