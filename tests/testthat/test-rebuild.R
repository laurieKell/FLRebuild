test_that("rebuild generic function works", {
  # Test that generic function exists
  expect_true(isGeneric("rebuild"))
})

test_that("rebuildTime generic function works", {
  # Test that generic function exists
  expect_true(isGeneric("rebuildTime"))
})

test_that("calculateRecoveryTime works correctly", {
  # Test basic functionality
  expect_equal(calculateRecoveryTime(0.5, 1.0, 0.1), log(2)/0.1)
  
  # Test when already at target
  expect_equal(calculateRecoveryTime(1.0, 1.0, 0.1), 0)
  
  # Test when above target
  expect_equal(calculateRecoveryTime(1.5, 1.0, 0.1), 0)
})

test_that("calculateTmax works with generation method", {
  # Mock object with generation time
  mock_object <- list()
  attr(mock_object, "generation_time") <- 10
  
  # Test generation method
  result <- calculateTmax(mock_object, method = "generation", tmin = 5)
  expect_equal(result, 15)
})

test_that("calculateTmax works with multiply method", {
  # Test multiply method
  result <- calculateTmax(NULL, method = "multiply", tmin = 5)
  expect_equal(result, 10)
})

test_that("calculateTmax throws error for invalid method", {
  expect_error(calculateTmax(NULL, method = "invalid"))
})

test_that("calculateTmax throws error for missing parameters", {
  expect_error(calculateTmax(NULL, method = "generation"))
  expect_error(calculateTmax(NULL, method = "multiply"))
})

test_that("abiAge, abiMsy, and abi methods exist and can be called", {
  expect_true(isGeneric("abiAge"))
  expect_true(isGeneric("abiMsy"))
  expect_true(isGeneric("abi"))
})

test_that("blim method exists and can be called", {
  expect_true(isGeneric("blim"))
})

# The following tests require FLBRP and FLStock objects, so we check for their existence
if (requireNamespace("FLBRP", quietly = TRUE) && requireNamespace("FLCore", quietly = TRUE)) {
  test_that("abiAge and abiMsy return FLQuant for FLBRP", {
    data(ple4brp, package = "FLBRP")
    expect_s4_class(abiAge(ple4brp), "FLQuant")
    expect_s4_class(abiMsy(ple4brp), "FLQuant")
  })

  test_that("abi returns FLQuant for FLStock and FLBRP", {
    data(ple4, package = "FLCore")
    data(ple4brp, package = "FLBRP")
    expect_s4_class(abi(ple4, ple4brp), "FLQuant")
  })

  test_that("blim returns FLPar for FLBRP", {
    data(ple4brp, package = "FLBRP")
    expect_s4_class(blim(ple4brp), "FLPar")
  })
} else {
  test_that("FLBRP and FLCore not available, skipping FLBRP/FLStock tests", {
    skip("FLBRP or FLCore not available")
  })
} 