# Test file for ageMax functions
context("ageMax functions")

test_that("ageMax works with simple FLQuant", {
  # Create a simple FLQuant for testing
  flq = FLQuant(c(1, 3, 2, 4, 5), 
                dimnames = list(age = 1:5, year = 2000, iter = 1))
  
  result = ageMax(flq)
  
  # Should return age index 4 (highest value excluding plusgroup)
  expect_equal(c(result), 4)
  
  # Should have 1 quant dimension
  expect_equal(dim(result)[1], 1)
})

test_that("ageMax works with multiple years and iterations", {
  # Create FLQuant with multiple years and iterations
  flq = FLQuant(c(1, 3, 2, 4, 5,  # year 1
                  5, 2, 4, 1, 3,  # year 2
                  2, 5, 1, 3, 4), # year 3
                dimnames = list(age = 1:5, year = 2000:2002, iter = 1))
  
  result = ageMax(flq)
  
  # Check results for each year (excluding plusgroup)
  expect_equal(c(result[,1]), 4)  # age 4 has max value in year 1 (excluding age 5)
  expect_equal(c(result[,2]), 1)  # age 1 has max value in year 2 (excluding age 5)
  expect_equal(c(result[,3]), 2)  # age 2 has max value in year 3 (excluding age 5)
  
  # Should have 1 quant dimension
  expect_equal(dim(result)[1], 1)
})

test_that("ageMax handles NA values correctly", {
  # Create FLQuant with NA values
  flq = FLQuant(c(1, NA, 2, 4, 5), 
                dimnames = list(age = 1:5, year = 2000, iter = 1))
  
  result = ageMax(flq, na.rm = TRUE)
  
  # Should return age index 4 (highest non-NA value, excluding plusgroup)
  expect_equal(c(result), 4)
})

test_that("ageMaxQ works similarly to ageMax", {
  # Create a simple FLQuant for testing
  flq = FLQuant(c(1, 3, 2, 4, 5), 
                dimnames = list(age = 1:5, year = 2000, iter = 1))
  
  result1 = ageMax(flq)
  result2 = ageMaxQ(flq)
  
  # Both should return the same result
  expect_equal(result1, result2)
})

test_that("ageMax works with FLQuant that has age dimension", {
  # Test specifically for the case that was causing the error
  # Create FLQuant with age dimension (like ssb.age output)
  flq = FLQuant(c(1, 3, 2, 4, 5), 
                dimnames = list(age = 1:5, year = 2000, iter = 1))
  
  # This should not throw an error about conflicting dimension names
  expect_no_error(result <- ageMax(flq))
  
  # Should return a valid result
  expect_true(is(result, "FLQuant"))
  expect_equal(dim(result)[1], 1)  # 1 quant dimension
  expect_equal(dim(result)[2], 1)  # 1 year
  expect_equal(dim(result)[6], 1)  # 1 iteration
}) 