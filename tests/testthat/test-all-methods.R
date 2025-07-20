# Test all S4 generics and methods

test_that("All S4 generics exist", {
  # Core rebuilding functions
  expect_true(isGeneric("rebuild"))
  expect_true(isGeneric("rebuildTime"))
  
  # Reference point functions
  expect_true(isGeneric("blim"))
  expect_true(isGeneric("msyVirgin"))
  
  # Age-based indicator functions
  expect_true(isGeneric("abiAge"))
  expect_true(isGeneric("abi"))
  expect_true(isGeneric("abiMsy"))
  
  # Stock-recruitment functions
  expect_true(isGeneric("invSRR"))
  expect_true(isGeneric("refCreate"))
  expect_true(isGeneric("rmax"))
  expect_true(isGeneric("rmsy"))
  expect_true(isGeneric("rvirgin"))
  
  # Utility functions
  expect_true(isGeneric("years"))
  expect_true(isGeneric("rboot"))
  expect_true(isGeneric("refptsEB"))
  
  # Pella-Tomlinson functions
  expect_true(isGeneric("rTime"))
  expect_true(isGeneric("brebuild"))
})

test_that("All exported functions exist", {
  # Regular functions
  expect_true(exists("calculateRecoveryTime"))
  expect_true(exists("calculateTmax"))
  expect_true(exists("gtLife"))
  expect_true(exists("gtR0"))
  expect_true(exists("plotRebuildTrajectories"))
  expect_true(exists("abistock"))
})

# Test with FLBRP objects if available
if (requireNamespace("FLBRP", quietly = TRUE) && requireNamespace("FLCore", quietly = TRUE)) {
  test_that("FLBRP methods work", {
    data(ple4brp, package = "FLBRP")
    
    # Test blim
    expect_s4_class(blim(ple4brp), "FLPar")
    
    # Test msyVirgin
    expect_true(is.numeric(msyVirgin(ple4brp)))
    
    # Test refCreate
    expect_s4_class(refCreate(ple4brp, "test", 1000), "FLPar")
    
    # Test rmax
    expect_true(is.numeric(rmax(ple4brp)) || is(rmax(ple4brp), "FLPar"))
    expect_s4_class(rmax(ple4brp, 0.5), "FLPar")
    
    # Test rmsy
    expect_s4_class(rmsy(ple4brp), "FLPar")
    expect_s4_class(rmsy(ple4brp, 0.5), "FLPar")
    
    # Test rvirgin
    expect_s4_class(rvirgin(ple4brp), "FLPar")
    expect_s4_class(rvirgin(ple4brp, 0.5), "FLPar")
    
    # Test invSRR
    rec <- FLQuant(1000)
    expect_s4_class(invSRR(ple4brp, rec), "FLQuant")
    
    # Test refptsEB
    expect_s4_class(refptsEB(ple4brp), "FLPar")
    
    # Test rebuild
    expect_s4_class(rebuild(ple4brp), "FLStock")
    
    # Test abiAge
    expect_s4_class(abiAge(ple4brp), "FLQuant")
    
    # Test abiMsy
    expect_s4_class(abiMsy(ple4brp), "FLQuant")
  })
}

# Test with biodyn objects if available
if (requireNamespace("biodyn", quietly = TRUE)) {
  test_that("biodyn methods work", {
    bd <- biodyn::biodyn(FLPar(r = 0.5, k = 1000, p = 1))
    
    # Test rebuild
    expect_s4_class(rebuild(bd), "biodyn")
    
    # Test rebuildTime
    expect_true(is.data.frame(rebuildTime(bd)))
  })
}

# Test with FLStock objects if available
if (requireNamespace("FLCore", quietly = TRUE)) {
  test_that("FLStock methods work", {
    data(ple4, package = "FLCore")
    
    # Test rebuildTime
    expect_true(is.data.frame(rebuildTime(ple4)))
  })
}

# Test utility functions
test_that("Utility functions work", {
  # Test calculateRecoveryTime
  expect_equal(calculateRecoveryTime(0.5, 1.0, 0.1), log(2)/0.1)
  
  # Test calculateTmax
  expect_equal(calculateTmax(NULL, method = "multiply", tmin = 5), 10)
  
  # Test gtLife
  ages <- 0:5
  lx <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
  mx <- c(0, 0, 2, 3, 1, 0)
  expect_true(is.numeric(gtLife(ages, lx, mx)))
  
  # Test gtR0
  expect_true(is.numeric(gtR0(2, 0.2)))
})

# Test Pella-Tomlinson functions
test_that("Pella-Tomlinson functions work", {
  # Test rTime
  expect_true(is.numeric(rTime(0.3, r = 0.2, p = 1.5)))
  
  # Test brebuild
  expect_true(is.numeric(brebuild(3, r = 0.2, p = 1.5)))
})

# Test FLQuant methods
if (requireNamespace("FLCore", quietly = TRUE)) {
  test_that("FLQuant methods work", {
    flq <- FLQuant(1:10, dimnames = list(year = 2000:2009))
    
    # Test years
    expect_s4_class(years(flq), "FLQuant")
    
    # Test rboot
    expect_s4_class(rboot(flq, n = 5), "FLQuant")
  })
} 