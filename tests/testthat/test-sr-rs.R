test_that("invSRR generic and methods exist", {
  expect_true(isGeneric("invSRR"))
})

test_that("rmax and rmsy generics exist", {
  expect_true(isGeneric("rmax"))
  expect_true(isGeneric("rmsy"))
})

test_that("refCreate generic exists", {
  expect_true(isGeneric("refCreate"))
})

if (requireNamespace("FLBRP", quietly = TRUE) && requireNamespace("FLCore", quietly = TRUE)) {
  test_that("rmax and rmsy return expected types for FLBRP", {
    data(ple4brp, package = "FLBRP")
    expect_true(is.numeric(rmax(ple4brp)) || is( rmax(ple4brp), "FLPar"))
    expect_true(is.numeric(rmsy(ple4brp)) || is( rmsy(ple4brp), "FLPar"))
  })

  test_that("rmax and rmsy with ratio argument work for FLBRP", {
    data(ple4brp, package = "FLBRP")
    out_rmax <- rmax(ple4brp, ratio = 0.5)
    out_rmsy <- rmsy(ple4brp, ratio = 0.5)
    expect_s4_class(out_rmax, "FLPar")
    expect_s4_class(out_rmsy, "FLPar")
    expect_true("rmax" %in% dimnames(out_rmax)$refpt)
    expect_true("rmsy" %in% dimnames(out_rmsy)$refpt)
  })

  test_that("refCreate returns FLPar for FLBRP and does not modify input", {
    data(ple4brp, package = "FLBRP")
    original_refpts <- refpts(ple4brp)
    out <- refCreate(ple4brp, "rec", 1000)
    expect_s4_class(out, "FLPar")
    expect_true("rec" %in% dimnames(out)$refpt)
    # Check that ple4brp is unchanged
    expect_identical(refpts(ple4brp), original_refpts)
  })

  test_that("refCreate returns FLPar for FLPar and does not modify input", {
    params <- FLPar(a = 1000, b = 2)
    original_params <- params
    out <- refCreate(params)
    expect_s4_class(out, "FLPar")
    # Check that params is unchanged
    expect_identical(params, original_params)
  })

  test_that("invSRR returns FLQuant for FLBRP", {
    data(ple4brp, package = "FLBRP")
    rec <- FLQuant(1000)
    out <- invSRR(ple4brp, rec)
    expect_s4_class(out, "FLQuant")
  })

  test_that("refptsEB returns FLPar with eb quant for FLBRP", {
    data(ple4brp, package = "FLBRP")
    out <- refptsEB(ple4brp)
    expect_s4_class(out, "FLPar")
    expect_true("eb" %in% dimnames(out)$quant)
    # Check that all original quants plus 'eb' are present
    orig_quants <- c("harvest", "yield", "rec", "ssb", "biomass", "revenue", "cost", "profit")
    expect_true(all(orig_quants %in% dimnames(out)$quant))
    expect_true("eb" %in% dimnames(out)$quant)
  })
} else {
  test_that("FLBRP or FLCore not available, skipping sr-rs tests", {
    skip("FLBRP or FLCore not available")
  })
}

test_that("bevholtInv and rickerInv return numeric or FLQuant", {
  params <- FLPar(a = 1000, b = 2)
  rec <- 500
  expect_true(is.numeric(bevholtInv(params, rec)) || is(bevholtInv(params, rec), "FLQuant"))
  expect_true(is.numeric(rickerInv(params, rec)) || is(rickerInv(params, rec), "FLQuant"))
}) 