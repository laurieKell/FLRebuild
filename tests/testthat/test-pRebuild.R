test_that("pRebuild reports rebuild times for an FLBRP and a list for FLBRPs", {
  library(FLCore)
  library(FLBRP)
  data(ple4)

  eq <- brp(FLBRP(ple4))
  rf <- pRebuild(eq, pRec = 0.20, maxYears = 80)

  expect_s4_class(rf, "FLPar")
  expect_true(all(c("B0", "Bmsy", "Fmsy", "GT", "SSB_pR",
                    "Trec_F0", "Trec_Fmsy", "Trec_HCR",
                    "Trec_F0_GT", "Bstar_F0", "Bstar_Fmsy") %in%
                  dimnames(rf)$params))
  expect_gt(c(rf["B0"]), c(rf["Bmsy"]))
  expect_gt(c(rf["Fmsy"]), 0)
  expect_gt(c(rf["GT"]), 1)
  expect_true(is.finite(c(rf["Trec_F0"])))
  expect_true(is.finite(c(rf["Trec_Fmsy"])))
  expect_true(is.finite(c(rf["Trec_HCR"])))
  expect_gte(c(rf["Trec_F0"]), 0)
  expect_lte(c(rf["Trec_F0"]), c(rf["Trec_HCR"]))
  expect_lte(c(rf["Trec_HCR"]), c(rf["Trec_Fmsy"]))
  expect_true(is.finite(c(rf["Bstar_F0"])))
  expect_true(is.finite(c(rf["Bstar_Fmsy"])))
  expect_gt(c(rf["Bstar_F0"]), 0)
  expect_lt(c(rf["Bstar_Fmsy"]), c(rf["Bmsy"]))
  expect_lte(c(rf["Bstar_F0"]), c(rf["Bstar_Fmsy"]))

  brps <- FLBRPs(plaice = eq, plaice2 = eq)
  rfs <- pRebuild(brps, pRec = 0.20, maxYears = 80)
  expect_type(rfs, "list")
  expect_named(rfs, c("plaice", "plaice2"))
  expect_s4_class(rfs[["plaice"]], "FLPar")

  rf2 <- plyr::ldply(rfs, model.frame)
  expect_true(all(c(".id", "Bmsy", "Trec_HCR_GT", "Bstar_F0") %in% names(rf2)))
  expect_equal(nrow(rf2), 2L)
})

test_that("SSB_pR is the Beverton-Holt biomass at pRec of virgin recruitment", {
  library(FLCore)
  library(FLBRP)
  data(ple4)

  eq <- FLBRP(ple4)
  model(eq) <- bevholt()$model
  params(eq) <- FLPar(a = 2, b = 1)
  eq <- brp(eq)
  rf <- pRebuild(eq, pRec = 0.20, maxYears = 80)

  B0 <- c(rf["B0"])
  b <- 1
  p <- 0.20
  expect_equal(c(rf["SSB_pR"]), p * b * B0 / (b + B0 * (1 - p)))
  expect_lt(c(rf["SSB_pR"]), 0.2 * B0)
  expect_lte(c(rf["Trec_F0"]), c(rf["Trec_HCR"]))
  expect_lte(c(rf["Trec_HCR"]), c(rf["Trec_Fmsy"]))
})
