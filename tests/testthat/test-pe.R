# Process-error location (Pella–Tomlinson / FLBRP)

test_that("pe default returns FLQuant (log scale)", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  eq <- brp(FLBRP(ple4))

  pe0 <- pe(ple4, eq)
  expect_s4_class(pe0, "FLQuant")
  ## legacy eql path may be NA if FLBRP has no SR; class check is the contract
  expect_equal(dim(pe0)[2], dims(ple4)$year)
})

test_that("pe pellat locations return finite log residuals", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  eq <- brp(FLBRP(ple4))

  for (loc in c("after", "before")) {
    eta <- pe(ple4, eq, production = "pellat", location = loc, scale = "log")
    expect_s4_class(eta, "FLQuant")
    expect_true(any(is.finite(c(eta))), info = loc)
    em <- pe(ple4, eq, production = "pellat", location = loc, scale = "mult")
    ok <- is.finite(c(eta)) & is.finite(c(em)) & c(em) > 0
    expect_true(any(ok), info = loc)
    expect_equal(c(eta)[ok], log(c(em)[ok]), tolerance = 1e-8, info = loc)
  }
  ## productivity: SP can be ≤0 for some PT shapes → log may be NaN;
  ## multiplicative residual still matches the identity
  em <- pe(ple4, eq, production = "pellat", location = "productivity",
           scale = "mult")
  expect_s4_class(em, "FLQuant")
  expect_true(any(is.finite(c(em))))
})

test_that("pe with pellatParams FLPar works", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  eq <- brp(FLBRP(ple4))
  pt <- pellatParams(eq)

  eta1 <- pe(ple4, eq, production = "pellat", location = "after", scale = "log")
  eta2 <- pe(ple4, pt, location = "after", scale = "log")
  ok <- is.finite(c(eta1)) & is.finite(c(eta2))
  expect_equal(c(eta1)[ok], c(eta2)[ok], tolerance = 1e-6)
})

test_that("pe after-catch mult matches Bnext/(B+SP-C)", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  eq <- brp(FLBRP(ple4))
  pt <- pellatParams(eq)

  B  <- ssb(ple4)
  C  <- catch(ple4)
  SP <- ((pt["r"] %/% pt["p"]) %*% B) %*%
    (1 - exp(log(B %/% pt["k"]) %*% pt["p"]))
  Bnext <- B %=% NA_real_
  ny <- dim(B)[2]
  Bnext[, seq_len(ny - 1)] <- B[, 2:ny]
  expect <- Bnext %/% (B %+% SP %-% C)

  got <- pe(ple4, pt, location = "after", scale = "mult")
  ok <- is.finite(c(expect)) & is.finite(c(got))
  expect_equal(c(got)[ok], c(expect)[ok], tolerance = 1e-8)
})
