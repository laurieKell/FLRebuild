# Closed-form Pella-Tomlinson shape / Trebuild (no FLStock required)

test_that("shape2p and p2shape are inverses at Fox and Schaefer", {
  expect_equal(p2shape(shape2p(fox_shape())), fox_shape(), tolerance = 1e-8)
  expect_equal(p2shape(shape2p(schaefer_shape())), schaefer_shape(),
               tolerance = 1e-8)
  expect_equal(p2shape(1), 0.5, tolerance = 1e-10)
  expect_equal(p2shape(0), exp(-1), tolerance = 1e-10)
})

test_that("impose_shape keeps BMSY and sets K = BMSY/shape", {
  imp <- impose_shape(fmsy = 0.2, bmsy = 100, shape = schaefer_shape())
  expect_equal(imp$k, 200)
  expect_equal(imp$p, 1, tolerance = 1e-8)
  expect_equal(imp$r, 0.4, tolerance = 1e-8)
})

test_that("Schaefer Trebuild is 7.6% slower than Fox from 0.30 BMSY", {
  fmsy <- 0.25
  fox <- impose_shape(fmsy, bmsy = 1, shape = fox_shape())
  sch <- impose_shape(fmsy, bmsy = 1, shape = schaefer_shape())
  t_fox <- trebuild_from_bbmsy(fox$r, fox$p, fox$shape, bbmsy0 = 0.30)
  t_sch <- trebuild_from_bbmsy(sch$r, sch$p, sch$shape, bbmsy0 = 0.30)
  expect_equal(as.numeric(t_sch / t_fox), 1.076, tolerance = 5e-3)
})

test_that("Euler and closed-form Trebuild agree", {
  imp <- impose_shape(0.2, 100, fox_shape())
  B1 <- 0.30 * imp$bmsy
  tc <- trebuild_pt(imp$r, imp$p, imp$k, B1, imp$bmsy)
  te <- trebuild_pt_euler(imp$r, imp$p, imp$k, B1, imp$bmsy, dt = 0.01)
  expect_equal(te, tc, tolerance = 0.05)
})

test_that("spict_n_from_shape maps Schaefer to n = 2", {
  expect_equal(spict_n_from_shape(0.5), 2)
  expect_equal(spict_n_from_shape(fox_shape()), shape2p(fox_shape()) + 1,
               tolerance = 1e-8)
})
