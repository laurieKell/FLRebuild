# Pella-Tomlinson shape, r, and closed-form Trebuild (F = 0).
# Fox convention in this project: BMSY/K = 0.4 (not strictly 1/e).
# Schaefer: BMSY/K = 0.5 (p = 1).
#
# These replace the unexported one-liners in flbrp_pellat.R / OMstats.R
# when this file is loaded (later in alphabetical Collate order).

#' Fox production-function shape (BMSY/K)
#'
#' RebuildBiomass / Paper 2 convention: \code{0.4} (Pella-Tomlinson
#' \eqn{p \approx 0.19}), not the strict Fox limit \eqn{1/e}.
#'
#' @return Numeric scalar.
#' @export
#' @examples
#' fox_shape()
fox_shape <- function() 0.4

#' Schaefer production-function shape (BMSY/K)
#'
#' Symmetric Pella-Tomlinson: \code{BMSY/K = 0.5}, \eqn{p = 1}.
#'
#' @return Numeric scalar.
#' @export
#' @examples
#' schaefer_shape()
schaefer_shape <- function() 0.5

#' Shape BMSY/K from Pella-Tomlinson exponent p
#'
#' \eqn{(1/(1+p))^{1/p}}; Fox limit \eqn{p \to 0} is \eqn{1/e}.
#'
#' @param p numeric Pella-Tomlinson exponent (may be negative).
#' @return Numeric vector of shapes in (0, 1).
#' @export
#' @seealso \code{\link{shape2p}}
p2shape <- function(p) {
  p <- as.numeric(p)
  ifelse(abs(p) < 1e-10, exp(-1), (1 / (1 + p))^(1 / p))
}

#' Pella-Tomlinson exponent p from shape BMSY/K
#'
#' @param shape numeric BMSY/K in (0, 1).
#' @return Numeric vector of exponents \code{p}.
#' @export
#' @seealso \code{\link{p2shape}}
shape2p <- function(shape) {
  vapply(shape, function(s) {
    if (!is.finite(s) || s <= 0 || s >= 1) return(NA_real_)
    fn <- function(p) {
      if (abs(p) < 1e-10) exp(-1) - s else (1 / (1 + p))^(1 / p) - s
    }
    uniroot(fn, lower = -0.999, upper = 100, tol = 1e-10, extendInt = "upX")$root
  }, numeric(1))
}

#' Intrinsic rate r given FMSY and BMSY/K
#'
#' \eqn{r = F_{\mathrm{MSY}}\, p / (1 - \mathrm{shape}^p)}; Fox limit
#' \eqn{p \to 0} is \eqn{r = F_{\mathrm{MSY}} / (-\log(\mathrm{shape}))}.
#'
#' @param fmsy numeric harvest rate MSY/BMSY.
#' @param shape numeric BMSY/K.
#' @param p optional exponent; computed from \code{shape} if \code{NULL}.
#' @return Numeric vector of \code{r}.
#' @export
r_from_fmsy_shape <- function(fmsy, shape, p = NULL) {
  n <- length(fmsy)
  shape <- rep(shape, length.out = n)
  if (is.null(p)) p <- shape2p(shape)
  p <- rep(p, length.out = n)
  out <- fmsy * p / (1 - shape^p)
  tiny <- abs(p) < 1e-8
  out[tiny] <- fmsy[tiny] / (-log(shape[tiny]))
  out
}

#' Closed-form years from B1 to B2 under F = 0 (Pella-Tomlinson)
#'
#' \deqn{T = (1/r) \log[(z_2/(1-z_2)) \cdot ((1-z_1)/z_1)],\quad z=(B/K)^p}
#' Fox limit \eqn{p \to 0}: \eqn{T = (1/r) \log[\log(K/B_1)/\log(K/B_2)]}.
#' For \eqn{p<0} (left-skewed shape \eqn{< 1/e}), \eqn{z>1} for \eqn{B<K};
#' the same formula is used provided the argument of the log is positive.
#'
#' @param r intrinsic growth rate.
#' @param p Pella-Tomlinson exponent.
#' @param K carrying capacity.
#' @param B1 starting biomass.
#' @param B2 target biomass; default BMSY.
#' @param shape optional BMSY/K (used if \code{bmsy} is missing).
#' @param bmsy optional BMSY (default \code{K * shape} or \code{K * p2shape(p)}).
#' @return Numeric vector of years (\code{0} if already at/above target;
#'   \code{Inf} if the target is at K).
#' @export
#' @seealso \code{\link{trebuild_pt_euler}}, \code{\link{trebuild_from_bbmsy}}
trebuild_pt <- function(r, p, K, B1, B2 = NULL, shape = NULL, bmsy = NULL) {
  if (is.null(bmsy)) {
    bmsy <- if (is.null(shape)) K * p2shape(p) else K * shape
  }
  if (is.null(B2)) B2 <- bmsy
  n <- max(length(r), length(p), length(K), length(B1), length(B2))
  r <- rep(r, length.out = n)
  p <- rep(p, length.out = n)
  K <- rep(K, length.out = n)
  B1 <- rep(B1, length.out = n)
  B2 <- rep(B2, length.out = n)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (!all(is.finite(c(r[i], p[i], K[i], B1[i], B2[i]))) ||
        r[i] <= 0 || K[i] <= 0 || B1[i] <= 0 || B2[i] <= 0) next
    if (B1[i] >= B2[i]) {
      out[i] <- 0
      next
    }
    if (B2[i] >= K[i] * 0.999999) {
      out[i] <- Inf
      next
    }
    if (abs(p[i]) < 1e-8) {
      out[i] <- (1 / r[i]) * log(log(K[i] / B1[i]) / log(K[i] / B2[i]))
    } else {
      z1 <- (B1[i] / K[i])^p[i]
      z2 <- (B2[i] / K[i])^p[i]
      if (!is.finite(z1) || !is.finite(z2) || z1 == 1 || z2 == 1 ||
          z1 <= 0 || z2 <= 0) next
      ratio <- (z2 / (1 - z2)) * ((1 - z1) / z1)
      if (!is.finite(ratio) || ratio <= 0) next
      out[i] <- (1 / r[i]) * log(ratio)
    }
  }
  out
}

#' Euler check of closed-form Pella-Tomlinson Trebuild
#'
#' @param r,p,K,B1,B2 as in \code{\link{trebuild_pt}}.
#' @param dt step size in years.
#' @param max_t maximum years.
#' @return Numeric years, or \code{Inf} if the target is not reached.
#' @export
trebuild_pt_euler <- function(r, p, K, B1, B2, dt = 0.01, max_t = 400) {
  if (abs(p) < 1e-8) {
    spfun <- function(B) r * B * log(K / B)
  } else {
    spfun <- function(B) (r / p) * B * (1 - (B / K)^p)
  }
  B <- B1
  t <- 0
  while (B < B2 && t < max_t) {
    B <- B + dt * spfun(B)
    t <- t + dt
    if (!is.finite(B) || B <= 0) return(NA_real_)
  }
  if (B < B2) Inf else t
}

#' Trebuild in years from a B/BMSY start (K scaled to 1)
#'
#' @param r,p,shape Pella-Tomlinson parameters.
#' @param bbmsy0 starting B/BMSY (default 0.30).
#' @return Numeric years.
#' @export
trebuild_from_bbmsy <- function(r, p, shape, bbmsy0 = 0.30) {
  r <- as.numeric(r)
  p <- as.numeric(p)
  shape <- as.numeric(shape)
  K <- 1
  bmsy <- shape * K
  trebuild_pt(r = r, p = p, K = K, B1 = bbmsy0 * bmsy, B2 = bmsy,
              shape = shape, bmsy = bmsy)
}

#' Impose Fox or Schaefer shape, keeping MSY and BMSY
#'
#' Resets \code{K = BMSY/shape} and \code{r} from \code{FMSY} and shape.
#'
#' @param fmsy numeric MSY/BMSY.
#' @param bmsy numeric BMSY.
#' @param shape numeric BMSY/K to impose (0.4 Fox, 0.5 Schaefer).
#' @return \code{data.frame} with \code{r}, \code{p}, \code{k}, \code{virgin},
#'   \code{shape}, \code{fmsy}, \code{bmsy}.
#' @export
impose_shape <- function(fmsy, bmsy, shape) {
  n <- length(fmsy)
  shape <- rep(shape, length.out = n)
  p <- shape2p(shape)
  r <- r_from_fmsy_shape(fmsy, shape, p)
  K <- bmsy / shape
  data.frame(r = r, p = p, k = K, virgin = K, shape = shape,
             fmsy = fmsy, bmsy = bmsy)
}

#' SPiCT production exponent n from BMSY/K
#'
#' SPiCT uses \eqn{n = p + 1}. Schaefer is \code{n = 2}. Fox / other shapes
#' map through \code{\link{shape2p}}.
#'
#' @param shape numeric BMSY/K.
#' @return Numeric vector of SPiCT \code{n}.
#' @export
spict_n_from_shape <- function(shape) {
  n <- shape2p(shape) + 1
  n[abs(shape - 0.5) < 1e-8] <- 2
  n
}

#' Pin JABBA shape prior (effectively fix BMSY/K)
#'
#' Sets the log-scale SD (or CV slot) of \code{shape} in a JABBA
#' \code{pr.sd} \code{FLPar} to a small value.
#'
#' @param pr.sd \code{FLPar} of prior standard deviations.
#' @param cv numeric SD for the shape parameter (default 0.01).
#' @return The modified \code{pr.sd}.
#' @export
jabba_fix_shape <- function(pr.sd, cv = 0.01) {
  if (!("shape" %in% dimnames(pr.sd)$params))
    stop("pr.sd must contain a 'shape' parameter")
  pr.sd["shape"] <- cv
  pr.sd
}

#' Catch and exploitable-biomass index for surplus-production fits
#'
#' Builds the JABBA/SPiCT observation tables from an \code{FLStock}
#' (or a \code{tseries} data.frame): catch plus selectivity-weighted
#' exploitable biomass as the abundance index.
#'
#' @param object \code{FLStock} or a data.frame from \code{\link{tseries}}.
#' @return List with \code{catch} and \code{index} data.frames
#'   (\code{year}, \code{catch} / \code{index}).
#' @export
#' @seealso \code{\link{ebiomass}}, \code{\link{jabbaData}}
catch_eb_index <- function(object) {
  if (is.data.frame(object)) {
    ts <- object
  } else if (exists("tseries", mode = "function")) {
    ts <- tseries(object)
  } else {
    yrs <- as.integer(dimnames(FLCore::catch(object))$year)
    return(list(
      catch = data.frame(year = yrs, catch = as.numeric(FLCore::catch(object))),
      index = data.frame(year = yrs, index = as.numeric(ebiomass(object)))
    ))
  }
  catch <- data.frame(year = ts$year, catch = ts$catch)
  if ("eb" %in% names(ts)) {
    index <- data.frame(year = ts$year, index = ts$eb)
  } else if ("index" %in% names(ts)) {
    index <- data.frame(year = ts$year, index = ts$index)
  } else if ("ssb" %in% names(ts)) {
    index <- data.frame(year = ts$year, index = ts$ssb)
  } else {
    stop("No eb/index/ssb column in tseries output")
  }
  list(catch = catch, index = index)
}

#' Interpolate rebuilding time at a given starting depletion
#'
#' @param initial numeric starting B/BMSY (or SSB/BMSY) grid.
#' @param t_recover numeric rebuilding times at those starts.
#' @param xout starting depletion to interpolate (default 0.30).
#' @return Numeric interpolated years, or \code{NA}.
#' @export
interp_t_at <- function(initial, t_recover, xout = 0.30) {
  ok <- is.finite(initial) & is.finite(t_recover)
  if (sum(ok) < 2) return(NA_real_)
  initial <- initial[ok]
  t_recover <- t_recover[ok]
  o <- order(initial)
  approx(initial[o], t_recover[o], xout = xout, rule = 1, ties = "ordered")$y
}
