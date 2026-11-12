#' Rebuilding diagnostics from equilibrium reference points
#'
#' Projects an age-structured stock from a depleted spawning biomass and
#' records how long it takes to reach SSB at MSY. Three harvest policies are
#' compared: no fishing, constant fishing at FMSY, and a hockey-stick rule
#' that sets F to \code{FMSY * min(1, SSB / BMSY)}. Generation time is the
#' unfished mean age of spawners. \code{Bstar} is the lowest SSB from which
#' BMSY is reached within one generation time.
#'
#' The rebuild clock starts at \code{pRec * B0}. \code{SSB_pR} is reported
#' separately: the SSB at which expected recruitment is \code{pRec} of virgin
#' recruitment.
#'
#' @param object An \code{FLBRP}, or an \code{FLBRPs} collection.
#' @param pRec Depletion used as the starting SSB, as a fraction of virgin
#'   SSB. Also the recruitment fraction used for \code{SSB_pR}.
#' @param maxYears Projection horizon. Times beyond this are \code{NA}.
#' @param tol Relative shortfall allowed when deciding that SSB has reached
#'   BMSY. The default \code{0.01} counts recovery at 99 percent of BMSY.
#'   Fishing exactly at FMSY approaches BMSY asymptotically, so an exact
#'   crossing is not a finite time.
#' @param ... Ignored. Present so methods share a signature.
#'
#' @return For an \code{FLBRP}, an \code{FLPar} with one column per iteration
#'   and one row per diagnostic (\code{B0}, \code{Bmsy}, \code{Fmsy},
#'   \code{GT}, \code{SSB_pR}, \code{SSB_pR_frac}, \code{Trec_F0},
#'   \code{Trec_Fmsy}, \code{Trec_HCR}, the same times in generation units,
#'   and \code{Bstar} under \code{F = 0} and \code{F = FMSY}). For an
#'   \code{FLBRPs}, a named list of those objects, one per stock, so
#'   \code{plyr::ldply(pRebuild(brps), model.frame)} is one row per stock.
#'
#' @examples
#' \dontrun{
#' data(ple4, package = "FLCore")
#' eq <- brp(FLBRP(ple4))
#' pRebuild(eq)
#' pRebuild(FLBRPs(A = eq, B = eq))
#' }
#' @export
setGeneric("pRebuild", function(object, ...) standardGeneric("pRebuild"))

#' @rdname pRebuild
#' @export
setMethod("pRebuild", signature(object = "FLBRP"),
          function(object, pRec = 0.20, maxYears = 200, tol = 0.01, ...) {
            if (!is.numeric(pRec) || length(pRec) != 1L ||
                !is.finite(pRec) || pRec <= 0 || pRec >= 1)
              stop("'pRec' must be a number between 0 and 1", call. = FALSE)
            if (!is.numeric(maxYears) || length(maxYears) != 1L || maxYears < 1)
              stop("'maxYears' must be a positive number", call. = FALSE)
            if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
                tol < 0 || tol >= 1)
              stop("'tol' must be a number in [0, 1)", call. = FALSE)

            ni <- dims(object)$iter
            param_names <- c("B0", "Bmsy", "Fmsy", "GT",
                             "SSB_pR", "SSB_pR_frac",
                             "Trec_F0", "Trec_Fmsy", "Trec_HCR",
                             "Trec_F0_GT", "Trec_Fmsy_GT", "Trec_HCR_GT",
                             "Bstar_F0", "Bstar_Fmsy",
                             "Bstar_F0_frac", "Bstar_Fmsy_frac")
            mat <- matrix(NA_real_, nrow = length(param_names), ncol = ni,
                          dimnames = list(params = param_names, iter = seq_len(ni)))
            for (i in seq_len(ni)) {
              obj_i <- if (ni == 1L) object else iter(object, i)
              mat[, i] <- .pRebuild_one(obj_i, pRec = pRec, maxYears = maxYears,
                                        tol = tol)
            }
            FLPar(mat, units = "NA")
          })

#' @rdname pRebuild
#' @export
setMethod("pRebuild", signature(object = "FLBRPs"),
          function(object, pRec = 0.20, maxYears = 200, tol = 0.01, ...) {
            if (!length(object))
              return(list())
            nms <- names(object)
            if (is.null(nms))
              nms <- character(length(object))
            blank <- !nzchar(nms)
            if (any(blank))
              nms[blank] <- as.character(which(blank))
            res <- lapply(seq_along(object), function(i)
              pRebuild(object[[i]], pRec = pRec, maxYears = maxYears, tol = tol))
            names(res) <- nms
            res
          })

# One iteration. Returns a named numeric vector aligned with the FLPar rows.
.pRebuild_one <- function(brp, pRec, maxYears, tol) {
  prep <- .pRebuild_prep(brp)
  rp <- refpts(prep$eql)
  B0 <- c(rp["virgin", "ssb"])[1]
  Bmsy <- c(rp["msy", "ssb"])[1]
  Fmsy <- c(rp["msy", "harvest"])[1]
  GT <- .generation_time(prep$eql)

  SSB_pR <- .ssb_at_p_recruit(prep$eql, pRec, B0)
  init <- pRec * B0

  T_F0 <- .safe_time(prep, init, Bmsy, targetF = 0, maxYears = maxYears, tol = tol)
  T_Fmsy <- .safe_time(prep, init, Bmsy, targetF = Fmsy, maxYears = maxYears,
                       tol = tol)
  T_HCR <- .safe_time(prep, init, Bmsy, targetF = "hcr", maxYears = maxYears,
                      Fmsy = Fmsy, tol = tol)

  Bstar_F0 <- .bstar_one_gen(prep, GT, Bmsy, targetF = 0, Fmsy = Fmsy, tol = tol)
  Bstar_Fmsy <- .bstar_one_gen(prep, GT, Bmsy, targetF = Fmsy, Fmsy = Fmsy,
                               tol = tol)

  c(B0 = B0,
    Bmsy = Bmsy,
    Fmsy = Fmsy,
    GT = GT,
    SSB_pR = SSB_pR,
    SSB_pR_frac = SSB_pR / B0,
    Trec_F0 = T_F0,
    Trec_Fmsy = T_Fmsy,
    Trec_HCR = T_HCR,
    Trec_F0_GT = T_F0 / GT,
    Trec_Fmsy_GT = T_Fmsy / GT,
    Trec_HCR_GT = T_HCR / GT,
    Bstar_F0 = Bstar_F0,
    Bstar_Fmsy = Bstar_Fmsy,
    Bstar_F0_frac = Bstar_F0 / B0,
    Bstar_Fmsy_frac = Bstar_Fmsy / B0)
}

# Equilibrium biology used by the projection. Selectivity is F-at-age / fbar
# so a target F is an fbar, not an apical F. Year 1 of an FLBRP is often F = 0,
# where harvest is zero and cannot define selectivity.
.pRebuild_prep <- function(brp) {
  eql <- FLBRP::brp(brp)
  M <- c(m(eql)[, 1])
  matv <- c(mat(eql)[, 1])
  swt <- c(stock.wt(eql)[, 1])
  mspwn <- c(m.spwn(eql)[, 1])
  nage <- length(M)
  if (length(matv) == 1L) matv <- rep(matv, nage)
  if (length(swt) == 1L) swt <- rep(swt, nage)
  if (length(mspwn) == 1L) mspwn <- rep(mspwn, nage)

  fb <- c(fbar(eql))
  yr0 <- which.min(fb)
  Neq <- c(stock.n(eql)[, yr0])
  ssb_eq <- sum(Neq * matv * swt * exp(-M * mspwn))

  yrF <- which.max(fb)
  sel <- if (length(yrF) == 1L && is.finite(fb[yrF]) && fb[yrF] > 0)
    c(harvest(eql)[, yrF]) / fb[yrF]
  else
    rep(0, nage)
  sel[!is.finite(sel)] <- 0

  list(eql = eql, M = M, mat = matv, swt = swt, mspwn = mspwn,
       Neq = Neq, ssb_eq = ssb_eq, sel = sel, nage = nage,
       sr = .make_sr(eql))
}

.generation_time <- function(eql) {
  gt <- tryCatch(c(genTimeFn(eql)), error = function(e) NA_real_)
  gt <- gt[is.finite(gt) & gt > 0]
  if (length(gt))
    return(gt[1])
  yr <- which.min(c(fbar(eql)))
  n <- c(stock.n(eql)[, yr])
  matv <- c(mat(eql)[, 1])
  swt <- c(stock.wt(eql)[, 1])
  age <- as.numeric(dimnames(stock.n(eql))$age)
  w <- n * matv * swt
  if (!length(w) || !is.finite(sum(w)) || sum(w) <= 0)
    return(NA_real_)
  sum(age * w) / sum(w)
}

.par1 <- function(brp, nm) {
  p <- dimnames(params(brp))$params
  if (is.null(p) || !nm %in% p)
    return(NA_real_)
  v <- c(params(brp)[nm])[1]
  if (!is.finite(v)) NA_real_ else v
}

# SSB at which expected recruitment is p * R0.
# Beverton-Holt R = a S / (b + S):
#   S = p * b * B0 / (b + B0 * (1 - p))
.ssb_at_p_recruit <- function(brp, p, B0) {
  if (!is.finite(B0) || B0 <= 0)
    return(NA_real_)
  nm <- tryCatch(SRModelName(model(brp)), error = function(e) "")
  a <- .par1(brp, "a")
  b <- .par1(brp, "b")
  if (grepl("bevholt", nm, ignore.case = TRUE) && is.finite(b))
    return(p * b * B0 / (b + B0 * (1 - p)))
  if (grepl("segreg|hockey", nm, ignore.case = TRUE) && is.finite(b))
    return(p * b)
  sr <- .make_sr(brp)
  r0 <- tryCatch(sr(B0), error = function(e) NA_real_)
  if (!is.finite(r0))
    return(NA_real_)
  target <- p * r0
  tryCatch(uniroot(function(s) sr(s) - target, c(1e-8, B0))$root,
           error = function(e) NA_real_)
}

.make_sr <- function(brp) {
  nm <- tryCatch(SRModelName(model(brp)), error = function(e) "")
  a <- .par1(brp, "a")
  b <- .par1(brp, "b")
  if (grepl("bevholt", nm, ignore.case = TRUE) && is.finite(a) && is.finite(b))
    return(function(ssb) a * ssb / (b + ssb))
  if (grepl("ricker", nm, ignore.case = TRUE) && is.finite(a) && is.finite(b))
    return(function(ssb) a * ssb * exp(-b * ssb))
  if (grepl("segreg|hockey", nm, ignore.case = TRUE) && is.finite(a) && is.finite(b))
    return(function(ssb) ifelse(ssb >= b, a * b, a * ssb))
  if (grepl("mean", nm, ignore.case = TRUE) && is.finite(a))
    return(function(ssb) a)
  function(ssb) {
    v <- tryCatch(c(predict(as(brp, "FLSR"), ssb = FLQuant(ssb))),
                  error = function(e) NA_real_)
    if (!length(v)) NA_real_ else v[1]
  }
}

.safe_time <- function(prep, init, target, targetF, maxYears, Fmsy = NULL,
                     tol = 0.01) {
  tryCatch(
    .rebuild_time(prep, init = init, target = target, targetF = targetF,
                  maxYears = maxYears, Fmsy = Fmsy, tol = tol),
    error = function(e) NA_real_)
}

# First year in which SSB reaches `target`, starting from `init`.
# Returns 0 if the stock is already at the target, and NA if it has not
# arrived within maxYears.
.rebuild_time <- function(prep, init, target, targetF, maxYears, Fmsy = NULL,
                        tol = 0.01) {
  if (!is.finite(init) || !is.finite(target) || init <= 0 || target <= 0)
    return(NA_real_)
  goal <- target * (1 - tol)
  if (init >= goal)
    return(0)
  if (!is.finite(prep$ssb_eq) || prep$ssb_eq <= 0)
    return(NA_real_)

  maxYears <- as.integer(ceiling(maxYears))
  N <- prep$Neq * (init / prep$ssb_eq)
  ssb_now <- init
  sr <- prep$sr

  for (y in seq_len(maxYears)) {
    Fy <- if (identical(targetF, "hcr")) {
      if (!is.finite(Fmsy)) return(NA_real_)
      Fmsy * min(1, ssb_now / target)
    } else {
      as.numeric(targetF)[1]
    }
    if (!is.finite(Fy) || Fy < 0)
      return(NA_real_)
    Fa <- prep$sel * Fy
    Za <- Fa + prep$M
    R <- sr(ssb_now)
    if (!is.finite(R) || R < 0)
      return(NA_real_)
    Nn <- numeric(prep$nage)
    Nn[1] <- R
    if (prep$nage > 1L) {
      surv <- exp(-Za[seq_len(prep$nage - 1L)])
      Nn[2:prep$nage] <- N[seq_len(prep$nage - 1L)] * surv
      Nn[prep$nage] <- Nn[prep$nage] + N[prep$nage] * exp(-Za[prep$nage])
    }
    N <- Nn
    ssb_now <- sum(N * prep$mat * prep$swt * exp(-prep$M * prep$mspwn))
    if (is.finite(ssb_now) && ssb_now >= goal)
      return(y)
  }
  NA_real_
}

# Lowest SSB from which BMSY is reached within one generation time.
.bstar_one_gen <- function(prep, GT, Bmsy, targetF, Fmsy = NULL, tol = 0.01) {
  if (!is.finite(GT) || GT <= 0 || !is.finite(Bmsy) || Bmsy <= 0)
    return(NA_real_)
  horizon <- max(ceiling(GT), 1)
  t_at <- function(s) {
    t <- .safe_time(prep, s, Bmsy, targetF = targetF, maxYears = horizon,
                    Fmsy = Fmsy, tol = tol)
    if (!is.finite(t)) Inf else t
  }
  lo <- Bmsy * 1e-4
  hi <- Bmsy * (1 - 1e-8)
  if (t_at(lo) <= GT)
    return(lo)
  if (t_at(hi) > GT)
    return(NA_real_)
  for (i in seq_len(16L)) {
    mid <- exp((log(lo) + log(hi)) / 2)
    if (t_at(mid) <= GT) hi <- mid else lo <- mid
  }
  hi
}
