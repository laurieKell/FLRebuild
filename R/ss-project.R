# Stochastic SS3 projections with FLasher::fwd and recruitment residuals.
# Optional MCMC conditioning: scale terminal SSB to posterior draws and use
# each draw's SSB_MSY / annF_MSY for Kobe ratios.

.ssOneUnit <- function(stk, unit = "F") {
  out <- stk
  for (sl in slotNames(stk)) {
    if (!is(slot(stk, sl), "FLQuant")) next
    q <- slot(stk, sl)
    un <- dimnames(q)$unit
    if (is.null(un) || length(un) == 1L) {
      dimnames(q)$unit <- "unique"
      slot(out, sl) <- q
      next
    }
    i <- if (unit %in% un) unit else un[1L]
    q1 <- q[, , i]
    dimnames(q1)$unit <- "unique"
    slot(out, sl) <- q1
  }
  out
}

.ssCatchShare <- function(stk, unit = "F", year = dims(stk)$maxyear) {
  q <- catch(stk)
  un <- dimnames(q)$unit
  if (is.null(un) || length(un) == 1L) return(1)
  if (!unit %in% un) return(1)
  as.numeric(q[, ac(year), unit]) / sum(as.numeric(q[, ac(year)]))
}

.ssExtend <- function(stk, years) {
  end <- max(years)
  out <- FLCore::window(stk, end = end)
  last <- dims(stk)$maxyear
  slots <- c("m", "mat", "stock.wt", "catch.wt", "landings.wt", "discards.wt",
             "m.spwn", "harvest.spwn", "harvest", "landings.n", "discards.n",
             "catch.n")
  for (sl in slots) {
    if (!sl %in% slotNames(out)) next
    q <- slot(out, sl)
    for (y in years)
      if (ac(y) %in% dimnames(q)$year)
        q[, ac(y)] <- q[, ac(last)]
    slot(out, sl) <- q
  }
  out
}

.ssResiduals <- function(sr, years, nits, method = c("hist", "none", "mean"),
                         seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  dn <- list(age = "0", year = ac(years), iter = seq_len(nits))
  if (identical(method, "none") || identical(method, "mean"))
    return(FLCore::FLQuant(1, dimnames = dn))
  ## Historical multiplicative residuals in year order (estimated years only).
  ## Resample contiguous blocks (circular) so each iteration keeps the
  ## autocorrelation of the recdevs; iid year draws make jagged trajectories.
  rq <- FLCore::residuals(sr)
  yrs <- as.numeric(dimnames(rq)$year)
  rs <- c(rq)
  ok <- is.finite(rs) & rs > 0 & abs(rs - 1) > 1e-6 & is.finite(yrs)
  rsv <- rs[ok][order(yrs[ok])]
  if (!length(rsv)) {
    rsv <- rs[is.finite(rs) & rs > 0]
    if (!length(rsv)) return(FLCore::FLQuant(1, dimnames = dn))
  }
  rsv <- exp(log(rsv) - mean(log(rsv)))
  nh <- length(rsv)
  ny <- length(years)
  out <- FLCore::FLQuant(1, dimnames = dn)
  for (i in seq_len(nits)) {
    start <- sample.int(nh, 1L)
    idx <- ((start - 1L + seq_len(ny) - 1L) %% nh) + 1L
    out[, ac(years), , , , i] <- rsv[idx]
  }
  out
}

## Posterior draws for terminal SSB and MSY refpts from derived_posteriors.sso
.ssMcmcStart <- function(mcmc, year, nits, seed = NULL) {
  dat <- if (exists("mcmcTable", mode = "function", inherits = TRUE))
    mcmcTable(mcmc)
  else
    getFromNamespace("mcmcTable", "FLRebuild")(mcmc)
  ssbCol <- paste0("SSB_", year)
  if (!ssbCol %in% names(dat))
    stop("no ", ssbCol, " in ", mcmc)
  need <- c("SSB_MSY", "annF_MSY")
  miss <- setdiff(need, names(dat))
  if (length(miss))
    stop("missing columns in ", mcmc, ": ", paste(miss, collapse = ", "))
  ok <- is.finite(dat[[ssbCol]]) & dat[[ssbCol]] > 0 &
    is.finite(dat$SSB_MSY) & dat$SSB_MSY > 0 &
    is.finite(dat$annF_MSY) & dat$annF_MSY > 0
  dat <- dat[ok, , drop = FALSE]
  if (!nrow(dat))
    stop("no finite MCMC rows for ", ssbCol, " / SSB_MSY / annF_MSY")
  if (!is.null(seed)) set.seed(seed)
  nits <- as.integer(nits)
  idx <- if (nrow(dat) >= nits)
    sample.int(nrow(dat), nits, replace = FALSE)
  else
    sample.int(nrow(dat), nits, replace = TRUE)
  list(
    iter = dat$Iter[idx],
    ssb = as.numeric(dat[[ssbCol]][idx]),
    bmsy = as.numeric(dat$SSB_MSY[idx]),
    fmsy = as.numeric(dat$annF_MSY[idx]))
}

## Scale numbers in the assessment terminal year so SSB matches each MCMC draw
.ssConditionMcmc <- function(stk, ssbTarget, year) {
  y <- ac(year)
  stk1 <- FLCore::iter(stk, 1)
  ssb0 <- FLCore::ssb(stk1)[, y]
  if (!is.finite(c(ssb0)) || c(ssb0) <= 0)
    stop("SSB in ", year, " of FLR stock is not positive")
  n1 <- FLCore::stock.n(stk1)[, y]
  scale <- FLCore::FLQuant(ssbTarget,
                           dimnames = list(iter = dimnames(FLCore::stock.n(stk))$iter)) %/%
    ssb0
  nY <- FLCore::propagate(n1, length(ssbTarget))
  dimnames(nY)$iter <- dimnames(scale)$iter
  for (i in seq_along(ssbTarget))
    FLCore::iter(nY, i) <- n1 * c(FLCore::iter(scale, i))
  FLCore::stock.n(stk)[, y] <- nY
  stk
}

#' @rdname ssProject
#' @param catch Numeric vector of constant annual catches (tonnes, total).
#' @param nits Number of iterations (MCMC draws and/or process-error paths).
#' @param end Last projection year.
#' @param residuals \code{"hist"} circular-block resamples the historical
#'   multiplicative recruitment residual series (preserves recdev
#'   autocorrelation); \code{"none"} / \code{"mean"} use 1.
#' @param unit Sex unit to project for two-sex stocks (default \code{"F"}).
#' @param seed Optional RNG seed for residual / MCMC draws.
#' @param bmsy,fmsy Optional reference points; default from SS3 MLE
#'   \code{SSB_MSY} and \code{annF_MSY}, or per MCMC draw when \code{mcmc} is set.
#' @param mcmc Optional path to \code{derived_posteriors.sso}. When supplied,
#'   each iteration is conditioned on a posterior draw of terminal SSB and
#'   uses that draw's \code{SSB_MSY} and \code{annF_MSY}.
#' @param ssbYear Assessment year whose \code{SSB_<year>} is used to condition
#'   the FLR stock (default: stock \code{maxyear}).
#' @export
setMethod("ssProject", signature(object = "character"),
  function(object, catch = 0, nits = 100L, end = NULL,
           residuals = c("hist", "none", "mean"), unit = "F",
           seed = NULL, bmsy = NULL, fmsy = NULL, mcmc = NULL,
           ssbYear = NULL, ...) {
  if (!requireNamespace("ss3om", quietly = TRUE))
    stop("ss3om is required for ssProject")
  if (!requireNamespace("FLasher", quietly = TRUE))
    stop("FLasher is required for ssProject")
  residuals <- match.arg(residuals)
  dir <- normalizePath(object, winslash = "/", mustWork = TRUE)
  ss  <- r4ss::SS_output(dir, verbose = FALSE, printstats = FALSE, covar = FALSE)
  stk <- ss3om::readFLSss3(dir)
  sr  <- ss3om::readFLSRss3(dir)
  stk1 <- .ssOneUnit(stk, unit = unit)
  share <- .ssCatchShare(stk, unit = unit, year = dims(stk)$maxyear)
  ## One-sex stock: bevholtss3 R0 is total recruitment; keep only this unit's
  ## share (sratio), otherwise female SSB rebuilds twice as fast as SS3.
  if ("sratio" %in% dimnames(FLCore::params(sr))$params) {
    sratio <- c(FLCore::params(sr)["sratio"])
    if (is.finite(sratio) && sratio > 0 && sratio < 1)
      FLCore::params(sr)["R0"] <- FLCore::params(sr)["R0"] * sratio
  }
  if (is.null(end))
    end <- ss$endyr + 46L
  endHist <- dims(stk1)$maxyear
  years <- seq(endHist + 1L, end)
  stkW <- .ssExtend(stk1, years)
  nits <- as.integer(nits)

  ## No process error and no MCMC → single MLE trajectory (fast check path)
  if (residuals %in% c("none", "mean") && is.null(mcmc) && nits > 1L)
    nits <- 1L

  mc <- NULL
  if (!is.null(mcmc)) {
    mcmc <- normalizePath(mcmc, winslash = "/", mustWork = TRUE)
    if (is.null(ssbYear)) ssbYear <- endHist
    mc <- .ssMcmcStart(mcmc, year = ssbYear, nits = nits, seed = seed)
    bmsy <- mc$bmsy
    fmsy <- mc$fmsy
  } else {
    if (is.null(bmsy))
      bmsy <- ss$derived_quants["SSB_MSY", "Value"]
    if (is.null(fmsy))
      fmsy <- ss$derived_quants["annF_MSY", "Value"]
    bmsy <- rep(as.numeric(bmsy), length.out = nits)
    fmsy <- rep(as.numeric(fmsy), length.out = nits)
  }

  stkP <- FLCore::propagate(stkW, nits)
  if (!is.null(mc))
    stkP <- .ssConditionMcmc(stkP, mc$ssb, year = ssbYear)
  ## Residual RNG: if MCMC already used seed, advance from there
  dev <- .ssResiduals(sr, years, nits, method = residuals,
                      seed = if (is.null(mc)) seed else NULL)

  rows <- vector("list", length(catch))
  for (i in seq_along(catch)) {
    tac <- catch[i]
    message("ssProject: catch = ", tac, " t (", residuals,
            if (!is.null(mc)) ", MCMC" else ", MLE",
            ", nits = ", nits, ")")
    ct <- FLCore::propagate(
      FLCore::FLQuant(tac * share, dimnames = list(year = ac(years)), units = "t"),
      nits)
    pr <- FLasher::fwd(stkP, catch = ct, sr = sr, residuals = dev)
    ssbQ <- FLCore::ssb(pr)[, ac(years)]
    fQ   <- FLCore::fbar(pr)[, ac(years)]
    ## Per-iter refpts (FLQuant iter dim); %/% recycles over years
    bmsyQ <- FLCore::FLQuant(bmsy, dimnames = list(iter = dimnames(ssbQ)$iter))
    fmsyQ <- FLCore::FLQuant(fmsy, dimnames = list(iter = dimnames(fQ)$iter))
    ## drop=FALSE keeps length-1 iter so rows stay aligned for nits = 1
    flqDf <- function(q) {
      d <- as.data.frame(q, drop = FALSE)
      data.frame(
        year = as.integer(as.character(d$year)),
        iter = as.integer(as.character(d$iter)),
        data = d$data)
    }
    dfB <- flqDf(ssbQ %/% bmsyQ)
    dfF <- flqDf(fQ %/% fmsyQ)
    rows[[i]] <- data.frame(
      catch = tac,
      year = dfB$year,
      iter = dfB$iter,
      stock = dfB$data,
      harvest = dfF$data,
      bmsy = bmsy[dfB$iter],
      fmsy = fmsy[dfF$iter],
      mcmcIter = if (is.null(mc)) NA_integer_ else as.integer(mc$iter[dfB$iter]))
  }
  kb <- do.call(rbind, rows)
  rownames(kb) <- NULL
  attr(kb, "bmsy") <- if (length(unique(round(bmsy, 6))) == 1L) bmsy[1] else bmsy
  attr(kb, "fmsy") <- if (length(unique(round(fmsy, 8))) == 1L) fmsy[1] else fmsy
  attr(kb, "catchShare") <- share
  attr(kb, "residuals") <- residuals
  attr(kb, "unit") <- unit
  attr(kb, "dir") <- dir
  attr(kb, "mcmc") <- mcmc
  attr(kb, "ssbYear") <- if (is.null(mcmc)) NA_integer_ else as.integer(ssbYear)
  kb
})
