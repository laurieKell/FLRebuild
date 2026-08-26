#' Recruitment residuals and regime mean / SD
#'
#' Extracts stock–recruit residuals and expands STARS regime mean and SD
#' (\code{\link{rod}}, \code{\link{rodMn}}, \code{\link{rodSD}}) to each
#' observation. The default scale is log (FLSR residuals). Optional
#' simulated multiplicative deviations follow the 02.0 construction:
#' \code{exp(log(rlnoise(nits, rsdl \%=\% 0)) \%*\% rodSD + rodMn)}.
#'
#' @param object An \code{FLSR}, \code{FLSRs}, \code{FLBRP} or \code{FLBRPs}.
#' @param scale Residual scale: \code{"log"} (default), \code{"mult"}
#'   (\eqn{R/\hat{R}}), or \code{"raw"} (\eqn{R-\hat{R}}).
#' @param nits Integer. If \code{> 0}, simulate that many iterations of
#'   regime-scaled recruitment deviations.
#' @param pe Logical. For \code{FLBRP} / \code{FLBRPs}, also return surplus-
#'   production process error and its \code{rod()} table.
#' @param stock Optional \code{FLStock} / \code{FLStocks} used if
#'   \code{processError(FLBRP)} fails and \code{pe = TRUE}.
#' @param ... Passed to \code{\link{rod}} (\code{n}, \code{sig}, \ldots).
#'
#' @return A named list:
#' \describe{
#'   \item{residuals}{\code{FLQuant} or \code{FLQuants} of residuals}
#'   \item{mn}{Regime mean by year (\code{FLQuant} / \code{FLQuants})}
#'   \item{sd}{Regime SD by year (\code{FLQuant} / \code{FLQuants})}
#'   \item{rod}{\code{data.frame} of regime polygons (\code{sid} on collections)}
#'   \item{sim}{Simulated deviations if \code{nits > 0}, else \code{NULL}}
#'   \item{pe, pe_rod}{If \code{pe = TRUE}}
#' }
#'
#' @export
#' @seealso \code{\link{peDevs}}, \code{\link{rod}}, \code{\link{rodMn}},
#'   \code{\link{rodSD}}, \code{\link{processError}}
#' @examples
#' \dontrun{
#' data(nsher, package = "FLCore")
#' recDevs(nsher)
#' recDevs(nsher, scale = "mult", nits = 20)
#' }
setGeneric("recDevs", function(object, ...) standardGeneric("recDevs"))

.srFromFLBRP <- function(object) {
  sr <- attributes(object)[["sr"]]
  if (is(sr, "FLSR"))
    return(sr)
  if (is.list(sr) && length(sr) && is(sr[[1L]], "FLSR"))
    return(sr[[1L]])
  rsdl <- attributes(object)[["rec.residuals"]]
  if (is(rsdl, "FLQuant"))
    return(rsdl)
  NULL
}

.logResidualsFLSR <- function(object) {
  r <- try(residuals(object), silent = TRUE)
  if (!inherits(r, "try-error") && is(r, "FLQuant") &&
      any(is.finite(c(r))))
    return(r)
  log(rec(object)) - log(fitted(object))
}

.scaleResiduals <- function(log_r, object, scale) {
  scale <- match.arg(scale, c("log", "mult", "raw"))
  if (identical(scale, "log"))
    return(log_r)
  if (identical(scale, "mult"))
    return(exp(log_r))
  rec(object) - fitted(object)
}

.recDevsFromQuant <- function(rsdl, nits = 0L, scale = "log", ...) {
  mn <- try(rodMn(rsdl, ...), silent = TRUE)
  sd <- try(rodSD(rsdl, ...), silent = TRUE)
  if (inherits(mn, "try-error")) mn <- rsdl %=% NA_real_
  if (inherits(sd, "try-error")) sd <- rsdl %=% NA_real_
  rd <- try(rod(rsdl, ...), silent = TRUE)
  if (inherits(rd, "try-error") || is.null(rd))
    rd <- data.frame(year = integer(), data = numeric(),
                     regime = integer(), stringsAsFactors = FALSE)
  sim <- .recDevsSim(rsdl, nits, scale = scale, ...)
  list(residuals = rsdl, mn = mn, sd = sd, rod = rd, sim = sim)
}

.recDevsSim <- function(rsdl, nits, scale = "log", ...) {
  nits <- as.integer(nits)
  if (!is.finite(nits) || nits <= 0L)
    return(NULL)
  mn <- try(rodMn(rsdl, ...), silent = TRUE)
  sd <- try(rodSD(rsdl, ...), silent = TRUE)
  if (inherits(mn, "try-error") || inherits(sd, "try-error"))
    return(NULL)
  x <- rlnoise(nits, rsdl %=% 0)
  y <- log(x) %*% sd
  if (identical(scale, "log"))
    exp(y %+% mn)
  else
    mn %+% y
}

.combineRecDevs <- function(res, nits = 0L, pe = FALSE) {
  nms <- names(res)
  pick <- function(slot) {
    out <- lapply(res, `[[`, slot)
    names(out) <- nms
    Filter(Negate(is.null), out)
  }
  rsdl <- pick("residuals")
  mn   <- pick("mn")
  sd   <- pick("sd")
  sim  <- pick("sim")
  rods <- pick("rod")
  peQ  <- pick("pe")
  pe_r <- pick("pe_rod")
  rod_df <- if (length(rods)) {
    plyr::ldply(rods, .id = "sid")
  } else {
    data.frame(sid = character(), year = integer(), data = numeric(),
               regime = integer(), stringsAsFactors = FALSE)
  }
  pe_rod_df <- if (length(pe_r)) {
    plyr::ldply(pe_r, .id = "sid")
  } else {
    data.frame(sid = character(), year = integer(), data = numeric(),
               regime = integer(), stringsAsFactors = FALSE)
  }
  rtn <- list(
    residuals = FLQuants(rsdl),
    mn        = FLQuants(mn),
    sd        = FLQuants(sd),
    rod       = rod_df,
    sim       = if (as.integer(nits) > 0L && length(sim)) FLQuants(sim) else NULL)
  if (isTRUE(pe)) {
    rtn$pe <- if (length(peQ)) FLQuants(peQ) else FLQuants()
    rtn$pe_rod <- pe_rod_df
  }
  rtn
}

#' @rdname recDevs
#' @export
setMethod("recDevs", signature(object = "FLSR"),
          function(object, scale = "log", nits = 0L, ...) {
            scale <- match.arg(scale, c("log", "mult", "raw"))
            log_r <- .logResidualsFLSR(object)
            rsdl <- .scaleResiduals(log_r, object, scale)
            .recDevsFromQuant(rsdl, nits = nits, scale = scale, ...)
          })

#' @rdname recDevs
#' @export
setMethod("recDevs", signature(object = "FLSRs"),
          function(object, scale = "log", nits = 0L, ...) {
            nms <- names(object)
            if (is.null(nms) || !nzchar(nms[1L]))
              nms <- as.character(seq_along(object))
            res <- lapply(seq_along(object), function(i) {
              tryCatch(recDevs(object[[i]], scale = scale, nits = nits, ...),
                       error = function(e) NULL)
            })
            names(res) <- nms
            res <- Filter(Negate(is.null), res)
            .combineRecDevs(res, nits = nits, pe = FALSE)
          })

#' @rdname recDevs
#' @export
setMethod("recDevs", signature(object = "FLBRP"),
          function(object, scale = "log", nits = 0L, pe = FALSE,
                   stock = NULL, ...) {
            scale <- match.arg(scale, c("log", "mult", "raw"))
            sr <- .srFromFLBRP(object)
            if (is.null(sr))
              stop("FLBRP has no attached FLSR (attributes(object)$sr) ",
                   "or rec.residuals", call. = FALSE)
            if (is(sr, "FLQuant")) {
              rsdl <- sr
              if (identical(scale, "mult"))
                rsdl <- exp(rsdl)
              else if (identical(scale, "raw"))
                stop("scale = 'raw' needs an FLSR (rec and fitted)", call. = FALSE)
              rtn <- .recDevsFromQuant(rsdl, nits = nits, scale = scale, ...)
            } else {
              rtn <- recDevs(sr, scale = scale, nits = nits, ...)
            }
            if (isTRUE(pe)) {
              pd <- try(peDevs(object, scale = "raw", nits = 0L,
                               stock = stock, ...), silent = TRUE)
              if (!inherits(pd, "try-error") && is.list(pd)) {
                rtn$pe <- pd$residuals
                rtn$pe_rod <- pd$rod
              }
            }
            rtn
          })

#' @rdname recDevs
#' @export
setMethod("recDevs", signature(object = "FLBRPs"),
          function(object, scale = "log", nits = 0L, pe = FALSE,
                   stock = NULL, ...) {
            nms <- names(object)
            if (is.null(nms) || !nzchar(nms[1L]))
              nms <- as.character(seq_along(object))
            res <- lapply(nms, function(nm) {
              stk <- NULL
              if (!is.null(stock)) {
                if (is(stock, "FLStock"))
                  stk <- stock
                else if (nm %in% names(stock))
                  stk <- stock[[nm]]
              }
              tryCatch(
                recDevs(object[[nm]], scale = scale, nits = nits, pe = pe,
                        stock = stk, ...),
                error = function(e) NULL)
            })
            names(res) <- nms
            res <- Filter(Negate(is.null), res)
            .combineRecDevs(res, nits = nits, pe = pe)
          })

.peQuants <- function(object, stock = NULL) {
  peq <- try(processError(object), silent = TRUE)
  if (!inherits(peq, "try-error") && is(peq, "FLQuants"))
    return(peq)
  if (!is.null(stock) && is(stock, "FLStock")) {
    p <- try(pe(stock, object), silent = TRUE)
    if (!inherits(p, "try-error") && is(p, "FLQuant"))
      return(FLQuants(pe = p))
  }
  NULL
}

.peLead <- function(x) {
  y <- x %=% NA_real_
  ny <- dim(x)[2L]
  if (ny > 1L)
    y[, seq_len(ny - 1L)] <- x[, seq.int(2L, ny)]
  y
}

.peSeries <- function(object, stock = NULL, scale = "log") {
  scale <- match.arg(scale, c("log", "mult", "raw"))
  peq <- .peQuants(object, stock)
  if (is.null(peq))
    stop("Could not compute process error for FLBRP", call. = FALSE)
  have <- all(c("ssb", "catch", "sp") %in% names(peq))
  if (have) {
    ssb   <- peq[["ssb"]]
    Bnext <- .peLead(ssb)
    Bpred <- ssb %-% peq[["catch"]] %+% peq[["sp"]]
    return(switch(scale,
                  log  = log(Bnext) %-% log(Bpred),
                  mult = Bnext %/% Bpred,
                  raw  = if ("pe" %in% names(peq))
                    peq[["pe"]]
                  else
                    (Bnext %-% Bpred) %/% ssb))
  }
  p <- if ("pe" %in% names(peq)) peq[["pe"]] else peq[[1L]]
  switch(scale,
         raw  = p,
         log  = log(p %+% 1),
         mult = p %+% 1)
}

#' Process-error residuals and regime mean / SD
#'
#' Companion to \code{\link{recDevs}} for surplus-production process error on
#' \code{FLBRP} / \code{FLBRPs}. Default scale is log
#' (\eqn{\log B_{t+1}-\log(B_t-C_t+SP_t)}). \code{"mult"} is
#' \eqn{B_{t+1}/(B_t-C_t+SP_t)}; \code{"raw"} is the usual relative PE
#' \eqn{(B_{t+1}-B_t+C_t-SP_t)/B_t} from \code{\link{processError}}.
#'
#' @param object An \code{FLBRP} or \code{FLBRPs}.
#' @param scale \code{"log"} (default), \code{"mult"}, or \code{"raw"}.
#' @param nits Integer. If \code{> 0}, simulate regime-scaled PE deviations.
#' @param stock Optional \code{FLStock} / \code{FLStocks} if
#'   \code{processError} fails.
#' @param ... Passed to \code{\link{rod}}.
#'
#' @return Same list shape as \code{\link{recDevs}}: \code{residuals},
#'   \code{mn}, \code{sd}, \code{rod}, and \code{sim} if \code{nits > 0}.
#'
#' @export
#' @seealso \code{\link{recDevs}}, \code{\link{processError}}, \code{\link{pe}}
#' @examples
#' \dontrun{
#' peDevs(eq)
#' peDevs(eqls[["bh1"]], scale = "raw", stock = oms)
#' }
setGeneric("peDevs", function(object, ...) standardGeneric("peDevs"))

#' @rdname peDevs
#' @export
setMethod("peDevs", signature(object = "FLBRP"),
          function(object, scale = "log", nits = 0L, stock = NULL, ...) {
            scale <- match.arg(scale, c("log", "mult", "raw"))
            rsdl <- .peSeries(object, stock = stock, scale = scale)
            .recDevsFromQuant(rsdl, nits = nits, scale = scale, ...)
          })

#' @rdname peDevs
#' @export
setMethod("peDevs", signature(object = "FLBRPs"),
          function(object, scale = "log", nits = 0L, stock = NULL, ...) {
            nms <- names(object)
            if (is.null(nms) || !nzchar(nms[1L]))
              nms <- as.character(seq_along(object))
            res <- lapply(nms, function(nm) {
              stk <- NULL
              if (!is.null(stock)) {
                if (is(stock, "FLStock"))
                  stk <- stock
                else if (nm %in% names(stock))
                  stk <- stock[[nm]]
              }
              tryCatch(
                peDevs(object[[nm]], scale = scale, nits = nits,
                       stock = stk, ...),
                error = function(e) NULL)
            })
            names(res) <- nms
            res <- Filter(Negate(is.null), res)
            .combineRecDevs(res, nits = nits, pe = FALSE)
          })
