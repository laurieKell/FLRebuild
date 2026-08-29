# ICES stock-attribute accessors. Previously these S4 methods lived in FLCandy;
# they are implemented here so FLRebuild does not need that package.

.parNames <- function(x) {
  if (methods::is(x, "FLPar"))
    dimnames(x)$params
  else
    names(x)
}

`.parNames<-` <- function(x, value) {
  if (methods::is(x, "FLPar")) {
    dimnames(x)$params <- value
    return(x)
  }
  names(x) <- value
  x
}

.attrAsFLPar <- function(x) {
  if (is.null(x)) return(NULL)
  if (methods::is(x, "FLPar")) return(x)
  nms <- names(x)
  if (is.numeric(x) || is.logical(x)) {
    nms <- if (is.null(nms) || !length(nms)) paste0("p", seq_along(x)) else nms
    return(FLPar(array(as.numeric(x), dim = c(length(x), 1L),
                       dimnames = list(params = nms, iter = "1"))))
  }
  FLPar(x)
}

.bindStockPars <- function(object, fun) {
  plyr::ldply(object, function(x) {
    val <- fun(x)
    if (is.null(val)) return(NULL)
    as.data.frame(t(val))
  })
}

.keepPars <- function(x, keep) {
  nms <- .parNames(x)
  x[nms[nms %in% keep], ]
}

#' @rdname benchmark
#' @export
setMethod("benchmark", signature(object = "FLStock"), function(object) {
  if (!("benchmark" %in% names(attributes(object)))) {
    warning("No benchmark attribute found for this FLStock object.")
    return(NULL)
  }

  bm <- attributes(object)$benchmark
  if (is.logical(bm))
    return(FLPar(fmsy = NA, flim = NA, fpa = NA, blim = NA, bpa = NA, btrigger = NA))
  bm <- .attrAsFLPar(bm)
  .parNames(bm) <- tolower(.parNames(bm))
  .keepPars(bm, c("fmsy", "flim", "fpa", "blim", "bpa", "btrigger"))
})

#' @rdname benchmark
#' @export
setMethod("benchmark", signature(object = "FLStocks"), function(object) {
  .bindStockPars(object, benchmark)
})

#' @rdname benchmark
#' @export
setMethod("benchmark", signature(object = "list"), function(object) {
  .bindStockPars(object, benchmark)
})

#' @rdname benchmark
#' @export
setMethod("benchmark", signature(object = "FLBRP"), function(object) {
  refs <- refpts(object)
  pars <- FLPar(
    fmsy     = refs["msy", "harvest"],
    flim     = refs["lim", "harvest"],
    fpa      = refs["pa", "harvest"],
    blim     = refs["lim", "ssb"],
    bpa      = refs["pa", "ssb"],
    btrigger = refs["trigger", "ssb"]
  )
  pars[!is.na(pars)]
})

#' @rdname fishlife
#' @export
setMethod("fishlife", signature(object = "FLStock"), function(object) {
  if (!("fishlife" %in% names(attributes(object)))) {
    warning("No fishlife attribute found for this FLStock object.")
    return(NULL)
  }

  fl <- attributes(object)$fishlife
  if (is.logical(fl))
    return(FLPar(Fmsy = NA, Flim = NA, Fpa = NA, Blim = NA, Bpa = NA, Btrigger = NA))
  .attrAsFLPar(fl)
})

#' @rdname fishlife
#' @export
setMethod("fishlife", signature(object = "FLStocks"), function(object) {
  .bindStockPars(object, fishlife)
})

#' @rdname fishlife
#' @export
setMethod("fishlife", signature(object = "list"), function(object) {
  .bindStockPars(object, fishlife)
})

.eqsimNames <- function(nms) {
  nms <- paste0(tolower(substr(nms, 1L, 1L)), substring(nms, 2L))
  gsub("MSY", "msy", nms, fixed = TRUE)
}

#' @rdname eqsim
#' @export
setMethod("eqsim", signature(object = "FLStock"), function(object) {
  if (!("eqsim" %in% names(attributes(object)))) {
    warning("No eqsim attribute found for this FLStock object.")
    return(NULL)
  }

  eq <- attributes(object)$eqsim
  if (is.logical(eq))
    return(FLPar(
      catchequi = NA, bmsy = NA, b0 = NA, fmsyMedianC = NA,
      fmsyMedianL = NA, f5percRiskBlim = NA, flimEqsim = NA, r0 = NA
    ))
  eq <- .attrAsFLPar(eq)
  .parNames(eq) <- .eqsimNames(.parNames(eq))
  .keepPars(eq, c(
    "catchequi", "bmsy", "b0", "fmsyMedianC", "fmsyMedianL",
    "f5percRiskBlim", "flimEqsim", "r0"
  ))
})

#' @rdname eqsim
#' @export
setMethod("eqsim", signature(object = "FLStocks"), function(object) {
  .bindStockPars(object, eqsim)
})

#' @rdname eqsim
#' @export
setMethod("eqsim", signature(object = "list"), function(object) {
  .bindStockPars(object, eqsim)
})
