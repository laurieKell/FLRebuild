# =============================================================================
# Stock-Recruitment Relationship Functions
# =============================================================================

# SSB for a given recruitment level (returns both lower and upper solutions)
require(gsl)

#' Generic for inverse stock-recruitment relationships
#' 
#' @param params FLPar object with SRR parameters
#' @param rec Recruitment values (numeric, FLPar, or FLQuant)
#' @param ... Additional arguments
#' @return SSB values corresponding to the given recruitment
#' @export
setGeneric("invSRR", function(params, rec, ...) standardGeneric("invSRR"))

#' @rdname invSRR
#' @export
setMethod("invSRR", signature(params = "FLBRP", rec = "FLQuant"),
  function(params, rec) {
    switch(SRModelName(model(params)),
           bevholt = bevholtInv(params(params), rec),
           ricker = rickerInv(params(params), rec))
  })

#' @rdname invSRR
#' @export
setMethod("invSRR", signature(params = "FLSR", rec = "FLQuant"),
  function(params, rec) {
    switch(SRModelName(model(params)),
           bevholt = bevholtInv(params(params), rec),
           ricker = rickerInv(params(params), rec))
  })

# =============================================================================
# Inverse SRR Helper Functions
# =============================================================================

#' Inverse Beverton-Holt stock-recruitment relationship
#' 
#' @param params FLPar with parameters 'a' and 'b'
#' @param rec Recruitment values
#' @return SSB values
#' @keywords internal
bevholtInv = function(params, rec) {
  return((rec %*% params["b"]) %/% (params["a"] %-% rec))
}

#' Inverse Ricker stock-recruitment relationship
#' 
#' @param params FLPar with parameters 'a' and 'b'
#' @param rec Recruitment values
#' @return SSB values
#' @keywords internal
rickerInv = function(params, rec) {
  arg = -(rec %*% params["b"]) %/% params["a"]
  arg[arg > 0] = NA
  arg[arg <= 0] = sapply(arg[arg <= 0], gsl::lambert_Wm1)
  
  SSB = -1 / (params["b"] %*% arg)
  
  return(SSB)
}

#' Maximum recruitment SSB for Ricker model
#' 
#' @param b Ricker parameter b
#' @return SSB at maximum recruitment
#' @keywords internal
rickerSSBMaxRec = function(b) {
  return(1 / b)
}

# =============================================================================
# Inverse SRR Generic Functions
# =============================================================================

#' Generic for SSB at a given recruitment level for Ricker model
#' @export
setGeneric("rickerRec", function(params, rec, ...) standardGeneric("rickerRec"))

#' Generic for SSB at a given recruitment level for Beverton-Holt model
#' @export
setGeneric("bevertonRec", function(params, rec, ...) standardGeneric("bevertonRec"))

#' Generic for SSB at a given recruitment level for Segmented Regression model
#' @export
setGeneric("segRegRec", function(params, rec, ...) standardGeneric("segRegRec"))

# =============================================================================
# Ricker Model Methods
# =============================================================================

#' Helper function for rickerRec with numeric inputs
#' @keywords internal
.rickerRec_numeric = function(params, rec) {
  rickerFun = function(S) as.numeric(params["a"]) * S * exp(-as.numeric(params["b"]) * S) - rec
  lower = tryCatch(uniroot(rickerFun, lower = .Machine$double.eps, upper = 1/as.numeric(params["b"]))$root, 
                   error = function(e) NA)
  upper = tryCatch(uniroot(rickerFun, lower = 1/as.numeric(params["b"]), upper = 10*(1/as.numeric(params["b"])))$root, 
                   error = function(e) NA)
  c(lower = lower, upper = upper)
}

#' @rdname rickerRec
#' @export
setMethod("rickerRec", signature(params = "FLPar", rec = "numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .rickerRec_numeric(params, r))
  })

#' @rdname rickerRec
#' @export
setMethod("rickerRec", signature(params = "FLPar", rec = "FLPar"),
  function(params, rec, ...) {
    n_iter = dims(rec)$iter
    res = sapply(seq_len(n_iter), function(i) .rickerRec_numeric(params, as.numeric(rec[i])))
    dimnames(res) = list(sol = c("lower","upper"), iter = seq_len(n_iter))
    FLPar(res)
  })

#' @rdname rickerRec
#' @export
setMethod("rickerRec", signature(params = "FLPar", rec = "FLQuant"),
  function(params, rec, ...) {
    n_iter = dim(rec)[6]
    res = sapply(seq_len(n_iter), function(i) .rickerRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    dimnames(res) = list(sol = c("lower","upper"), iter = seq_len(n_iter))
    FLQuant(res)
  })

# =============================================================================
# Beverton-Holt Model Methods
# =============================================================================

#' Helper function for bevertonRec with numeric inputs
#' @keywords internal
.bevertonRec_numeric = function(params, rec) {
  denom = as.numeric(params["a"]) - as.numeric(params["b"]) * rec
  if (denom <= 0) NA else rec / denom
}

#' @rdname bevertonRec
#' @export
setMethod("bevertonRec", signature(params = "FLPar", rec = "numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .bevertonRec_numeric(params, r))
  })

#' @rdname bevertonRec
#' @export
setMethod("bevertonRec", signature(params = "FLPar", rec = "FLPar"),
  function(params, rec, ...) {
    n_iter = dims(rec)$iter
    res = sapply(seq_len(n_iter), function(i) .bevertonRec_numeric(params, as.numeric(rec[i])))
    names(res) = seq_len(n_iter)
    FLPar(res)
  })

#' @rdname bevertonRec
#' @export
setMethod("bevertonRec", signature(params = "FLPar", rec = "FLQuant"),
  function(params, rec, ...) {
    n_iter = dim(rec)[6]
    res = sapply(seq_len(n_iter), function(i) .bevertonRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    names(res) = seq_len(n_iter)
    FLQuant(res)
  })

# =============================================================================
# Segmented Regression Model Methods
# =============================================================================

#' Helper function for segRegRec with numeric inputs
#' @keywords internal
.segRegRec_numeric = function(params, rec) {
  m1 = as.numeric(params["m1"])
  m2 = as.numeric(params["m2"])
  b1 = as.numeric(params["b1"])
  c = as.numeric(params["c"])
  b2 = m1 * c + b1 - m2 * c
  if (rec <= m1 * c + b1) (rec - b1) / m1 else (rec - b2) / m2
}

#' @rdname segRegRec
#' @export
setMethod("segRegRec", signature(params = "FLPar", rec = "numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .segRegRec_numeric(params, r))
  })

#' @rdname segRegRec
#' @export
setMethod("segRegRec", signature(params = "FLPar", rec = "FLPar"),
  function(params, rec, ...) {
    n_iter = dims(rec)$iter
    res = sapply(seq_len(n_iter), function(i) .segRegRec_numeric(params, as.numeric(rec[i])))
    names(res) = seq_len(n_iter)
    FLPar(res)
  })

#' @rdname segRegRec
#' @export
setMethod("segRegRec", signature(params = "FLPar", rec = "FLQuant"),
  function(params, rec, ...) {
    n_iter = dim(rec)[6]
    res = sapply(seq_len(n_iter), function(i) .segRegRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    names(res) = seq_len(n_iter)
    FLQuant(res)
  })

# =============================================================================
# Forward SRR Functions
# =============================================================================

#' Generic for expected recruitment given SSB and an SRR object
#' 
#' @param object FLBRP or FLSR object
#' @param ssb Spawning stock biomass (numeric or FLPar)
#' @return Expected recruitment
#' @export
setGeneric("recHat", function(object, ssb) standardGeneric("recHat"))

#' @rdname recHat
#' @export
setMethod("recHat", signature(object = "FLBRP", ssb = "numeric"),
  function(object, ssb) {
    model_name = as.character(model(object))
    pars = params(object)
    
    if (model_name == "ricker") {
      return(as.numeric(pars["a"]) * ssb * exp(-as.numeric(pars["b"]) * ssb))
    } else if (model_name == "bevholt") {
      return((as.numeric(pars["a"]) * ssb) / (1 + as.numeric(pars["b"]) * ssb))
    } else if (model_name == "segreg") {
      m1 = as.numeric(pars["m1"])
      m2 = as.numeric(pars["m2"])
      b1 = as.numeric(pars["b1"])
      c = as.numeric(pars["c"])
      b2 = m1 * c + b1 - m2 * c
      return(ifelse(ssb <= c, m1 * ssb + b1, m2 * ssb + b2))
    } else {
      stop("Unknown or unsupported SRR model in FLBRP: ", model_name)
    }
  })

#' @rdname recHat
#' @export
setMethod("recHat", signature(object = "FLSR", ssb = "numeric"),
  function(object, ssb) {
    model_name = as.character(model(object))
    pars = params(object)
    
    if (model_name == "ricker") {
      return(as.numeric(pars["a"]) * ssb * exp(-as.numeric(pars["b"]) * ssb))
    } else if (model_name == "bevholt") {
      return((as.numeric(pars["a"]) * ssb) / (1 + as.numeric(pars["b"]) * ssb))
    } else if (model_name == "segreg") {
      m1 = as.numeric(pars["m1"])
      m2 = as.numeric(pars["m2"])
      b1 = as.numeric(pars["b1"])
      c = as.numeric(pars["c"])
      b2 = m1 * c + b1 - m2 * c
      return(ifelse(ssb <= c, m1 * ssb + b1, m2 * ssb + b2))
    } else {
      stop("Unknown or unsupported SRR model in FLSR: ", model_name)
    }
  })

#' @rdname recHat
#' @export
setMethod("recHat", signature(object = "FLBRP", ssb = "FLPar"),
  function(object, ssb) {
    model_name = as.character(model(object))
    pars = params(object)
    n_iter = dims(ssb)$iter
    res = rep(NA_real_, n_iter)
    
    if (model_name == "ricker") {
      for (i in seq_len(n_iter)) {
        res[i] = as.numeric(pars["a"]) * as.numeric(ssb[i]) * exp(-as.numeric(pars["b"]) * as.numeric(ssb[i]))
      }
    } else if (model_name == "bevholt") {
      for (i in seq_len(n_iter)) {
        res[i] = (as.numeric(pars["a"]) * as.numeric(ssb[i])) / (1 + as.numeric(pars["b"]) * as.numeric(ssb[i]))
      }
    } else if (model_name == "segreg") {
      m1 = as.numeric(pars["m1"])
      m2 = as.numeric(pars["m2"])
      b1 = as.numeric(pars["b1"])
      c = as.numeric(pars["c"])
      b2 = m1 * c + b1 - m2 * c
      for (i in seq_len(n_iter)) {
        res[i] = if (as.numeric(ssb[i]) <= c) m1 * as.numeric(ssb[i]) + b1 else m2 * as.numeric(ssb[i]) + b2
      }
    } else {
      stop("Unknown or unsupported SRR model in FLBRP: ", model_name)
    }
    
    return(FLPar(res, dimnames = list(iter = seq_len(n_iter))))
  })

#' @rdname recHat
#' @export
setMethod("recHat", signature(object = "FLSR", ssb = "FLPar"),
  function(object, ssb) {
    model_name = as.character(model(object))
    pars = params(object)
    n_iter = dims(ssb)$iter
    res = rep(NA_real_, n_iter)
    
    if (model_name == "ricker") {
      for (i in seq_len(n_iter)) {
        res[i] = as.numeric(pars["a"]) * as.numeric(ssb[i]) * exp(-as.numeric(pars["b"]) * as.numeric(ssb[i]))
      }
    } else if (model_name == "bevholt") {
      for (i in seq_len(n_iter)) {
        res[i] = (as.numeric(pars["a"]) * as.numeric(ssb[i])) / (1 + as.numeric(pars["b"]) * as.numeric(ssb[i]))
      }
    } else if (model_name == "segreg") {
      m1 = as.numeric(pars["m1"])
      m2 = as.numeric(pars["m2"])
      b1 = as.numeric(pars["b1"])
      c = as.numeric(pars["c"])
      b2 = m1 * c + b1 - m2 * c
      for (i in seq_len(n_iter)) {
        res[i] = if (as.numeric(ssb[i]) <= c) m1 * as.numeric(ssb[i]) + b1 else m2 * as.numeric(ssb[i]) + b2
      }
    } else {
      stop("Unknown or unsupported SRR model in FLSR: ", model_name)
    }
    
    return(FLPar(res, dimnames = list(iter = seq_len(n_iter))))
  })
