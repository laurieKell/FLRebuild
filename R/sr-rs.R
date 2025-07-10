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
rickerSSBMaxRec=function(b) return(1 / b)


