# =============================================================================
# Stock-Recruitment Relationship Functions
# =============================================================================

require(gsl)

# =============================================================================
# S4 Generics
# =============================================================================

#' Generic for inverse stock-recruitment relationships
#' 
#' @param params FLPar object with SRR parameters
#' @param rec Recruitment values (numeric, FLPar, or FLQuant)
#' @param ... Additional arguments
#' @return SSB values corresponding to the given recruitment
#' @export
setGeneric("invSRR", function(params, rec, ...) standardGeneric("invSRR"))

#' Create reference point FLPar object
#' 
#' @param object An object (typically FLBRP)
#' @param ref Reference point name (character)
#' @param value Value for recruitment (numeric, default NA)
#' @param ... Additional arguments
#' @return FLPar object with reference point structure
#' @export
setGeneric("refCreate", function(object, ref, value = NA, ...) standardGeneric("refCreate"))

#' Maximum recruitment for a stock-recruitment model
#' 
#' @param object An object (typically FLBRP)
#' @param ... Additional arguments
#' @return Maximum recruitment (numeric or FLPar)
#' @export
setGeneric("rmax", function(object, ratio) standardGeneric("rmax"))

#' Recruitment at MSY for a stock-recruitment model
#' 
#' @param object An object (typically FLBRP)
#' @param ... Additional arguments
#' @return Recruitment at MSY (numeric or FLPar)
#' @export
setGeneric("rmsy", function(object, ...) standardGeneric("rmsy"))

# =============================================================================
# S4 Methods
# =============================================================================

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

#' @rdname refCreate
#' @export
setMethod("refCreate", signature(object = "FLBRP"), function(object, ref, value = NA, ...) {
  rtn = list(refpt=ref,value=ref,
             quant = c("harvest", "yield", "rec",
                       "ssb", "biomass", "revenue",
                       "cost", "profit"),
             iter = ifelse(is.na(value), 1, seq(value)))
  rtn = FLPar(NA, dimnames = rtn)
  rtn[ref, value] = value
  rtn})

#' @rdname rmax
#' @export
setMethod("rmax", signature(object="FLBRP",ratio="missing"), 
          function(object) {
  switch(SRModelName(model(object)),
         bevholt = params(object)["a"],
         ricker  = 1/params(object)["b"],
         segreg  = params(object)["a"] %/% params(object)["a"])})

setMethod("rmax", signature(object="FLBRP",ratio="numeric"), 
     function(object,ratio) {
        #invSRR
        #refpts(object)=refCreate(object, "rmax", "rec", rmax(object)*ratio)
        refpts(object)=refCreate(object, "rmax","ssb", srrInv(object,rmax(object)*ratio))
        computeRefpts(object)})

#' @rdname rmsy
#' @export
setMethod("rmsy", signature(object="FLBRP"), 
      function(object,ratio=1.0) {
        #refpts(object)=refCreate(object, "rmsy","rec", refpts(object)["msy","rec"]*ratio)
        refpts(object)=refCreate(object, "rmsy","ssb", srrInv(object,refpts(object)["msy","rec"]*ratio))
        computeRefpts(object)})

# =============================================================================
# Helper Functions (not S4)
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
  SSB = -1 %/% (params["b"] %*% arg)
  return(SSB)
}

#' Maximum recruitment SSB for Ricker model
#' 
#' @param b Ricker parameter b
#' @return SSB at maximum recruitment
#' @keywords internal
rickerMaxRec = function(b) return(1/b)

