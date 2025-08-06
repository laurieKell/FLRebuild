setGeneric("lag", function(object, ...)
  standardGeneric("lag"))
# }}}


# lag {{{

#' @title Estimate population dynamics lag
#' 
#' @description 
#' This method estimates the lag in population dynamics by finding the age with maximum SSB contribution.
#' Excludes the plusgroup from the calculation.
#' 
#' @param object An object of class `FLBRP`
#' @param ... Additional arguments (not currently used)
#' 
#' @return An FLQuant object with 1 quant dimension, containing the age indices
#' where SSB contribution is maximized, indicating the population dynamics lag
#' 
#' @docType methods
#' @rdname lag
#' 
#' @seealso [FLBRP] 
#' 
#' @examples
#' data(ple4brp)
#'
#' # Estimate lag using SSB age distribution
#' lag_time=lag(ple4brp)

setMethod('lag', signature(object='FLBRP'),
          function(object, ...){
            
            ftrend=fbar(object)%/%refpts(object)["msy","harvest"]
            
            # Remove plusgroup if it exists
            if (!is.na(range(object)["plusgroup"])){
              range(object)["plusgroup"] = NA
              object = object[-dim(m(object))[1]]
            }
            
            # Use ageMax to find age with maximum SSB contribution
            result = ageMax(ssb.age(object))
            
            return(mcf(FLQuants(lag=result,fbar=ftrend)))}) # }}}
