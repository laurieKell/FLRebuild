#' Calculate recovery time to BMSY under Pella-Tomlinson dynamics
#'
#' @description Numerically integrates the Pella-Tomlinson production function to estimate the time required 
#' for a population to recover from an initial biomass to BMSY, with optional fishing mortality during recovery.
#'
#' @param object Initial biomass (must be less than BMSY) or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param F Fishing mortality rate during recovery (default=0)
#' @param ... Additional arguments
#' @return Recovery time in years
#' @export
#' @examples
#' rTime(0.3, r=0.2, p=1.5)
#' rTime(0.3, r=0.2, p=1.5, F=0.05)
#' 
#' # With FLPar
#' library(FLCore)
#' params=FLPar(initial=0.3, r=0.2, p=1.5, bmsy=1)
#' rTime(params, F=0.05)
#' 
#' # With FLBRP object
#' data(ple4brp)
#' rTime(ple4brp, initial=0.3, F=0.05)
setGeneric("rTime", function(object, ...) standardGeneric("rTime"))

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="numeric"),
  function(object, r, p, bmsy=1, F=0, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    # Find the maximum length among all arguments
    n = max(length(object), length(r), length(p), length(bmsy), length(F))
    # Recycle all arguments to the same length
    initial = rep(object, length.out=n)
    r = rep(r, length.out=n)
    p = rep(p, length.out=n)
    bmsy = rep(bmsy, length.out=n)
    F = rep(F, length.out=n)

    # Check for NA values in any argument
    naMask = is.na(initial) | is.na(r) | is.na(p) | is.na(bmsy) | is.na(F)
    if (all(naMask)) {
      return(rep(NA_real_, n))
    }
 
    # Calculate carrying capacity (virgin biomass) from bmsy and p
    K = calculateCarryingCapacity(p, bmsy)
     
    # Effective growth rate under fishing mortality
    rEff = r - F
 
    # Check if recovery is possible under current fishing mortality
    invalid = rEff <= 0
    if (any(invalid, na.rm=TRUE)) {
      warning("Fishing mortality too high: recovery not possible for some elements")
    }

    # Initial and target biomass
    B0 = initial
    B1 = bmsy

    # Calculate recovery time
    t = calculateRecoveryTime(B0, B1, K, rEff, p, invalid, naMask)
    
    return(as.numeric(t))
  })

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="FLPar"),
  function(object, F=0, initial = NULL, ...) {
    # Get parameter names from FLPar
    param_names = dimnames(object)$params
    
    # Check if this is a pt() output (has r, p, bmsy but no initial)
    if (all(c("r", "p", "bmsy") %in% param_names) && !("initial" %in% param_names)) {
      if (is.null(initial)) {
        stop("FLPar from pt() requires 'initial' parameter to be provided")
      }
      
      nits = dims(object)$iter
      result = numeric(nits)
      
      for (i in seq(nits)) {
        result[i] = rTime(
          initial,
          r = as.numeric(object["r", i]),
          p = as.numeric(object["p", i]),
          bmsy = as.numeric(object["bmsy", i]),
          F = F
        )
      }
      
      return(result)
    }
    
    # Original behavior for FLPar with initial, r, p
    if (!all(c("initial", "r", "p") %in% param_names)) {
      stop("FLPar must contain either 'initial', 'r', 'p' (for recovery time) or 'r', 'p', 'bmsy' (from pt)")
    }
    
    bmsy = if ("bmsy" %in% param_names) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    nits = dims(object)$iter
    result = numeric(nits)
    
    for (i in seq(nits)) {
      result[i] = rTime(
        as.numeric(object["initial", i]),
        r = as.numeric(object["r", i]),
        p = as.numeric(object["p", i]),
        bmsy = as.numeric(bmsy[i]),
        F = F
      )
    }
    
    return(result)
  })



#' @rdname rTime
#' @export
setMethod("rTime", signature(object="FLBRP"),
  function(object, initial, biomass="ssb", F=0, ...) {
    if (missing(initial)) {
      stop("initial biomass must be provided")
    }
    
    # Extract parameters using pt() function from pt.R
    params = pt(object, biomass = biomass)
    
    # Calculate recovery time
    rTime(params, initial = initial, F = F)
  })

#' Calculate carrying capacity from BMSY and shape parameter
#' 
#' @description Uses the Pella-Tomlinson relationship: BMSY/K = (1/(1+p))^(1/p)
#' This is the inverse of the relationship used in pt.R for consistency.
#' 
#' @param p Shape parameter
#' @param bmsy Biomass at MSY
#' @return Carrying capacity K
#' @export
calculateCarryingCapacity = function(p, bmsy) {
  # Handle vector inputs
  if (length(p) > 1 || length(bmsy) > 1) {
    # Vector case
    result = numeric(max(length(p), length(bmsy)))
    # Use the inverse of the BMSY/K ratio from pt.R: K = BMSY / ((1/(1+p))^(1/p))
    result = bmsy / ((1/(1 + p))^(1/p))
    return(result)
  } else {
    # Scalar case
    # Use the inverse of the BMSY/K ratio from pt.R: K = BMSY / ((1/(1+p))^(1/p))
    bmsy / ((1/(1 + p))^(1/p))
  }
}

#' Calculate recovery time using Pella-Tomlinson dynamics
#' 
#' @param B0 Initial biomass
#' @param B1 Target biomass
#' @param K Carrying capacity
#' @param rEff Effective growth rate
#' @param p Shape parameter
#' @param invalid Logical vector of invalid cases
#' @param naMask Logical vector of NA cases
#' @return Recovery time vector
#' @export
calculateRecoveryTime = function(B0, B1, K, rEff, p, invalid, naMask) {
  t = rep(NA_real_, length(B0))
  
  # Pella-Tomlinson model with fishing mortality
  valid_cases = !invalid & !naMask
  if (any(valid_cases)) {
    idx = which(valid_cases)
    t[idx] = calculatePellatRecoveryTime(
      B0[idx], B1[idx], K[idx], rEff[idx], p[idx]
    )
  }
  
  # Set t to NA for invalid r_eff and NA inputs
  t[invalid | naMask] = NA_real_
  
  return(t)
}

#' Calculate recovery time for Pella-Tomlinson model
#' 
#' @param B0 Initial biomass
#' @param B1 Target biomass
#' @param K Carrying capacity
#' @param rEff Effective growth rate
#' @param p Shape parameter
#' @return Recovery time
#' @export
calculatePellatRecoveryTime = function(B0, B1, K, rEff, p) {
  B0Kp = (B0/K)^p
  B1Kp = (B1/K)^p
  denom0 = 1 - B0Kp
  denom1 = 1 - B1Kp
  
  validDomain = denom0 > 0 & denom1 > 0
  
  t = rep(NA_real_, length(B0))
  t[validDomain] = (1/(rEff[validDomain] * p[validDomain])) * 
    log(denom0[validDomain]/denom1[validDomain])
  
  return(t)
}

#' Calculate rebuild time using Pella-Tomlinson dynamics
#'
#' @description Calculate the time required to rebuild a stock from current biomass to target biomass
#' using Pella-Tomlinson surplus production dynamics.
#'
#' @param object Current biomass or FLPar object with parameters
#' @param target Target biomass (default is BMSY)
#' @param r Intrinsic growth rate (if object is numeric)
#' @param p Shape parameter (if object is numeric)
#' @param k Virgin biomass (if object is numeric)
#' @param bmsy Biomass at MSY (if object is numeric)
#' @param F Fishing mortality during rebuild (default=0)
#' @param ... Additional arguments
#' @return Rebuild time in years
#' @export
setGeneric("rebuildTime", function(object, ...) standardGeneric("rebuildTime"))

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object="numeric"),
  function(object, target = NULL, r, p, k, bmsy = NULL, F = 0, ...) {
    if (missing(r) || missing(p) || missing(k)) {
      stop("r, p, and k must be provided when object is numeric")
    }
    
    if (is.null(target)) {
      if (is.null(bmsy)) {
        stop("Either target or bmsy must be provided")
      }
      target = bmsy
    }
    
    # Use rTime for the calculation
    rTime(object, r = r, p = p, bmsy = target, F = F)
  })



#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object="FLPar"),
  function(object, initial, target = NULL, F = 0, ...) {
    if (missing(initial)) {
      stop("initial biomass must be provided")
    }
    
    # Get parameter names from FLPar
    param_names = dimnames(object)$params
    
    # Check if this is a pt() output (has r, p, bmsy, virgin)
    if (all(c("r", "p", "bmsy") %in% param_names)) {
      if (is.null(target)) {
        target = object["bmsy"]
      }
      
      # Use rTime for the calculation
      rTime(object, initial = initial, F = F)
    } else {
      # Original behavior for FLPar with initial, r, p
      if (!all(c("initial", "r", "p") %in% param_names)) {
        stop("FLPar must contain either 'initial', 'r', 'p' (for recovery time) or 'r', 'p', 'bmsy' (from pt)")
      }
      
      if (is.null(target)) {
        if (!("bmsy" %in% param_names)) {
          stop("Either target or bmsy must be provided")
        }
        target = object["bmsy"]
      }
      
      # Use rTime for the calculation
      rTime(object, F = F)
    }
  })

#' Calculate initial biomass required for recovery to BMSY in given time
#'
#' @description Numerically integrates the Pella-Tomlinson production function to estimate the initial biomass 
#' required for a population to recover to BMSY within a specified number of years, with optional fishing mortality during recovery.
#'
#' @param nyrs Number of years for recovery
#' @param r Intrinsic rate of population increase
#' @param p Shape parameter of the Pella-Tomlinson model
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param F Fishing mortality rate during recovery (default=0)
#' @param ... Additional arguments
#' @return Initial biomass required for recovery
#' @export
#' @examples
#' brebuild(10, r=0.2, p=1.5)
#' brebuild(10, r=0.2, p=1.5, F=0.05)
#' 
#' # With FLPar
#' library(FLCore)
#' params=FLPar(nyrs=10, r=0.2, p=1.5, bmsy=1)
#' brebuild(params, F=0.05)
#' 
#' # With FLBRP object
#' data(ple4brp)
#' brebuild(ple4brp, nyrs=10, F=0.05)
setGeneric("brebuild", function(object, ...) standardGeneric("brebuild"))

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="numeric"),
  function(object, r, p, bmsy=1, F=0, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    nyrs = object
    
    # Calculate carrying capacity (virgin biomass) from bmsy and p
    K = calculateCarryingCapacity(p, bmsy)
    
    # Effective growth rate under fishing mortality
    rEff = r - F
    
    # Check if recovery is possible under current fishing mortality
    if (rEff <= 0) {
      stop("Fishing mortality too high: recovery not possible")
    }
    
    # Target biomass
    B1 = bmsy
    
    # Pella-Tomlinson model with fishing mortality
    # Solve for B0 in t = (1/(r_eff * p)) * log((1-(B0/K)^p)/(1-(B1/K)^p))
    # Rearranging: B0 = K * (1 - (1-(B1/K)^p) * exp(r_eff * p * t))^(1/p)
    B1Kp = (B1/K)^p
    B0 = K * (1 - (1 - B1Kp) * exp(rEff * p * nyrs))^(1/p)
    
    return(as.numeric(B0))
  })

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="FLPar"),
  function(object, F=0, ...) {
    if (!all(c("nyrs", "r", "p") %in% names(object)))
      stop("FLPar must contain 'nyrs', 'r', and 'p' parameters")
    
    bmsy=if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    nits=dims(object)$iter
    result=numeric(nits)
    
    for (i in 1:nits) {
      result[i]=brebuild(as.numeric(object["nyrs", i]), 
                           as.numeric(object["r", i]), 
                           as.numeric(object["p", i]), 
                           as.numeric(bmsy[i]), 
                           F=F)
    }
    
    return(FLPar(result, dimnames=list(iter=1:nits)))
  })

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="FLBRP"),
  function(object, nyrs, F=0, biomass="ssb", ...) {
    # Extract Pella-Tomlinson parameters from FLBRP object using pt() function from pt.R
    params = pt(object, biomass=biomass, ...)
    
    if (is.null(params)) {
      stop("Could not extract Pella-Tomlinson parameters from FLBRP object")
    }
    
    # Extract parameters - pt() returns normalized BMSY directly
    r = params["r"]
    p = params["p"]
    bmsy = params["bmsy"]  # pt() already returns normalized BMSY
    
    # Calculate initial biomass for given recovery time
    return(brebuild(nyrs, r=r, p=p, bmsy=bmsy, F=F, ...))
  })

#' Calculate generation time from r and R0
#'
#' @description Calculates generation time using the intrinsic rate of increase and the net reproductive rate.
#'
#' @param R0 Net reproductive rate (mean number of offspring per individual per generation)
#' @param r Intrinsic rate of population increase
#' @return Generation time in years
#' @export
#' @examples
#' gt(R0=2, r=0.2)
gtR0=function(R0, r) {
  return(log(R0) / r)
}

#' Calculate generation time from a life table
#'
#' @description Calculates the mean generation time from age-specific survivorship and fecundity.
#'
#' @param ages Vector of ages
#' @param lx Vector of age-specific survivorship (proportion surviving to each age)
#' @param mx Vector of age-specific fecundity (mean number of offspring at each age)
#' @return Generation time in years
#' @export
#' @examples
#' ages<-0:5
#' lx<-c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
#' mx<-c(0, 0, 2, 3, 1, 0)
#' gtLife(ages, lx, mx)
gtLife=function(ages, lx, mx) {
  numerator=sum(ages * lx * mx)
  denominator=sum(lx * mx)
  return(numerator / denominator)
}


