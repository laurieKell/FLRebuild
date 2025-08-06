#' Pella-Tomlinson Utility Functions
#' 
#' @description Common utility functions for Pella-Tomlinson production model calculations.
#' These functions are used across multiple pellat* files to avoid code duplication.

#' Extract reference points from FLBRP object
#' 
#' @param x FLBRP object
#' @param biomass Character string specifying biomass type ("ssb" or "eb")
#' @return Named vector with reference points
#' @export
refs = function(x, biomass = "ssb") {
  res = refptsEB(x)[c("msy", "crash", "virgin"),
                    c("harvest", "yield", "ssb", "eb"), drop = TRUE]
  
  c(msy = res[1, "yield"],
    bmsy = res[1, biomass],
    virgin = res[3, biomass],
    fcrash = res[2, "harvest"],
    fmsy = res[1, "harvest"])
}

#' Calculate BMSY/K ratio from shape parameter p
#' 
#' @param p Shape parameter of Pella-Tomlinson model
#' @return BMSY/K ratio
#' @export
bmsyK = function(p) {
  if (p <= -1) {
    stop("Shape parameter p must be > -1.")
  }
  (1 / (1 + p))^(1 / p)
}

#' Convert shape parameter p to BMSY/K ratio
#' 
#' @param p Shape parameter
#' @return BMSY/K ratio
#' @export
p2shape = function(p) {
  (1 / (1 + p))^(1 / p)
}

#' Robust root-finding function for p given a valid BMSY/K ratio
#' 
#' @param target_ratio BMSY/K ratio (must be between 0 and 1)
#' @param lower Lower bound for p search
#' @param upper Upper bound for p search
#' @param tol Tolerance for root finding
#' @return Shape parameter p
#' @export
bmsyK2p = function(target_ratio, lower = 0.05, upper = 5, tol = 1e-8) {
  if (!is.numeric(target_ratio) || length(target_ratio) != 1 || target_ratio <= 0 || target_ratio >= 1) {
    warning("target_ratio must be a scalar between 0 and 1, got: ", target_ratio)
    return(NA)
  }
  
  f = function(p) {
    if (abs(p) < tol) {
      # Analytic limit as p -> 0
      exp(-1) - target_ratio
    } else if (p <= -1) {
      1000
    } else {
      (1 / (1 + p))^(1 / p) - target_ratio
    }
  }
  
  tryCatch({
    res = uniroot(f, lower = lower, upper = upper, tol = tol)
    p = res$root
    
    # Validate that p is within valid range
    if (p <= -1) {
      warning("Calculated p <= -1 (p = ", p, "), which is invalid for Pella-Tomlinson model")
      return(NA)
    }
    
    # Bound p at 1e-4 to prevent it from getting too close to zero
    if (p < 1e-4) {
      p = 1e-4
    }
    
    return(p)
  }, error = function(e) {
    warning("Failed to find valid p for bmsyK ratio ", target_ratio, ": ", e$message, ". Setting p to lower bound (", lower, ")")
    return(lower)
  })
}

#' Calculate shape parameter p from r and FMSY
#' 
#' @param r Intrinsic growth rate
#' @param Fmsy Fishing mortality at MSY
#' @return Shape parameter p
#' @export
rFmsy2p = function(r, Fmsy) {
  if (!is.numeric(r) || !is.numeric(Fmsy) || length(r) != 1 || length(Fmsy) != 1) {
    warning("r and Fmsy must be single numeric values, got r=", r, ", Fmsy=", Fmsy)
    return(NA)
  }
  if (r <= 0 || Fmsy <= 0) {
    warning("r and Fmsy must be positive for valid p calculation, got r=", r, ", Fmsy=", Fmsy)
    return(NA)
  }
  
  p = (r / Fmsy) - 1
  
  # Validate that p is within valid range
  if (p <= -1) {
    warning("Calculated p <= -1 (p = ", p, "), which is invalid for Pella-Tomlinson model")
    return(NA)
  }
  
  # Bound p at 1e-4 to prevent it from getting too close to zero
  if (p < 1e-4) {
    p = 1e-4
  }
  
  return(p)
}



#' Validate Pella-Tomlinson parameters
#' 
#' @param r Intrinsic growth rate
#' @param p Shape parameter
#' @param virgin Virgin biomass
#' @param bmsy Biomass at MSY
#' @param fmsy Fishing mortality at MSY
#' @param msy Maximum sustainable yield
#' @return Logical indicating if parameters are valid
#' @export
validatePellatParams = function(r = NULL, p = NULL, virgin = NULL, bmsy = NULL, fmsy = NULL, msy = NULL) {
  errors = character()
  
  if (!is.null(r) && (is.na(r) || r <= 0)) {
    errors = c(errors, "r must be positive")
  }
  
  if (!is.null(p) && (is.na(p) || p <= -1)) {
    errors = c(errors, "p must be > -1")
  }
  
  if (!is.null(virgin) && (is.na(virgin) || virgin <= 0)) {
    errors = c(errors, "virgin must be positive")
  }
  
  if (!is.null(bmsy) && (is.na(bmsy) || bmsy <= 0)) {
    errors = c(errors, "bmsy must be positive")
  }
  
  if (!is.null(fmsy) && (is.na(fmsy) || fmsy <= 0)) {
    errors = c(errors, "fmsy must be positive")
  }
  
  if (!is.null(msy) && (is.na(msy) || msy <= 0)) {
    errors = c(errors, "msy must be positive")
  }
  
  # Check logical relationships
  if (!is.null(bmsy) && !is.null(virgin) && bmsy >= virgin) {
    errors = c(errors, "bmsy cannot be greater than or equal to virgin biomass")
  }
  
  if (length(errors) > 0) {
    warning(paste("Parameter validation errors:", paste(errors, collapse = "; ")))
    return(FALSE)
  }
  
  return(TRUE)
} 