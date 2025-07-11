# =============================================================================
# Helper Functions
# =============================================================================

#' Interpolate rebuilding times
#' 
#' @param df Data frame with columns: initial, year, ssb
#' @return Data frame with columns: initial, year
#' @keywords internal
interp = function(df) {
  # Input validation
  if (!is.data.frame(df) || !all(c("initial", "year", "ssb") %in% names(df))) {
    stop("df must be a data frame with columns: initial, year, biomass")
  }
  
  # Get unique initial values and sort them
  initials = sort(unique(df$initial))
  Yr = rep(NA, length(initials))
  
  # Find the first year where SSB >= 1 for each initial value
  for (i in seq_along(initials)) {
    sub_df = df[df$initial == initials[i], ]
    idx = which(sub_df$biomass >= 1)
    if (length(idx) > 0) {
      Yr[i] = sub_df$year[min(idx)]
    } else {
      Yr[i] = NA
    }
  }
  
  # Return result
  result = data.frame(
    initial = initials,
    year = Yr
  )
  
  return(result)
}

# =============================================================================
# rebuild Methods
# =============================================================================

#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object = "FLBRP"),
  function(object, targetF = computeRefpts(object)["msy","harvest"] * 0,
           targetSSB = computeRefpts(object)["msy","ssb"],
           nInitial = 100, growthRate = 0.25, minVal = 1e-6, maxVal = 1,
           burnin = 20, truncate = TRUE) {
    
    # Input validation
    if (!is(object, "FLBRP")) {
      stop("object must be an FLBRP object")
    }
    
    if (!all(sapply(list(nInitial, burnin), function(x) is.numeric(x) && x > 0))) {
      stop("nInitial and burnin must be positive integers")
    }
    
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) {
      stop("growthRate, minVal, and maxVal must be numeric")
    }
    
    if (minVal >= maxVal) {
      stop("minVal must be less than maxVal")
    }
    
    # Generate SSB sequence
    targetSSB = c(targetSSB) * seq(minVal^growthRate, maxVal^growthRate, 
                                  length.out = nInitial)^(1/growthRate)
    targetF = c(targetF)
    
    # Setup equilibrium
    eql = object
    fbar(eql)[] = 0.2
    
    # Create target biomass array
    btar = FLQuant(rep(targetSSB, each = dim(fbar(eql))[2]),
                   dimnames = dimnames(propagate(ssb(eql), nInitial)))
    
    # Project stock
    stk = propagate(as(eql, "FLStock"), nInitial)
    stk = fwd(stk, ssb_end = btar[,-seq(dims(stk)[["min"]]+2)], sr = eql)
    
    # Apply target F
    ftar = fbar(stk) %=% targetF
    stk = fwd(stk, f = ftar[,-seq(burnin)], sr = eql)
    
    # Post-process
    if (truncate) {
      stk = stk[,-seq(burnin)]
    }
    
    # Rename years
    stk = qapply(stk, function(x) {
      dimnames(x)$year = seq(length(dimnames(x)$year))
      x
    })
    
    return(stk)
  })

# =============================================================================
# rebuildTime Methods
# =============================================================================

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object = "FLStock"),
  function(object) {

    # Extract SSB data
    df   =FLCore:::as.data.frame(ssb(object), drop = TRUE)
    iters=sort(an(unique(df$iter)))

    # Calculate BMSY and scale SSB
    bmsy      =c(ssb(object)[,1,,,,dim(object)[6]])
    df$biomass=df$data/bmsy
    df$initial=rep(c(ssb(object[,1]))[iters]/bmsy,each=dim(object)[2])
    df = na.omit(df)

    # Interpolate results
    return(interp(df))
  })

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object = "biodyn"), 
  function(object, target = refpts(object)["bmsy"], nInitial = 100, 
           growthRate = 0.3, minVal = 1e-6, maxVal = 1, nx = 101) {
    
    # Input validation
    if (!is(object, "biodyn")) {
      stop("object must be a biodyn object")
    }
    
    if (!is.numeric(target) || length(target) != 1) {
      stop("target must be a single numeric value")
    }
    
    if (!all(sapply(list(nInitial, nx), function(x) is.numeric(x) && x > 0))) {
      stop("nInitial and nx must be positive integers")
    }
    
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric))) {
      stop("growthRate, minVal, and maxVal must be numeric")
    }

    if (minVal >= maxVal) {
      stop("minVal must be less than maxVal")
    }
    
    # Get BMSY
    bmsy = c(refpts(object)["bmsy"])
    
    # Create propagated object
    rtn = propagate(object, nInitial)
    
    # Create stock projection
    target_seq = c(target) * seq(minVal^growthRate, maxVal^growthRate, 
                                length.out = nInitial)^(1/growthRate)
    rtn@stock = FLQuant(rep(target_seq, each = dim(rtn)[2]), 
                       dimnames = dimnames(stock(rtn)))
    rtn = fwd(rtn, catch = catch(rtn)[,-1] %=% 0.0)
    
    # Transform data
    df = as.data.frame(stock(rtn), drop = TRUE)
    df$initial = c(stock(rtn)[,1])[an(df$iter)]
    df = df[,-2]
    names(df) = c("year","biomass","initial")
    
    # Interpolate results
    return(interp(df))
  })

