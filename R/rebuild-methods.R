#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object="FLBRP"),
  function(object, targetF=computeRefpts(object)["msy","harvest"]*0,
                   targetSSB=computeRefpts(object)["msy","ssb"],
                   nInitial=100, growthRate=0.25, minVal=1e-6, maxVal=1,
                   burnin=20, truncate=TRUE) {
            
    # Input validation
    if (!is(object, "FLBRP"))
      stop("object must be an FLBRP object")
            
    if (!all(sapply(list(nInitial, burnin), function(x) is.numeric(x) && x > 0)))
      stop("nInitial and burnin must be positive integers")
            
    if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric)))
      stop("growthRate, minVal, and maxVal must be numeric")
            
    if (minVal >= maxVal)
      stop("minVal must be less than maxVal")
            
    # Generate SSB sequence
    targetSSB=c(targetSSB) * seq(minVal^growthRate, maxVal^growthRate, 
                                            length.out=nInitial)^(1/growthRate)
    targetF=c(targetF)
            
    # Setup equilibrium
    eql=object
    fbar(eql)[]=0.2
            
    # Create target biomass array
    btar=FLQuant(rep(targetSSB, each=dim(fbar(eql))[2]),
                dimnames=dimnames(propagate(ssb(eql), nInitial)))
        
    # Project stock
    stk=propagate(as(eql, "FLStock"), nInitial)
    stk=fwd(stk, ssb_end=btar[,-seq(dims(stk)[["min"]]+2)], sr=eql)
            
    # Apply target F
    ftar=fbar(stk) %=% targetF
    stk=fwd(stk, f=ftar[,-seq(burnin)], sr=eql)
            
    # Post-process
    if (truncate) 
        stk=stk[,-seq(burnin)]
            
    stk=qapply(stk, function(x) {
        dimnames(x)$year=seq(length(dimnames(x)$year))
        x})
            
    return(stk)})

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object="biodyn"), 
          function(object, target=refpts(object)["bmsy"], nInitial=100, 
                   growthRate=0.3, minVal=1e-6, maxVal=1, nx=101) {
            if (!is(object, "biodyn"))
              stop("object must be a biodyn object")
            
            if (!is.numeric(target) || length(target) != 1)
              stop("target must be a single numeric value")
            
            if (!all(sapply(list(nInitial, nx), function(x) is.numeric(x) && x > 0)))
              stop("nInitial, and nx must be positive integers")
            
            if (!all(sapply(list(growthRate, minVal, maxVal), is.numeric)))
              stop("growthRate, minVal, and maxVal must be numeric")
            
            minVal=1e-6
            if (minVal >= maxVal)
              stop("minVal must be less than maxVal")
          
            bmsy=c(refpts(object)["bmsy"])
               
            rtn=propagate(object, nInitial)
            
            # Create stock projection
            target_seq=c(target)*seq(minVal^growthRate, maxVal^growthRate, length.out=nInitial)^(1/growthRate)
            rtn@stock=FLQuant(rep(target_seq, each=dim(rtn)[2]), dimnames=dimnames(stock(rtn)))
            rtn=fwd(rtn, catch=catch(rtn)[,-1]%=%0.0)

            # Transform data
            dat=as.data.frame(stock(rtn), drop=TRUE)
            dat$initial=c(stock(rtn)[,1])[an(dat$iter)]
            dat=dat[,-2]
            
            # Interpolate results
            dat=as.data.frame(with(dat, akima::interp(initial, 
                        data, year, yo=bmsy, duplicate="mean", nx=nx, jitter=1e-6)))[,c(3,1)]
            names(dat)=c("year", "initial")
            
            transform(dat,initial=initial/bmsy)})

#' @rdname rebuildTime
#' @export
setMethod("rebuildTime", signature(object="FLStock"),
          function(object, nx=101) {
            if (!is(object, "FLStock"))
              stop("object must be an FLStock object")
            
            if (!is.numeric(nx) || nx <= 0)
              stop("nx must be a positive integer")
            
            bmsy=c(ssb(object)[,1,,,,dim(object)[6]])
            
            dat=transmute(as.data.frame(ssb(object), drop=TRUE),
                             ssb=data/bmsy,
                             initial=c(ssb(object[,1]))[an(ac(iter))]/bmsy,
                             year=year)
            
            rtn=suppressWarnings(
              as.data.frame(with(dat, 
                                 akima::interp(x=initial, y=ssb, z=year, yo=1,
                                               duplicate="mean", nx=nx, jitter=1e-6)))[,c(3,1)])
            
            names(rtn)=c("year", "initial")
            return(data.frame(rtn))
          }) 