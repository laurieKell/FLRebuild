interp<-function(df){
  initials=sort(unique(df$initial))
  Yr      =rep(NA, length(initials))
  
  for (i in seq_along(initials)) {
    sub_df = df[df$initial == initials[i], ]
    idx = which(sub_df$ssb >= 1)
    if (length(idx) > 0) {
      Yr[i] = sub_df$year[min(idx)]
    } else {
      Yr[i] = NA}}
  
  result=data.frame(initial=initials,
                    Yr     =Yr)
  result}


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
setMethod("rebuildTime", signature(object="FLStock"),
  function(object, nx=101) {
    if (!is(object, "FLStock"))
      stop("object must be an FLStock object")
    if (!is.numeric(nx) || nx <= 0)
      stop("nx must be a positive integer")

    df <- as.data.frame(ssb(object), drop=TRUE)
    iters <- sort(unique(df$iter))
    # Name the bmsy and initial vectors by iter
    bmsy_vec <- setNames(as.numeric(ssb(object)[,1,,,,iters]), iters)
    initial_vec <- setNames(as.numeric(ssb(object[,1])), iters)
    # Use the iter column to match
    df$ssb <- df$data / bmsy_vec[as.character(df$iter)]
    df$initial <- initial_vec[as.character(df$iter)] / bmsy_vec[as.character(df$iter)]
    df <- na.omit(df)
    return(interp(df))
  })

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
            df=as.data.frame(stock(rtn), drop=TRUE)
            df$initial=c(stock(rtn)[,1])[an(df$iter)]
            df=df[,-2]
            
            # Interpolate results
            #df=as.data.frame(with(df, akima::interp(initial, 
            #            data, year, yo=bmsy, duplicate="mean", nx=nx, jitter=1e-6)))[,c(3,1)]
            #names(df)=c("year", "initial")
            
            #transform(df,initial=initial/bmsy)
            
            return(interp(df))})

