#' @rdname rebuild
#' @export
# Duplicate setMethod for 'rebuild' with signature(object="FLBRP") removed to avoid conflicts.

#' Calculate rebuilding time
#'
#' @param object An FLStock object
#' @param nx Number of interpolation points (default = 101)
#' @return A data frame with columns year and initial
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
