#' Rebuild Class
#'
#' @name Rebuild-class
#' @docType class
#' @slot params An \code{FLPar} object containing model parameters
#' @slot nyrs Numeric, number of years for projection
#' @slot niters Numeric, number of iterations
#'
#' @exportClass Rebuild
setClass("Rebuild",
         slots=c(params="FLPar",
                 nyrs="numeric",
                 niters="numeric"))
#' Rebuild Analysis
#'
#' @name rebuild
#' @description Performs rebuilding analysis for fish populations
#'
#' @param r Numeric, intrinsic rate of increase
#' @param p Numeric, shape parameter
#' @param k Numeric, carrying capacity (default: 1e3)
#' @param b0 Numeric, initial biomass ratio (default: 1)
#' @param nyrs Numeric, number of projection years (default: 50)
#' @param niters Numeric, number of iterations (default: 101)
#'
#' @return A data frame containing:
#'   \item{shape}{Shape parameter at BMSY}
#'   \item{year}{Year when target is reached}
#'   \item{initial}{Initial biomass ratio}
#'
#' @importFrom plyr ddply
#' @importFrom mpb biodyn
#' @importFrom FLCore window propagate fwd FLQuant params refpts stock catch
#'
#' @examples
#' \dontrun{
#' res <- rebuild(r=0.5, p=3, k=1000, b0=1)
#' }
#'
#' @export

#' Generation Time
#'
#' @name gt
#' @description Calculates generation time as 1/r
#'
#' @param x An object of class \code{Rebuild}
#'
#' @return Numeric value representing generation time
#'
#' @export
setGeneric("gt", function(x) standardGeneric("gt"))

#' @rdname gt
setMethod("gt", signature(x="Rebuild"),
          function(x) 1/params(x)["r"])

#' @rdname rebuild
setMethod("rebuild", signature(object="numeric", p="numeric"),
          function(object, p, k=1e3, b0=1, nyrs=50, niters=101) {
            r <- object
            object=biodyn(params=FLPar(r=r,p=p,k=k,B0=b0))
            shape=c(refpts(object)["bmsy"]%/%params(object)["k"])
            object=window(object,end=nyrs)
            object@stock[]=refpts(object)["bmsy"]
            object@catch[]=0
            target=c(refpts(object)["bmsy",1])
            object=propagate(object,niters)   
            object@stock=object@stock%*%FLQuant(rep(seq(0,1,length.out=niters),
                                                    each=nyrs),
                                                dimnames=dimnames(stock(object)))
            object=fwd(object,harvest=stock(object)[,-1]%=%0)
            dat=cbind(target=target,as.data.frame(stock(object),drop=TRUE))
            dat=ddply(dat,.(iter), with, {
              rtn=try(data.frame(year=year[(data-target)^2==min((data-target)^2)][1]))
              if ("try-error"%in%is(rtn)) return(NULL)
              return(rtn)
            })
            dat=cbind(dat,initial=c(stock(object)[,1,,,,dat$iter])/target)
            dat=dat[order(dat$initial),c("year","initial")][-1,]
            cbind(shape=shape,dat[dat$year<nyrs,])
          })


