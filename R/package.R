#' FLRebuild: Fisheries Stock Rebuilding Analysis
#' 
#' @description A package for analyzing fisheries stock rebuilding trajectories 
#' and calculating rebuilding times. Provides methods for projecting stock 
#' rebuilding from different initial depletion levels using FLR 
#' (Fisheries Library for R) objects.
#' 
#' @docType package
#' @name FLRebuild
#' @aliases FLRebuild-package
#' 
#' @importFrom methods setGeneric setMethod standardGeneric
#' @importFrom akima interp
#' @importFrom plyr ddply dlply ldply mdply
#' 
#' @examples
#' # Load the package
#' library(FLRebuild)
#' 
#' # Example with biodyn object
#' # bd <- biodyn(FLPar(r=0.5, k=1000, p=1))
#' # rebuild_data <- rebuildTime(bd)
#' 
#' # Example with FLBRP object  
#' # eq <- lhEql(lhPar(FLPar(linf=250, s=0.9)))
#' # stk <- rebuild(eq)
NULL 