#' Calculate Biomass Limit Reference Points (Blim)
#' 
#' @description
#' Calculates biomass limit reference points (Blim) based on a ratio of virgin recruitment
#' for an FLBRP object.
#'
#' @param object An FLBRP object containing stock-recruitment parameters.
#' @param ratio Numeric value between 0 and 1 indicating the fraction of virgin recruitment 
#' to use for Blim calculation (default: 0.3).
#'
#' @details
#' The method:
#' - Extracts virgin recruitment from reference points
#' - Applies the specified ratio to calculate Blim recruitment
#' - Computes associated reference points (harvest, yield, SSB, biomass, etc.)
#'
#' @return
#' An FLPar object containing Blim reference points with the following quantities:
#' - harvest: Fishing mortality at Blim
#' - yield: Catch at Blim
#' - rec: Recruitment at Blim
#' - ssb: Spawning Stock Biomass at Blim
#' - biomass: Total biomass at Blim
#' - revenue: Revenue at Blim
#' - cost: Cost at Blim
#' - profit: Profit at Blim
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' brp <- FLBRP(ple4)
#' blim_ref <- blim(brp, ratio=0.3)
#' }
#'
#' @seealso
#'   \code{\link{FLBRP}} \code{\link{refpts}} \code{\link{computeRefpts}}
#' @importFrom FLBRP refpts computeRefpts 
#' @importFrom FLCore FLPar
#' @export
setMethod("blim", signature(object="FLBRP"),
       function(object, ratio=0.3) {
         # Get virgin recruitment 
         rec = refpts(object)["virgin","rec",drop=TRUE]*ratio
         
         # Create new FLPar with NAs
         refpts(object) = FLPar(NA,
                                dimnames=list(refpt="rec",
                                              quant=c("harvest","yield","rec",
                                                      "ssb","biomass","revenue",
                                                      "cost","profit"),
                                              iter=dimnames(refpts(object))$iter))
         
         # Set recruitment
         refpts(object)[,"rec"] = rec
         
         # Return computed reference points
         rtn = computeRefpts(object)["rec"]
         dimnames(rtn)$refpt = "blim"
         
         return(rtn)})

setMethod("blim", signature(object="FLBRPs"), function(object, ratio=0.3) {
  
  rtn=ldply(object, blim, ratio=ratio)

  names(rtn)[c(2:9)]=c("harvest","yield","rec","ssb","biomass","revenue","cost","profit")
  rtn})
          
#' @seealso
#' \code{\link{FLBRP}} \code{\link{refpts}} \code{\link{computeRefpts}}

#' @rdname msyVirgin
#' @export
setMethod("msyVirgin", signature(object="FLBRP"),
  function(object, ...) {
    # Get reference points
    refs <- refpts(object)
    
    # Extract key metrics
    msy_metrics <- c(
      F_msy = c(refs["msy", "harvest"]),
      SSB_msy = c(refs["msy", "ssb"]),
      catch_msy = c(refs["msy", "yield"]),
      ebiomass_msy = c(ebiomass(object))
    )
    
    virgin_metrics <- c(
      F_virgin = c(refs["virgin", "harvest"]),
      SSB_virgin = c(refs["virgin", "ssb"]),
      catch_virgin = c(refs["virgin", "yield"]),
      ebiomass_virgin = c(ebiomass(object))
    )
    
    # Combine all metrics
    all_metrics <- c(msy_metrics, virgin_metrics)
    
    return(all_metrics)
  })
