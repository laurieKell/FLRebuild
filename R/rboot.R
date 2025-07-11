# years {{{

setGeneric("years", function(object, ...)
  standardGeneric("years"))

#' @examples
#' years(m(ple4))
#' # Seasonal objects get decimal years
#' years(expand(m(ple4), season=1:4))[,,,1]
#' years(expand(m(ple4), season=1:4))[,,,2]

setMethod("years", signature(object="FLQuant"),
          function(object) {
            
            res=FLQuant(an(dimnames(object)$year), dimnames=dimnames(object),
                           units="")
            
            # TODO: DEAL with (spawning) units + seasons (e.g. SKJ)
            
            # SEASONAL object
            if(dim(res)[4] > 1) {
              
              nseas=dim(object)[4]
              
              seas=do.call(expand, c(list(x=FLQuant(seq(0, length=nseas, by=1 / nseas),
                                                       dimnames=list(season=dimnames(object)$season), quant=quant(object), 
                                                       units="")), dimnames(object)[-4]))
              
              res=res + seas
            }
            
            return(res)
          }
)
# }}}


setGeneric("rboot",
  function(object, n, len = 3, ...) standardGeneric("rboot"))

setMethod("rboot",
  signature(object = "FLQuant"),
  function(object, n, len = 3, ...) {
    # Number of blocks to sample
    num=ceiling(dim(object)[2] * n / len)
    # Sample starting indices for blocks
    its=sample(seq(dim(object)[2] - len + 1), num, replace = TRUE)
    # Expand to full block indices
    its=c(sapply(its, function(i) i + seq(0, len - 1)))
    # Set up new dimnames, updating 'iter'
    dmns=dimnames(object)
    dmns$iter=seq_len(n)
    # Create new FLQuant with resampled data
    rtn=FLQuant(c(object[, its]), dimnames = dmns)
    return(rtn)})

setMethod("rboot",
          signature(object = "FLQuant"),
          function(object, n, len = 3, ...) {
            # Number of blocks to sample
            num=ceiling(dim(object)[2] * n / len)
            # Sample starting indices for blocks
            its=sample(seq(dim(object)[2] - len + 1), num, replace = TRUE)
            # Expand to full block indices
            its=c(sapply(its, function(i) i + seq(0, len - 1)))
            # Set up new dimnames, updating 'iter'
            dmns=dimnames(object)
            dmns$iter=seq_len(n)
            # Create new FLQuant with resampled data
            rtn=FLQuant(c(object[, its]), dimnames = dmns)
            return(rtn)})

#' Resample FLQuant values from specified years
#'
#' @param object An FLQuant object
#' @param years A vector of years to sample from (must be in dimnames(object)$year)
#' @param n Number of resampled iterations to return
#' @return An FLQuant with resampled values in the iter dimension
#' @export
rboot_years<-function(object, years, n = 100) {
  # Check that years are in the FLQuant
  all_years=dimnames(object)$year
  if (!all(years %in% all_years)) {
    stop("All specified years must be present in the FLQuant's year dimension.")
  }
  # Indices of years to sample from
  year_idx=match(years, all_years)
  # Prepare output dimnames
  dmns=dimnames(object)
  dmns$iter=seq_len(n)
  dmns$year = years
  dmns$iter = as.character(seq_len(n))
  out =FLQuant(NA,dimnames=dmns)
  
  # For each iteration, sample years with replacement
  for (i in seq_len(n)) {
    sampled_idx=sample(year_idx, length(years), replace = TRUE)
    out[,,,,,i]=object[, sampled_idx, , , , drop = FALSE]
  }
  # Return as FLQuant
  FLQuant(c(out), dimnames=dmns)
}


