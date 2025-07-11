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



