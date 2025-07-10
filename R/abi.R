#' Calculate Reference Age for FLBRP Object
#'
#' @description Calculates the reference age for an FLBRP object, i.e., the age at which a given proportion (p) of the cumulative stock numbers is reached at a specified reference point (e.g., MSY).
#'
#' @param object An FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#' @return An FLQuant object containing reference ages.
#' @export
setMethod("abiAge",
          signature(object = "FLBRP"),
          function(object, ref = "msy", p = 0.9) {
            fbar(object) = as.FLQuant(computeRefpts(object)[ref, "harvest", drop = TRUE], dimnames = list(iter = seq(dim(object)[6])))
            stk.n = stock.n(object)[-1]
            # Calculate cumulative sum and normalize
            cumN = apply(stk.n, c(2, 6), cumsum)
            totN = quantSums(stk.n)
            # Avoid division by zero
            totN[totN == 0] = NA
            propN = cumN / rep(totN, each = dim(stk.n)[1])
            ages_arr = array(as.numeric(ages(stk.n)), dim = dim(stk.n))
            ages_arr[propN <= p] = NA
            # Find the minimum age above threshold for each year/iter
            ref_age = apply(ages_arr, c(2, 6), function(x) min(c(x, dims(object)$max), na.rm = TRUE))
            # Return as FLQuant with correct dimnames
            FLQuant(ref_age, dimnames = list(year = dimnames(stk.n)$year, iter = dimnames(stk.n)$iter))
          }
)

#' Calculate Proportion Above Reference Age at MSY for FLBRP Object
#'
#' @description Calculates the proportion of stock numbers above the reference age at MSY for an FLBRP object.
#'
#' @param object An FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#' @return An FLQuant object containing the proportion above the reference age at MSY.
#' @export
setMethod("abiMsy",
          signature(object = "FLBRP"),
          function(object, ref = "msy", p = 0.9) {
            fbar(object) = as.FLQuant(computeRefpts(object)[ref, "harvest", drop = TRUE], dimnames = list(iter = seq(dim(object)[6])))
            A = abiAge(object, ref, p)
            stk.n = stock.n(object)[-1]
            # Expand A to match stk.n dimensions (age, year, unit, season, area, iter)
            A_exp = FLCore:::expand(A, year = dimnames(stk.n)$year, iter = dimnames(stk.n)$iter)
            # Create a flag array for ages >= reference age
            flag = FLQuant(NA, dimnames = dimnames(stk.n))
            for (i in seq_along(dim(stk.n)[6])) {
              flag[,,, , ,i] = FLQuant(ages(stk.n) >= c(A_exp[,,, , ,i]), dimnames = dimnames(stk.n)[1:6])
            }
            apply(stk.n * flag, c(2, 6), sum, na.rm = TRUE) / apply(stk.n, c(2, 6), sum, na.rm = TRUE)
          }
)

#' Helper: Proportion Above Reference Age for FLStock
#'
#' @description Calculates the proportion of stock numbers above a given reference age for an FLStock object.
#' @param x An FLStock object.
#' @param A Reference age (FLQuant).
#' @return An FLQuant object with the proportion above the reference age.
#' @keywords internal
abistock = function(x, A) {
  stk.n = stock.n(x)[-1]
  if (dim(x)[6] > 1 & dim(A)[6] == 1)
    A = propagate(A, dim(x)[6])
  # Expand A to match stk.n dimensions
  A_exp = FLCore:::expand(A, year = dimnames(stk.n)$year, iter = dimnames(stk.n)$iter)
  flag = FLQuant(NA, dimnames = dimnames(stk.n))
  for (i in seq_along(dim(stk.n)[6])) {
    flag[,,, , ,i] = FLQuant(ages(stk.n) >= c(A_exp[,,, , ,i]), dimnames = dimnames(stk.n)[1:6])
  }
  apply(stk.n * flag, c(2, 6), sum, na.rm = TRUE) / apply(stk.n, c(2, 6), sum, na.rm = TRUE)
}

#' Calculate Proportion Above Reference Age for FLStock
#'
#' @description Calculates the observed proportion above the reference age for an FLStock object, relative to the MSY reference.
#' @param object An FLStock object.
#' @param age An FLBRP object or FLQuant reference age.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#' @return An FLQuant object with the observed proportion above the reference age.
#' @export
setMethod("abi",
          signature(object = "FLStock", age = "FLBRP"),
          function(object, age, ref = "msy", p = 0.9) {
            pmsy = abiMsy(age, ref, p)
            age = abiAge(age, ref, p)
            pt  = abistock(object, age)
            pt %/% pmsy
          })

#' @rdname abi
#' @export
setMethod("abi",
          signature(object = "FLStock", age = "FLQuant"),
          function(object, age) {
            stk.n = stock.n(object)[-1]
            if (dim(object)[6] > 1 & dim(age)[6] == 1) {
              age = propagate(age, dim(object)[6])
            }
            amsy = FLQuant(rep(c(age), each = prod(dim(stk.n)[-6])), dimnames = dimnames(stk.n))
            flag = FLQuant(ages(stk.n) >= amsy)
            apply(stk.n %*% flag, c(2, 6), sum) %/% apply(stk.n, c(2, 6), sum)
          }
)

#' @examples
#' \dontrun{
#' library(FLCore)
#' library(FLBRP)
#' data(ple4)
#' data(ple4brp)
#' abiAge(ple4brp)
#' abiMsy(ple4brp)
#' abi(ple4, ple4brp)
#' }
