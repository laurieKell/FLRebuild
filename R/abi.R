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
              # Ensure we have valid ages and reference ages before comparison
              ages_vec = as.numeric(ages(stk.n))
              ref_age = c(A_exp[,,, , ,i])
              if (length(ref_age) == 1) {
                ref_age = rep(ref_age, length(ages_vec))
              }
              # Create flag for ages >= reference age, handling NA values
              age_flag = !is.na(ages_vec) & !is.na(ref_age) & ages_vec >= ref_age
              flag[,,, , ,i] = FLQuant(age_flag, dimnames = dimnames(stk.n)[1:6])
            }
            # Calculate proportions, handling NA values
            numerator = apply(stk.n * flag, c(2, 6), sum, na.rm = TRUE)
            denominator = apply(stk.n, c(2, 6), sum, na.rm = TRUE)
            # Avoid division by zero
            denominator[denominator == 0] = NA
            numerator / denominator
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
    # Ensure we have valid ages and reference ages before comparison
    ages_vec = as.numeric(ages(stk.n))
    ref_age = c(A_exp[,,, , ,i])
    if (length(ref_age) == 1) {
      ref_age = rep(ref_age, length(ages_vec))
    }
    # Create flag for ages >= reference age, handling NA values
    age_flag = !is.na(ages_vec) & !is.na(ref_age) & ages_vec >= ref_age
    flag[,,, , ,i] = FLQuant(age_flag, dimnames = dimnames(stk.n)[1:6])
  }
  # Calculate proportions, handling NA values
  numerator = apply(stk.n * flag, c(2, 6), sum, na.rm = TRUE)
  denominator = apply(stk.n, c(2, 6), sum, na.rm = TRUE)
  # Avoid division by zero
  denominator[denominator == 0] = NA
  numerator / denominator
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
            # Ensure dimensions match before division
            if (dim(pt)[6] == 1 && dim(pmsy)[6] > 1) {
              pt = propagate(pt, dim(pmsy)[6])
            } else if (dim(pmsy)[6] == 1 && dim(pt)[6] > 1) {
              pmsy = propagate(pmsy, dim(pt)[6])
            }
            # Match year dimensions if needed
            if (dim(pt)[2] == 1 && dim(pmsy)[2] > 1) {
              pt = propagate(pt, dim(pmsy)[2])
            } else if (dim(pmsy)[2] == 1 && dim(pt)[2] > 1) {
              pmsy = propagate(pmsy, dim(pt)[2])
            }
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
