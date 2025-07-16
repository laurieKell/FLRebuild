# Examples for sr-rs.R methods

library(FLRebuild)
library(FLBRP)
library(FLCore)

# Example FLBRP object (replace with your own if needed)
if (requireNamespace("FLBRP", quietly = TRUE) && requireNamespace("FLCore", quietly = TRUE)) {
  data(ple4brp, package = "FLBRP")

  # Example: Calculate maximum recruitment (rmax)
  max_rec <- rmax(ple4brp)
  cat("Maximum recruitment (rmax):", max_rec, "\n")

  # Example: Calculate recruitment at MSY (rmsy)
  rec_msy <- rmsy(ple4brp)
  cat("Recruitment at MSY (rmsy):", rec_msy, "\n")

  # Example: rmax with ratio argument
  max_rec_80 <- rmax(ple4brp, ratio = 0.8)
  cat("80% of maximum recruitment (rmax):\n")
  print(max_rec_80)
  str(max_rec_80)

  # Example: rmsy with ratio argument
  rec_msy_90 <- rmsy(ple4brp, ratio = 0.9)
  cat("90% of recruitment at MSY (rmsy):\n")
  print(rec_msy_90)
  str(rec_msy_90)

  # Example: Create a reference point FLPar object
  ref_par <- refCreate(ple4brp, "rec", 1000)
  print(ref_par)
  str(ref_par)

  # Example: Inverse SRR (invSRR) for a given recruitment
  rec_val <- FLQuant(1000)
  ssb_val <- invSRR(ple4brp, rec_val)
  cat("SSB for recruitment 1000:", as.numeric(ssb_val), "\n")

  # Example: Helper functions (direct usage)
  params <- FLPar(a = 1000, b = 2)
  rec <- 500
  ssb_bh <- bevholtInv(params, rec)
  ssb_rk <- rickerInv(params, rec)
  cat("Beverton-Holt SSB for rec=500:", ssb_bh, "\n")
  cat("Ricker SSB for rec=500:", ssb_rk, "\n")
} else {
  cat("FLBRP or FLCore not available, skipping examples.\n")
} 