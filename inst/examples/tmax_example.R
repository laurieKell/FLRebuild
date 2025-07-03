# Example: Tmax Calculation Methods
# This script demonstrates the three methods for calculating Tmax
# as described in National Standard 1 guidelines

library(FLRebuild)
library(FLife)

# Example 1: Generation Time Method
# Tmax = Tmin + generation time
tmin <- 5  # Minimum rebuilding time
tmax_gen <- calculateTmax(object = NULL, method = "generation", tmin = tmin)
cat("Tmax using generation method:", tmax_gen, "years\n")

# Example 2: Multiply Method  
# Tmax = Tmin * 2
tmax_mult <- calculateTmax(object = NULL, method = "multiply", tmin = tmin)
cat("Tmax using multiply method:", tmax_mult, "years\n")

# Example 3: Recovery Time Calculation
# Calculate time to recovery for different initial biomass levels
initial_levels <- c(0.1, 0.3, 0.5, 0.7)
target_level <- 1.0
growth_rate <- 0.1

recovery_times <- sapply(initial_levels, function(x) {
  calculateRecoveryTime(x, target_level, growth_rate)
})

cat("\nRecovery times for different initial biomass levels:\n")
for(i in seq_along(initial_levels)) {
  cat("Initial biomass:", initial_levels[i], "-> Recovery time:", recovery_times[i], "years\n")
}

# Example 4: Working with FLR objects
# Create a simple life history equilibrium
par <- lhPar(FLPar(linf = 250, s = 0.9))
eq <- lhEql(par)

# Get generation time
gen_time <- gt(eq)
cat("\nGeneration time for this stock:", gen_time, "years\n")

# Calculate Tmax using generation time method with actual object
tmax_with_object <- calculateTmax(object = eq, method = "generation", tmin = tmin)
cat("Tmax with actual object:", tmax_with_object, "years\n") 