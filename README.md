# FLRebuild: Fisheries Stock Rebuilding Analysis

[![R-CMD-check](https://github.com/username/FLRebuild/workflows/R-CMD-check/badge.svg)](https://github.com/username/FLRebuild/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/FLRebuild)](https://CRAN.R-project.org/package=FLRebuild)

A comprehensive R package for analyzing fisheries stock rebuilding trajectories and calculating rebuilding times. FLRebuild provides methods for projecting stock rebuilding from different initial depletion levels using FLR (Fisheries Library for R) objects.

## Features

- **Multiple Object Support**: Works with FLBRP, biodyn, and FLStock objects
- **Flexible Rebuilding Analysis**: Projects trajectories from various initial depletion levels
- **Stock-Recruitment Analysis**: Forward and inverse SRR functions for multiple models
- **Age-Based Indicators**: ABI methods for analyzing stock age structure
- **Reference Point Calculations**: Blim and MSY calculations
- **Tmax Calculation**: Implements different methods for calculating maximum rebuilding time
- **NS1 Guidelines Compliance**: Follows National Standard 1 guidelines for rebuilding analysis
- **Visualization Tools**: Includes plotting functions for rebuilding trajectories
- **Comprehensive Documentation**: Vignettes and examples for all major functions

## Installation

### From GitHub
```r
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install from GitHub
remotes::install_github("lauriekell/FLRebuild")
```

### From Source
```r
# Clone the repository and install
git clone https://github.com/lauriekell/FLRebuild.git
install.packages("FLRebuild", repos = NULL, type = "source")
```

## Quick Start

```r
# Load required packages
library(FLCore)
library(FLBRP)
library(FLife)
library(FLasher)
library(FLRebuild)

# Create a life history equilibrium object
eq = lhEql(lhPar(FLPar(linf = 250, s = 0.9)))

# Run rebuilding analysis
stk = rebuild(eq)

# Plot results
ggplot(ssb(stk)) +
  geom_line(aes(year, data, group = iter))

# Calculate rebuilding times
rT = rebuildTime(stk)
```

## Core Functions

### Rebuilding Analysis

- `rebuild(object, ...)`: Projects rebuilding trajectories from different initial SSB levels
- `rebuildTime(object, ...)`: Calculates rebuilding time for given trajectories
- `calculateRecoveryTime(initial_biomass, target_biomass, growth_rate)`: Calculates time to recovery using population growth rate
- `calculateTmax(object, method, tmin, mfmt)`: Calculates maximum rebuilding time using various methods

### Stock-Recruitment Relationships

- `recHat(object, ssb)`: Predicts recruitment from SSB for FLBRP/FLSR objects
- `rickerRec(params, rec)`: Calculates SSB for given recruitment (Ricker model)
- `bevertonRec(params, rec)`: Calculates SSB for given recruitment (Beverton-Holt model)
- `segRegRec(params, rec)`: Calculates SSB for given recruitment (Segmented Regression model)

### Reference Points

- `blim(object, ratio = 0.3)`: Calculates biomass limit reference points (Blim) for an FLBRP object. **FLBRPs method returns a data frame with one row per iteration.**
- `msyVirgin(object)`: Calculates MSY and virgin state metrics
- `lag(object, refpt = "msy")`: Estimates population dynamics lag by simulating pulse recruitment events
- `refs(object)`: Extracts reference points from FLBRP objects. **FLBRPs method returns a data frame with one row per iteration and columns for all reference point parameters.**
- `ptMle(object, biomass = "ssb", ...)`: Fits Pella-Tomlinson surplus production model parameters using maximum likelihood estimation. Optimized for performance with vectorized calculations and efficient parameter estimation. Automatically determines p parameter bounds based on BMSY/K ratio. **FLBRPs method returns a data frame with one row per iteration and columns for all parameters.**

- `fitPellaTDirect(MSY, BMSY, K)`: Direct calculation of Pella-Tomlinson parameters from known values. Automatically determines p parameter bounds based on BMSY/K ratio.

- `createModifiedPellatProduction(r, p, k, bmsy, msy, biomass)`: Creates a modified Pella-Tomlinson production curve that ensures the right-hand limb has an increasing slope, avoiding the issue where p=0.5 creates a curve that peaks too late and then declines.

- `productionModified(params, biomass)`: Calculates production using the modified Pella-Tomlinson model that guarantees proper curve shape and increasing slopes on the right-hand limb.

- `validateModifiedProductionCurve(biomass, production, bmsy, k)`: Validates that the modified production curve has the correct shape with an increasing slope on the right-hand limb.

### Pella-Tomlinson Surplus Production Model

The package includes a comprehensive, refactored implementation of the Pella-Tomlinson surplus production model with the following key functions:

#### Core Methods
- `pt(object, biomass = "ssb", use_scaling = FALSE, ...)`: Main method for fitting Pella-Tomlinson parameters to FLBRP, FLPar, or FLBRPs objects. **Reference point priority: MSY (1st), FMSY (2nd), BMSY (3rd)**. **Automatically prevents decreasing slopes on the right-hand limb.**
  - **FLPar method**: Now accepts any three of `fmsy`, `fcrash`, `bmsy`, `k` (or `virgin`), `msy` to yield a unique solution. **Note: `virgin` is synonymous with `k`** - if `virgin` is supplied but not `k`, then `virgin` is used in place of `k`.
  - **FLBRPs method**: Returns a data frame with one row per iteration and columns for all parameters (r, p, virgin, bmsy, fmsy, msy, fcrash, scaling).


#### Utility Functions
- `calculatePFromBmsyK(bmsy_k_ratio)`: **Simple function to calculate p parameter from BMSY/K ratio** - provides straightforward calculation of the shape parameter p for any valid BMSY/K ratio

- `calculateGrowthRate(p, bmsy, k, msy)`: Calculate intrinsic growth rate for both Pella-Tomlinson and Fox models
- `calculateTheoreticalMSY(r, p, bmsy, k)`: Calculate theoretical MSY from parameters
- `calculateScalingParameter(bmsy_k, r, p, bmsy, k, msy, use_scaling)`: Calculate scaling parameter for low BMSY/K ratios

#### Decreasing Slope Protection
- **Automatic Prevention**: The package automatically detects and prevents any fits that would create decreasing slopes on the right-hand limb
- **Custom Production Functions**: When problematic parameters are detected, the system automatically switches to custom production functions that guarantee proper curve shape
- **Parameter Validation**: Enhanced validation ensures only biologically realistic production curves are accepted
- **Smart Fallbacks**: For low BMSY/K ratios that commonly cause issues, the system uses proven custom production functions

#### Recovery Time Calculations
- `rTime(object, ...)`: Calculate recovery time to BMSY under Pella-Tomlinson dynamics
- `rebuildTime(object, ...)`: Calculate rebuild time using Pella-Tomlinson dynamics
- `calculateCarryingCapacity(p, bmsy)`: Calculate carrying capacity from BMSY and shape parameter
- `calculateRecoveryTime(B0, B1, K, rEff, p, invalid, naMask)`: Calculate recovery time using Pella-Tomlinson dynamics

#### FLPar Management
- `getParamNames(use_scaling)`: Get parameter names based on scaling option
- `createEmptyFLPar(param_names, dims)`: Create empty FLPar with correct structure
- `processPellatIteration(object, i, use_scaling)`: Process single iteration for pellat
- `fillFLPar(res, params, i, use_scaling)`: Fill FLPar with calculated parameters

### Age-Based Indicators (ABI)

- `abiAge(object, ref = "msy", p = 0.9)`: Calculates reference age for FLBRP object
- `abiMsy(object, ref = "msy", p = 0.9)`: Calculates proportion above reference age at MSY
- `abi(object, age, ...)`: Calculates observed proportion above reference age for FLStock
- `ageMax(object, na.rm = TRUE)`: Finds age with maximum value for each year and iteration in FLQuant (returns 1 quant dimension)
- `ageMaxQ(object, na.rm = TRUE)`: Alternative implementation using qapply for better performance

## Examples

### Basic Rebuilding Analysis
```r
# Create equilibrium object
eq = lhEql(lhPar(FLPar(linf = 250, s = 0.9)))

# Run rebuilding analysis
stk = rebuild(eq, targetF = 0, nInitial = 50)

# Calculate rebuilding times
rT = rebuildTime(stk)

# Plot trajectories
plotRebuildTrajectories(as.data.frame(ssb(stk)))
```

### Stock-Recruitment Analysis
```r
# Create SRR parameters
params = FLPar(a = 2000, b = 0.001)

# Predict recruitment from SSB
rec = recHat(eq, ssb = 1000)

# Calculate SSB for given recruitment (Ricker model)
ssb_vals = rickerRec(params, rec = 500)
```

### Age-Based Indicators
```r
# Calculate reference age
ref_age = abiAge(eq)

# Calculate proportion at MSY
p_msy = abiMsy(eq)

# Calculate observed proportion
p_obs = abi(stk, eq)
```

### Reference Points
```r
# Calculate Blim
blim_val = blim(eq, ratio = 0.3)

# Calculate MSY and virgin metrics
msy_virgin = msyVirgin(eq)
```

### Modified Pella-Tomlinson Production Model
```r
# Set up parameters for problematic case (p = 0.5)
r <- 0.5        # Growth rate
p <- 0.5        # Shape parameter (causes issues in standard model)
k <- 1000       # Virgin biomass
bmsy <- 200     # Target BMSY (BMSY/K = 0.2)
msy <- 50       # Target MSY

# Create biomass range
biomass <- seq(0, k, length.out = 100)

# Calculate production using modified model
modified_prod <- createModifiedPellatProduction(r, p, k, bmsy, msy, biomass)

# Compare with standard model
standard_prod <- r * biomass * (1 - (biomass / k)^p) / p

# Plot comparison
plot(biomass, standard_prod, type = "l", col = "red", lwd = 2,
     xlab = "Biomass", ylab = "Production",
     main = "Pella-Tomlinson Comparison (p = 0.5)")
lines(biomass, modified_prod, col = "blue", lwd = 2)
abline(v = bmsy, lty = 2, col = "green")
legend("topright", legend = c("Standard (problematic)", "Modified (fixed)", "Target BMSY"),
       col = c("red", "blue", "green"), lty = c(1, 1, 2), lwd = c(2, 2, 1))

# Validate the modified curve
is_valid <- validateModifiedProductionCurve(biomass, modified_prod, bmsy, k)
cat("Modified curve is valid:", is_valid, "\n")

# Use with FLPar objects
params <- FLPar(r = r, p = p, virgin = k, bmsy = bmsy, msy = msy)
biomass_flq <- FLQuant(biomass)
modified_result <- productionModified(params, biomass_flq)
```

## Dependencies

### Required
- **FLCore** (>= 2.6.0): Core FLR functionality
- **FLBRP** (>= 2.5.0): Biological reference points
- **FLlife** (>= 2.6.0): Life history parameters
- **ggplotFL** (>= 0.2.0): Plotting for FLR objects
- **data.table** (>= 1.14.0): Fast data manipulation

### Suggested
- **testthat** (>= 3.0.0): Unit testing
- **knitr** (>= 1.30): Vignette building
- **rmarkdown** (>= 2.0): Documentation
- **akima** (>= 0.6.0): Interpolation functions

## Documentation

- **Vignettes**: Run `browseVignettes("FLRebuild")` to view package vignettes
- **Help**: Use `?rebuild`, `?rebuildTime`, etc. for function documentation
- **Examples**: See the examples section in each function's help page

## Contributing

We welcome contributions! Please feel free to:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## Issues

If you encounter any issues or have questions, please:

1. Check the existing issues on GitHub
2. Create a new issue with a clear description
3. Include a minimal reproducible example if possible

## License

This package is licensed under the GPL-3 license. See the [LICENSE](LICENSE) file for details.

## Citation

If you use FLRebuild in your research, please cite it as:

```
FLR Team (2024). FLRebuild: Fisheries Stock Rebuilding Analysis. 
R package version 0.1.5. https://github.com/lauriekell/FLRebuild
```

## Acknowledgments

This package builds upon the FLR (Fisheries Library for R) framework and implements methods described in National Standard 1 guidelines for fisheries management. 
