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

- `blim(object, ratio = 0.3)`: Calculates biomass limit reference points (Blim) for an FLBRP object
- `msyVirgin(object)`: Calculates MSY and virgin state metrics
- `lag(object, refpt = "msy")`: Estimates population dynamics lag by simulating pulse recruitment events
- `fitPellaT(object, method = "mle")`: Fits Pella-Tomlinson surplus production model parameters using multiple methods (MLE, Method of Moments, Direct calculation, Reference points). Automatically determines p parameter bounds based on BMSY/K ratio.

- `fitPellaTDirect(MSY, BMSY, K)`: Direct calculation of Pella-Tomlinson parameters from known values. Automatically determines p parameter bounds based on BMSY/K ratio.

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

## Dependencies

### Required
- **FLCore** (>= 2.6.0): Core FLR functionality
- **FLBRP** (>= 2.5.0): Biological reference points
- **FLife** (>= 2.6.0): Life history parameters
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
