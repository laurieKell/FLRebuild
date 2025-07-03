# FLRebuild: Fisheries Stock Rebuilding Analysis

[![R-CMD-check](https://github.com/username/FLRebuild/workflows/R-CMD-check/badge.svg)](https://github.com/username/FLRebuild/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/FLRebuild)](https://CRAN.R-project.org/package=FLRebuild)

A comprehensive R package for analyzing fisheries stock rebuilding trajectories and calculating rebuilding times. FLRebuild provides methods for projecting stock rebuilding from different initial depletion levels using FLR (Fisheries Library for R) objects.

## Features

- **Multiple Object Support**: Works with FLBRP, biodyn, and FLStock objects
- **Flexible Rebuilding Analysis**: Projects trajectories from various initial depletion levels
- **Tmax Calculation**: Implements different methods for calculating maximum rebuilding time
- **NS1 Guidelines Compliance**: Follows National Standard 1 guidelines for rebuilding analysis
- **Visualization Tools**: Includes plotting functions for rebuilding trajectories
- **Comprehensive Documentation**: Vignettes and examples for all major functions

## Installation

### From CRAN (when available)
```r
install.packages("FLRebuild")
```

### From GitHub
```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install from GitHub
devtools::install_github("username/FLRebuild")
```

### From Source
```r
# Clone the repository and install
git clone https://github.com/username/FLRebuild.git
install.packages("FLRebuild", repos = NULL, type = "source")
```

## Quick Start

```r
# Load the package
library(FLRebuild)

# Create a life history equilibrium object
library(FLife)
eq <- lhEql(lhPar(FLPar(linf=250, s=0.9)))

# Run rebuilding analysis
stk <- rebuild(eq, nInitial = 50)

# Calculate rebuilding times
rebuild_times <- rebuildTime(stk)

# Plot results
plotRebuildTrajectories(as.data.frame(ssb(stk), drop = TRUE))
```

## Main Functions

### `rebuild()`
Projects rebuilding trajectories from different initial SSB levels.

```r
rebuild(object, targetF = NULL, targetSSB = NULL, nInitial = 100, 
        growthRate = 0.25, minVal = 1e-6, maxVal = 1, 
        burnin = 20, truncate = TRUE)
```

### `rebuildTime()`
Calculates rebuilding time for different initial conditions.

```r
rebuildTime(object, nx = 101)
```

### `calculateTmax()`
Calculates maximum rebuilding time using different methods as described in NS1 guidelines.

```r
calculateTmax(object, method = "generation", tmin = NULL, mfmt = NULL)
```

### `calculateRecoveryTime()`
Calculates time to recovery using population growth rate.

```r
calculateRecoveryTime(initial_biomass, target_biomass, growth_rate)
```

### `plotRebuildTrajectories()`
Creates ggplot visualizations of rebuilding trajectories.

```r
plotRebuildTrajectories(rebuild_data, target_line = 1, title = "Rebuilding Trajectories")
```

## ABI and Blim Methods

The package provides age-based index (ABI) methods for analyzing the age structure of stocks and a method for calculating the biomass limit reference point (Blim):

### ABI Methods

- `abiAge(object, ref = "msy", p = 0.9)`: Calculates the reference age for an FLBRP object at which a given proportion (p) of the cumulative stock numbers is reached.
- `abiMsy(object, ref = "msy", p = 0.9)`: Calculates the proportion of stock numbers above the reference age at MSY.
- `abi(object, age, ...)`: Calculates the observed proportion above the reference age for an FLStock object, relative to the MSY reference.

**Example:**
```r
library(FLCore)
library(FLBRP)
data(ple4)
data(ple4brp)
abiAge(ple4brp)
abiMsy(ple4brp)
abi(ple4, ple4brp)
```

### Blim Method

- `blim(object, ratio = 0.3)`: Calculates biomass limit reference points (Blim) for an FLBRP object based on a ratio of virgin recruitment.

**Example:**
```r
library(FLBRP)
data(ple4)
brp <- FLBRP(ple4)
blim_ref <- blim(brp, ratio = 0.3)
```

## Examples

### Basic Rebuilding Analysis

```r
library(FLRebuild)
library(FLife)

# Create equilibrium object
eq <- lhEql(lhPar(FLPar(linf=250, s=0.9)))

# Run rebuilding analysis
stk <- rebuild(eq, 
               targetF = 0,           # No fishing during rebuild
               targetSSB = 1000,      # Target SSB
               nInitial = 100,        # Number of initial conditions
               burnin = 20)           # Burn-in period

# Calculate rebuilding times
rebuild_times <- rebuildTime(stk, nx = 101)

# Plot results
plotRebuildTrajectories(as.data.frame(ssb(stk), drop = TRUE))
```

### Working with biodyn Objects

```r
library(biodyn)

# Create a biodyn object
bd <- biodyn(FLPar(r=0.5, k=1000, p=1))

# Calculate rebuilding time
rebuild_data <- rebuildTime(bd, target = 1000)
```

### Tmax Calculation Methods

The package implements three methods for calculating Tmax as described in NS1 guidelines:

1. **Generation Method**: Tmin + one generation time
2. **Fishing Method**: Time to rebuild at 75% of Maximum Fishing Mortality Threshold (MFMT)
3. **Multiply Method**: Tmin multiplied by 2

```r
# Calculate Tmax using different methods
tmax_gen <- calculateTmax(object = NULL, method = "generation", tmin = 10)
tmax_mult <- calculateTmax(object = NULL, method = "multiply", tmin = 10)
```

## Dependencies

- **FLBRP**: Fisheries Library for R - Biological Reference Points
- **FLife**: Fisheries Library for R - Life History
- **ggplotFL**: ggplot2 extensions for FLR objects
- **plyr**: Tools for splitting, applying and combining data
- **popbio**: Population Biology
- **interp**: Interpolation methods
- **akima**: Interpolation of irregularly and regularly spaced data

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
Author (Year). FLRebuild: Fisheries Stock Rebuilding Analysis. 
R package version 0.1.0. https://github.com/username/FLRebuild
```

## Acknowledgments

This package builds upon the FLR (Fisheries Library for R) framework and implements methods described in National Standard 1 guidelines for fisheries management. 