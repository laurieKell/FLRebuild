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
library(plyr)
library(FLCore)
library(FLBRP)
library(FLife)
library(FLasher)
library(FLRebuild)
```

```r
# Create a life history equilibrium object

eq=lhEql(lhPar(FLPar(linf=250, s=0.9)))
```

```r
# Run rebuilding analysis
stk=rebuild(eq)
```

```r
# Calculate rebuilding times
rT=rebuildTime(stk)
```

```r
# Plot results
ggplot(ssb(stk))+
   geom_line(aes(year,data,group=iter))
```

### Blim Method

- `blim(object, ratio = 0.3)`: Calculates biomass limit reference points (Blim) for an FLBRP object based on a ratio of virgin recruitment.

**Example:**
```r
blim(eq, ratio=0.3)
```

## ABI and Blim Methods

The package provides age-based index (ABI) methods for analysing the age structure of stocks and a method for calculating the biomass limit reference point (Blim):

### ABI Methods

- `abiAge(object, ref = "msy", p = 0.9)`: Calculates the reference age for an FLBRP object at which a given proportion (p) of the cumulative stock numbers is reached.
- `abiMsy(object, ref = "msy", p = 0.9)`: Calculates the proportion of stock numbers above the reference age at MSY.
- `abi(object, age, ...)`: Calculates the observed proportion above the reference age for an FLStock object, relative to the MSY reference.


```r
abiAge(eq)
abiMsy(eq)
abi(stk, eq)
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
