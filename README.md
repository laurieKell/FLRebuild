# FLRebuild: Rebuilding Analysis Functions

Consolidated R package containing functions from FLRebuild and rebuild packages needed for fisheries stock rebuilding analysis workflows.

## Installation

`FLRebuild` depends on [FLR](https://flrproject.org/) packages that are not on CRAN. Install those from the FLR R universe, then install this package. On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first so the C++ code (Rcpp and TMB) can compile.

``` r
repos <- c(
  FLR  = "https://flr.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
)

install.packages(
  c("FLCore", "FLBRP", "FLife", "FLasher", "ggplotFL", "mpb"),
  repos = repos
)

install.packages("remotes")
remotes::install_github(
  "laurieKell/FLRebuild",
  dependencies = TRUE,
  repos = repos
)
```

From a local clone of this repository:

``` r
remotes::install_local("path/to/FLRebuild", dependencies = TRUE, repos = repos)
```

Load the package with:

``` r
library(FLRebuild)
```

## Purpose

This package consolidates all functions needed for rebuilding analysis, eliminating the need to source multiple files from different packages. All functions previously sourced in analysis scripts are now available through a single package.

## Functions Included

### Core Functions:

- `toLogits()`, `fromLogits()` - Steepness parameter transformations
- `ftmb()`, `ftmb2()`, `ftmb3()` - Stock-recruitment relationship fitting
- `updateRefs()` - Update reference points
- `pellatParams()`, `pellatParamFn()` - Pella-Thompson production model parameters
- `abi()`, `abiAge()`, `abiMsy()` - Age-based indicators
- `processError()` - Process error calculation
- `covarFn()`, `leslieFn()` - Life history and demographic covariates
- `invSRR()`, `refCreate()`, `rmax()`, `rmsy()`, `rvirgin()` - Stock-recruitment reference points

### Rebuilding Analysis:

- Rebuilding trajectory functions (`rebuild()`, `rebuildTime()`, etc.)
- Pella-Tomlinson model functions
- Bootstrap and rebuilding analysis functions

## Usage

``` r
library(FLRebuild)
```

## Breaking Changes

⚠️ **Version 0.1.0 includes breaking changes:** - `from_logits()` → `fromLogits()` - `to_logits()` → `toLogits()`

See `NEWS.md` for migration guide.

## Package Status

This package consolidates functions from FLRebuild and rebuild packages. As of version 0.1.0, this package is self-contained and has no dependency on FLCandy.
