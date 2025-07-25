% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\name{FLRebuild}
\alias{FLRebuild}
\alias{FLRebuild-package}
\title{FLRebuild: Fisheries Stock Rebuilding Analysis}
\description{
A comprehensive package for analyzing fisheries stock rebuilding trajectories 
and calculating rebuilding times. Provides methods for projecting stock 
rebuilding from different initial depletion levels using FLR 
(Fisheries Library for R) objects.

The package supports multiple stock assessment models:
\itemize{
  \item \code{FLBRP} objects for equilibrium-based rebuilding analysis
  \item \code{biodyn} objects for biomass dynamic models
  \item \code{FLStock} objects for time series analysis
}

Key features include:
\itemize{
  \item Rebuilding trajectory projections
  \item Recovery time calculations
  \item Age-based indicators (ABI)
  \item Stock-recruitment relationship analysis
  \item Reference point calculations (Blim, MSY)
}
}
\examples{
# Load the package
library(FLRebuild)

# Example with biodyn object
# bd <- biodyn(FLPar(r=0.5, k=1000, p=1))
# rebuild_data <- rebuildTime(bd)

# Example with FLBRP object  
# eq <- lhEql(lhPar(FLPar(linf=250, s=0.9)))
# stk <- rebuild(eq)
}
