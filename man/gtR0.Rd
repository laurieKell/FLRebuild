% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rebuild-pellaT.R
\name{gtR0}
\alias{gtR0}
\title{Calculate generation time from r and R0}
\usage{
gtR0(R0, r)
}
\arguments{
\item{R0}{Net reproductive rate (mean number of offspring per individual per generation)}

\item{r}{Intrinsic rate of population increase}
}
\value{
Generation time in years
}
\description{
Calculates generation time using the intrinsic rate of increase and the net reproductive rate.
}
\examples{
gt(R0=2, r=0.2)
}
