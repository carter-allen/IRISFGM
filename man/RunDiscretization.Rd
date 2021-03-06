% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/Bicluster.R
\name{RunDiscretization}
\alias{RunDiscretization}
\alias{.runDiscretization}
\alias{RunDiscretization,IRISFGM-method}
\title{RunDiscretization}
\usage{
RunDiscretization(object, ...)

.runDiscretization(object = NULL, q = 0.06)

\S4method{RunDiscretization}{IRISFGM}(object = NULL, q = 0.06)
}
\arguments{
\item{object}{input IRIS-FGM object}

\item{...}{other arguments passed to methods}

\item{q}{quantile number which is used as discretized cutoff. The bigger q means more cells will be categorized into 1 in terms of binarizing one gene.}
}
\value{
It will generate quantile based binary matrix.
}
\description{
Run discretization based on Quantile method
}
