% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/LTMG.R
\name{CalBinaryMultiSignal}
\alias{CalBinaryMultiSignal}
\alias{.CalBinaryMultiSignal}
\alias{CalBinaryMultiSignal,IRISFGM-method}
\title{CalBinaryMultiSignal}
\usage{
CalBinaryMultiSignal(object)

.CalBinaryMultiSignal(object = NULL)

\S4method{CalBinaryMultiSignal}{IRISFGM}(object = NULL)
}
\arguments{
\item{object}{Input IRIS-FGM}
}
\value{
It will return a binary matrix based on LTMG signal matrix.
}
\description{
This function is for calculating multisignal from LTMG signaling matrix.
}
\examples{
\dontrun{object <- CalBinaryMultiSignal(object)}
}
