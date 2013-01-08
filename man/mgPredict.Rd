\name{mgPredict}
\alias{mgPredict}
\alias{mgPredict,mgNetwork-method}

\title{Prediction}

\description{Prediction based on a network and partially available data.}

\usage{
mgPredict(object, data)
}
\arguments{
\item{object}{a \code{mgNetwork}}
\item{data}{a \code{numeric} data}
}
\details{\code{Data} should be a \code{matrix} or \code{data.frame} of numerical values. Prediction is made for those nodes that are marked as not-available (\code{NA}) in the data and is based on maximum probability criterion.}

\value{A \code{numeric} \code{matrix} or \code{data.frame}, depending on the input data format.}

\author{N. Balov}

\keyword{methods}
\keyword{graphs}

