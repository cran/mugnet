\name{mgSigma-method}
\alias{mgSigma}
\alias{mgSigma,mgNetwork-method}

\title{Sigma Parameters}

\description{Returns the values of the \code{sigma}s of a \code{mgNetwork} object.}

\usage{
mgSigma(object)
}

\arguments{
\item{object}{a \code{mgNetwork} object}
}

\details{In the Gaussian mixture model each node has one sigma parameter corresponding to its standard deviation.}

\value{A \code{numeric} vector.}

\author{N. Balov}
\examples{
cnet <- cnRandomCatnet(numnodes=10, maxParents = 3, numCategories = 2)
nnet <- new("mgNetwork", cnet, "Gaus", NULL, NULL)
mgSigma(nnet)
}
\keyword{methods}
\keyword{graphs}
