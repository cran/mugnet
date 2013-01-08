\name{mgBeta-method}
\alias{mgBeta}
\alias{mgBeta,mgNetwork-method}

\title{Beta Parameters}

\description{Returns the values of the \code{beta} slot of a \code{mgNetwork} object.}

\usage{
mgBeta(object)
}
\arguments{
\item{object}{a \code{mgNetwork} object}
}

\details{Each node has the number of categories beta values, which in the Gaussian mixture model correspond to the mean values for the node categories.}

\value{A list of \code{numeric} vectors.}

\author{N. Balov}
\examples{
cnet <- cnRandomCatnet(numnodes=10, maxParents = 3, numCategories = 2)
nnet <- new("mgNetwork", cnet, "Gaus", NULL, NULL)
mgBeta(nnet)
}
\keyword{methods}
\keyword{graphs}
