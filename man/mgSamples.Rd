\name{mgSamples-method}
\alias{mgSamples}
\alias{mgSamples,mgNetwork-method}

\title{Samples from Network}

\description{Generates samples from of a \code{mgNetwork} object.}

\usage{
mgSamples(object, numsamples = 1, output="frame")
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{numsamples}{an \code{integer}, the number of samples to be generated}
\item{output}{a \code{character}, the output format. Can be a \code{data.frame} or \code{matrix}.}
}

\details{If the output format is "matrix" then the resulting sample matrix is in row-node format - the rows correspond to the object's nodes while 
the individual samples are represented by columns. 
If the output format is "frame", which is by default, 
the result is a data frame with columns representing the nodes.
}

\value{A \code{matrix} or \code{data.frame} of \code{numeric} values}

\author{N. Balov}

\seealso{\code{\link{mgPredict}} }

\keyword{methods}
\keyword{graphs}
