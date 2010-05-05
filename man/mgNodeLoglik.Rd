\name{mgNodeLoglik}
\alias{mgNodeLoglik}
\alias{mgNodeLoglik,mgNetwork-method}

\title{Log-Likelihood of Network Nodes}

\description{For a network and some data the function calculates the log-likelihood of the network  nodes.}

\usage{
mgNodeLoglik(object, nodes, data)
}
\arguments{
 	\item{object}{a \code{mgNetwork}}
	\item{nodes}{A \code{vector} of \code{integers}}
	\item{data}{A \code{numeric} data}
}
\details{\code{nodes} is a vector specifying the nodes which log-likelihood is requested.}

\value{A \code{numeric} \code{vector}.}

\author{N. Balov}

\keyword{methods}
\keyword{graphs}

