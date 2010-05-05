\name{mgFromCatnet}
\alias{mgFromCatnet}
\alias{mgFromCatnet,mgNetwork-method}

\title{mgCatnet from catNetwork}

\description{From existing \code{catNetwork} object and additional model parameters(beta and sigma), the function creates a new \code{mgNetwork} object}

\usage{
mgFromCatnet(object, cnet, betas, sigmas)
}
\arguments{
 	\item{object}{a \code{mgNetwork}}
	\item{cnet}{a \code{catNetwork}}
	\item{betas}{a \code{list} of the beta parameters of the nodes}
	\item{sigmas}{a \code{numeric} vector of the sigma parameters of the nodes}
}

\value{A \code{mgNetwork}.}

\author{N. Balov}

\keyword{methods}
\keyword{graphs}

