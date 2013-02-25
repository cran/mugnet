\name{mgFindAIC-method}
\alias{mgFindAIC}
\alias{mgFindAIC,catNetworkEvaluate-method}

\title{Find Network by AIC}

\description{
	This is a model selection routine that finds a network in a set of networks using the AIC criteria.
}

\usage{
 mgFindAIC(object)
}
\arguments{
 \item{object}{\code{catNetworkEvaluate}}
}

\details{The function returns the network with maximal AIC value from a list of networks 
as return by \code{mgSearchOrder}. The formula used for the AIC is \code{log(Likelihood) - Complexity}.}

\value{A \code{mgNetwork} object with optimal AIC value.}

\author{N. Balov}

\keyword{methods}
