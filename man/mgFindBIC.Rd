\name{mgFindBIC-method}
\alias{mgFindBIC}
\alias{mgFindBIC,catNetworkEvaluate-method}

\title{Find Network by BIC}

\description{
	This is a model selection routine that finds a network in a set of networks using the BIC criteria.
}

\usage{
 mgFindBIC(object)
}
\arguments{
 \item{object}{\code{catNetworkEvaluate}}
}

\details{The function returns the network with maximal BIC value from a list of networks 
as obtained from \code{mgSearchOrder}. 
The formula used for the BIC is \code{log(Likelihood) - 0.5*Complexity*log(numNodes)}.}

\value{A \code{mgNetwork} object with optimal BIC value.}

\author{N. Balov}

\keyword{methods}

