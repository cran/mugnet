\name{mgFind-method}
\alias{mgFind}
\alias{mgFind,catNetworkEvaluate-method}

\title{Find Network by Complexity}

\description{
	This is a model selection routine that finds a network in a set of networks for a given complexity.
}

\usage{
 mgFind(object, complexity = 0, alpha=0, factor=1)
}

\arguments{
 \item{object}{\code{catNetworkEvaluate}}
 \item{complexity}{an \code{integer}, target complexity}
 \item{alpha}{a \code{character} or \code{numeric}}
 \item{factor}{a \code{numeric}}
}

\details{The complexity must be at least the number of nodes of the networks. If no network with the requested complexity exists in the list, then the one with the closest complexity is returned. 
Alternatively, one can apply some standard model selection with \code{alpha}="BIC" and \code{alpha}=AIC. 
}

\value{A \code{mgNetwork} object.}

\author{N. Balov}

\keyword{methods}
