\name{mgSetProb-method}
\alias{mgSetProb}
\alias{mgSetProb,mgNetwork-method}

\title{Set Probability from Data}

\description{The function sets the probability structure of a network from data according to the Maximum Likelihood criterion.}

\usage{
 mgSetProb(object, data)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a data \code{matrix} or \code{data.frame}}
}

\details{
The function generates a new probability table for \code{object} and returns an updated \code{mgNetwork}.
The graph structure of the \code{object} is kept unchanged. 
}

\value{\code{mgNetwork}}

\author{N. Balov}

\keyword{methods}
