\name{mgNew}
\alias{mgNew}

\title{Creating Random mgNetwork}

\description{Creates a new \code{mgNetwork} with specified nodes, categories, parent sets, probability structure and error model parameters.}

\usage{
mgNew(nodes, cats, parents, probs, model, betas, sigmas)
}

\arguments{
\item{nodes}{a \code{vector} of nodes names}
\item{cats}{a \code{list} of node categories}
\item{parents}{a \code{list} of node parents}
\item{probs}{a \code{list} of probabilities}
\item{model}{a \code{character}, mixture model}
\item{betas}{a \code{list} of beta parameters}
\item{sigmas}{a \code{vector} of sigma parameters}
}

\value{A \code{mgNetwork} object.}

\author{N. Balov}

\examples{
cnet <- mgNew(
	nodes=c("a", "b", "c"),
	cats=list(c("1","2"), c("1","2"), c("1","2")), 
	parents=list(NULL, c(1), c(1,2)), 
	probs=list(	c(0.2,0.8), 
		list(c(0.6,0.4),c(0.4,0.6)), 
		list(list(c(0.3,0.7),c(0.7,0.3)), 
		list(c(0.9,0.1),c(0.1,0.9)))),
	model="Gaus", 
	betas=list(c(0,1), c(0,1), c(0,1)),  
	sigmas=c(1, 1, 1))
}
\keyword{classes}
\keyword{graphs}
