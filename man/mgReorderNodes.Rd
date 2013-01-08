\name{mgReorderNodes-method}
\alias{mgReorderNodes}
\alias{mgReorderNodes,mgNetwork,vector-method}

\title{Reorder Network Nodes}

\description{The function rearranges the nodes of a network according to a new order.}

\usage{
mgReorderNodes(object, nodeIndices)
}

\arguments{
\item{object}{a \code{mgNetwork}}
\item{nodeIndices}{a \code{vector} representing the new node order}
}

\value{A \code{mgNetwork} object.}

\details{Node reordering affects the list of node names, parents and probabilities. 
It is a useful operation in cases when comparison of two networks is needed. }

\author{N. Balov}

\keyword{methods}
\keyword{graphs}
