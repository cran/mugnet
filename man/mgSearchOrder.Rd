\name{mgSearchOrder}
\alias{mgSearchOrder}

\title{Fitting Networks to Data for Given Node Order}

\description{For some data, assumed a sample from an unknown categorical network with Gaussian noise, and an order of the node-variables, the function returns a list of  \code{mgNetwork} objects with complexities up to some maximal value that best fit the data.}

\usage{
mgSearchOrder(data, perturbations, 
		maxCategories, nodeCategories,
		maxParentSet, maxComplexity, 
		nodeOrder, 
		parentsPool=NULL, fixedParentsPool=NULL, 
		emIter=5,
		stopDelta=0,
		emStartIter=2,
		selectMode="BIC",
		model="Gaus",
		echo=FALSE)
}
\arguments{
 	\item{data}{A \code{numeric} data}
	\item{perturbations}{A binary matrix with the dimensions of \code{data}. A value 1 marks that the node in the corresponding sample as perturbed}
	\item{maxCategories}{an \code{integer}, the maximum number of node categories}
	\item{nodeCategories}{an \code{integer} vector, the number of categories for each node}
	\item{maxParentSet}{an \code{integer}, the maximal number of parents per node}
	\item{maxComplexity}{an \code{integer}, the maximal network complexity for the search}
	\item{nodeOrder}{a \code{vector} specifying a node order. The search is among the networks consistent with this topological order}
	\item{parentsPool}{a \code{list} of node parents to choose from}
	\item{fixedParentsPool}{a \code{list} of mandatory node parents}
	\item{emIter}{an \code{integer}, number of iterations for the EM algorithm. The default value is 5.}
	\item{stopDelta}{a \code{numeric}, convergence stopping threshold}
	\item{emStartIter}{an \code{integer}, the number of EM iterations to be performed before checking for convergence}
	\item{selectMode}{a \code{character}, network selection criterion such as "AIC" and "BIC"}
	\item{model}{a \code{character}, network node model such as "Gaus", "Pois" and "Exp"}
	\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{The data is a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format. 
The function returns a list of networks, one for any possible complexity within the specified range that is the MLE in the class of networks having the given topological order of the nodes and complexity.}

\value{A list of \code{mgNetwork} objects.}

\author{N. Balov}
\examples{
cnet <- cnRandomCatnet(numnodes=12, maxParents=2, numCategories=3)
psamples <- cnSamples(cnet, 1000, output="matrix", as.index=TRUE)
fsamples <- psamples + matrix(rnorm(length(psamples), 0.0, 0.05), nrow=dim(psamples)[1])
nodeOrder <- cnOrder(cnet)
nodecats <- rep(2, 12)
nodecats[1] <- 3
nets <- mgSearchOrder(fsamples, NULL, 3, nodecats, 2, 0, nodeOrder)
## next we find the resulting network with complexity of the original one and plot it
## cc <- cnComplexity(cnet) + 12*(1+2)
## cnPlot(cnFind(nets, cc))
}

\keyword{methods}
\keyword{graphs}

