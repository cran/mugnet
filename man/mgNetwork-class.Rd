\name{mgNetwork-class}
\docType{class}
\alias{mgNetwork}
\alias{mgNetwork-class}

\title{The class "mgNetwork"}

\description{This is the base class implementing Bayesian Networks with categorical values and normal errors}

\section{Objects from the Class}{
  Objects can be created by calls of 

  \code{mgNew(nodes, cats, parents, probs, model, betas, sigmas)}  

}
\section{Slots}{
  \describe{
	\item{\code{objectName}}{an optional object name of class \code{character}.}
	\item{\code{numnodes}:}{an integer, the number of nodes in the object.}

	\item{\code{nodes}:}{a \code{vector} specifying the node names.}
	\item{\code{parents}:}{a \code{list} specifying the node parents. 
		The list \code{parents} must be the same length as \code{nodes}. 
		Parents are kept as indices in the \code{nodes} vector. 
    	}

	\item{\code{categories}:}{a \code{list} of characters specifying a set of categorical values for each node.
    	}
	
	\item{\code{probabilities}:}{a numerical \code{list} that for each node specifies a discrete probability distribution - 
		the distribution of the node conditional on its parent set. 
		The elements of \code{probabilities} are lists themselves. See \code{cnProb} function for more details.
	    }
	
	\item{\code{maxParents}:}{an integer, the maximum number of node parents.}
	\item{\code{maxCategories}:}{an integer, the maximum number of node categories.}
	\item{\code{meta}:}{an object of class \code{character} storing some meta-data information.}

	\item{\code{nodeComplexity}:}{a numerical \code{vector}, the node complexities.} 
	\item{\code{nodeLikelihood}:}{a numerical \code{vector}, the node likelihoods of the sample being used for estimation.}

	\item{\code{complexity}:}{an \code{integer}, the network complexity} 
	\item{\code{likelihood}:}{a \code{numerical}, the total likelihood of the sample being used for estimation}

	\item{\code{model}}{a \code{character}, the parametric node error model; currently: Gaus, Pois, Exp.}
	\item{\code{betas}}{a \code{list}, the beta-parameters for the model.}
	\item{\code{sigmas}}{a \code{vector}, the sigma-parameters for the model.}
  }
}

\section{Methods}{
  \describe{
	\item{mgNew}{\code{signature}(nodes="\code{vector}", 
		cats="\code{list}", parents="\code{list}", probs="\code{list}", model="\code{character}", 	betas="\code{list}", sigmas="\code{vector}"):
		Creating a new class object.}
  }
}

\details{
   The \code{mgNetwork} class provides a general structure for representing discrete Bayesian networks. The class is designed such as to provide suitable tools for statistical inference, in particular, maximum likelihood based network estimation from a sample. 
 }

\author{N. Balov}

\seealso{\code{\link{mgNew}}
\code{\link{mgBeta}}, \code{\link{mgSigma}}, \code{\link{mgSearchOrder}}}

\keyword{classes}
\keyword{graphs}

