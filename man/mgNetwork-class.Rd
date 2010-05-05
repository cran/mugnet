\name{mgNetwork-class}
\docType{class}
\alias{mgNetwork}
\alias{mgNetwork-class}

\title{The class "mgNetwork"}

\description{This is the base class implementing Bayesian Networks with categorical values and normal errors}

\section{Objects from the Class}{
  Objects can be created by calls of 

  \code{mgNew(nodes, cats, parents, probs, betas, sigmas)}  

}
\section{Slots}{
  \describe{
	\item{\code{objectName}}{a \code{character}, an optional object name of class.}
	\item{\code{numnodes}:}{an \code{integer}, the number of nodes.}

	\item{\code{nodes}:}{a \code{character} vector specifying the node names.}
	\item{\code{parents}:}{a \code{list} specifying the node parents. 
		The list \code{parents} must be the same length as \code{nodes}. Parents are given as indices in the \code{nodes} vector. 
    }
	\item{\code{categories}:}{a list of \code{character} vectors specifying the categorical sets of the nodes.
    }
	\item{\code{probabilities}:}{a \code{list} specifying the conditional distributions of the nodes. 
    }
	
	\item{\code{maxParents}}{an \code{integer}, the maximum number of node parents.}
	\item{\code{maxCategories}}{an \code{integer}, the maximum number of node categories.}
	\item{\code{meta}}{a \code{character}, some additional information related to the object.}

	\item{\code{betas}}{a \code{list}, the beta-parameters for the model.}
	\item{\code{sigmas}}{a \code{vector}, the sigma-parameters for the model.}
  }
}

\section{Methods}{
  \describe{
	\item{mgNew}{\code{signature}(nodes="\code{vector}", 
		cats="\code{list}", parents="\code{list}", probs="\code{list}", betas="\code{list}", sigmas="\code{vector}"):
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

