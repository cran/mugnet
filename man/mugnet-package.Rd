\name{mugnet-package}
\alias{mugnet-package}
\alias{mugnet}

\docType{package}

\title{mugnet: Mixture of Gaussian Bayesian Network Model}

\description{

\code{mugnet} package provides tools for learning categorical Bayesian networks from continuous data.

A Bayesian network is defined by a graphical structure in form of directed acyclic graph and a probability model given as a set of conditional distributions, one for each node in the network. 
Considered in the package are categorical Bayesian networks with Gaussian error, namely, networks which nodes represent discrete random variables with additional zero-mean Gaussian random terms distributed independently among the nodes.

Implemented in \code{mugnet} is a learning function that for some data and given node order, finds the best fitting networks according to the MLE criterion. 
These optimal networks is believed to explain and represent the relations between the node-variables. 

Since \code{mugnet} inherits the \code{catnet} package, we suggest the user to familiarize with the main features of \code{catnet} - manipulating network attributes, visualizing networks and performing model selection, to mention a few.

}

\author{N. Balov}

