\name{mgLoglik-method}
\alias{mgLoglik}

\title{Sample Log-likelihood}

\description{Calculate the log-likelihood of a sample with respect to a given \code{mgNetwork} object}

\usage{
mgLoglik(object, data, bysample=FALSE)
}

\arguments{
\item{object}{a \code{mgNetwork} object}
\item{data}{a data \code{matrix} given in the column-sample format, or a \code{data.frame} in the row-sample format}
\item{bysample}{a \code{logical}}
}

\details{If \code{bysample} is set to TRUE, the function output is a vector of length the sample size with values the log-likelihood of each individual sample instance. Otherwise, the total sum of the log-likelihood of the sample instances is reported.
}

\value{A \code{numeric} value}

\author{N. Balov}

\seealso{\code{\link{mgNodeLoglik}} }

\keyword{methods}
\keyword{graphs}
\keyword{distribution}
