\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\alias{normalize.Plob.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, ...)}
\usage{normalize.Plob.vsn(aplob, ...)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{aplob}{A \code{Plob}. }
  \item{...}{Other arguments, currently ignored. \code{\link{vsn}} is
    called with default parameters.}
}

\details{
\code{normalize.AffyBatch.vsn}
is experimental, and has not yet been tested properly.
If you have questions, please contact the author.

\code{normalize.Plob.vsn}
will be obsolete with the release of package affy >1.2.
}

\value{  
An object of of the same class as the first argument
(i.e. \code{AffyBatch} or \code{Plob}). 
}
\author{Wolfgang Huber \url{mailto:w.huber@dkfz.de}}

\seealso{\code{\link{vsn}}}
\examples{
  ## coming soon
}

\keyword{robust}
