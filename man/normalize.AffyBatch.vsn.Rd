\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, ...)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{...}{Other arguments, currently ignored. \code{\link{vsn}} is
    called with default parameters.}
}

\details{
\code{normalize.AffyBatch.vsn}
is experimental, and has not yet been tested properly.
If you have questions, please contact the author.
}

\value{  
An object of class \code{AffyBatch}. 
}
\author{Wolfgang Huber \url{mailto:w.huber@dkfz.de}}

\seealso{\code{\link{vsn}}}
\examples{
  ## coming soon
}

\keyword{robust}
