\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method with expresso}
\description{Wrapper for \code{\link{vsn}} to be used as a normalization
  method with the expresso function of the package affy. This is
  deprecated, plus use \code{\link{justvsn}} instead.}
\usage{normalize.AffyBatch.vsn(abatch, subsample=20000, niter = 4, strata, ...)}
\arguments{
  \item{abatch}{An object of type \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{subsample}{The number of probes to be sampled for the fit of the transformation parameters.}
  \item{niter}{Parameter passed on to \code{\link{vsn}}.}
  \item{strata}{Integer vector of length nrow(y). See 
    \code{\link{vsn}} for details.}
  \item{...}{Further parameters for \code{\link{vsn}}.}
}
 
\details{Please refer to the "details" and "references" sections of the man page for
  \code{\link{vsn}} for more details about this method.

  \bold{Important note}: after calling \code{\link{vsn}}, the function
  \code{normalize.AffyBatch.vsn} \bold{exponentiates} the data. This is done in
  order to make the behavior of this function similar to the other
  normalization methods in affy. There, it is assumed that in subsequent
  analysis steps (e.g. in \code{\link[stats]{medpolish}}), the logarithm to
  base 2 needs to be taken.
}

\value{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\seealso{\code{\link{vsn}}}
\examples{
## Please see vignette.
}

\keyword{robust}
