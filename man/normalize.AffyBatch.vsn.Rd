\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, subsample=20000, niter = 4, ...)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{subsample}{The number of probe to be sampled for the fitting of the transformation parameters. 
    If you want all probes to be used, set this argument to NA. Be prepared for a long computation time.}
  \item{niter}{Parameter passed on to \code{\link{vsn}}.}
  \item{...}{Further parameters for \code{\link{vsn}}.}
}

\details{
  Processing of Affymetrix cell intensities involves the following steps:
  calibration, transformation, probe set summary. \code{\link{vsn}} performs the first 
  two of these, but not the probe set summary.

  Since the function \code{\link{medianpolish}} from the library
  \code{affy} always take the logarithm before 

  The package vsn provides a function \code{\link{truemedianpolish}}

  Please note that there is a , please make sure that it does 
  not take the logarithm, since the result of vsn is already on a transformed scale!

  The cooperation of \code{\link{vsn}} with the package \code{affy} is still experimental. 
  If you have questions, please contact the author.
}

\value{  
An object of class \code{AffyBatch}. 
}
\author{Wolfgang Huber \url{mailto:w.huber@dkfz.de}}

\seealso{\code{\link{vsn}}}
\examples{
  library(affy)

  normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods, "vsn")

  data(Dilution)

  es1 = express(Dilution,
        normalize.method = "vsn", niter=4,
        bg.correct       = FALSE,
        summary.stat     = truemedianpolish)

  es2 = express(Dilution)

  if(interactive()) {
    x11(width=9, height=3.5)
    par(mfrow=c(1,3))
    pch="."
    arrays = c(1,4)
    plot(exprs(Dilution)[, arrays], pch=pch, log="xy", main="raw")
    plot(exprs(es1)[, arrays]/log(2), pch=pch, main="vsn")
    plot(exprs(es2)[, arrays],        pch=pch, main="rma")
  }
}

\keyword{robust}
