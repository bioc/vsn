\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, subsample=20000, ...)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{subsample}{The number of probe to be sampled for the fitting of the transformation parameters. 
    If you want all probes to be used, set this argument to NA. Be prepared for a long computation time.}
  \item{...}{Parameters for \code{\link{vsn}} (see example).}
}

\details{
  Processing of Affymetrix cell intensities involves the following steps:
  calibration, transformation, probe set summary. \code{\link{vsn}} performs the first 
  two of these, but not the probe set summary. Use one of the methods provided
  in the affy package. The example below uses \code{avgdiff}. If you want to use
  another summary method, like \code{medianpolish}, please make sure that it does 
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

  data(affybatch.example)
  es = expresso(affybatch.example,
  	pmcorrect.method = "pmonly",
        bgcorrect.method = "none",
        normalize.method = "vsn",
        normalize.param  = list(lts.quantile=0.5, niter=4),
        summary.method   = "avgdiff", verbose = TRUE)

  plot(exprs(es))
}

\keyword{robust}
