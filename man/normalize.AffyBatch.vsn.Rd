\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, param)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{param}{A list with parameters for \code{\link{vsn}} (see example).}
}

\details{
\code{normalize.AffyBatch.vsn} is experimental.
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
  	pmcorrect.method = "subtractmm",
        bgcorrect.method = "none",
        normalize.method = "vsn",
        normalize.param  = list(lts.quantile=0.55, niter=4),
        summary.method   = "avgdiff", verbose = TRUE)

  plot(exprs(es))
}

\keyword{robust}
