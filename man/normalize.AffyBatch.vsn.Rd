\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for vsn to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, ...)}
\arguments{
  \item{abatch}{An \code{AffyBatch} (see documentation of package \code{affy}). }
  \item{...}{Parameters for \code{\link{vsn}} (see example).}
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
  	pmcorrect.method = "pmonly",
        bgcorrect.method = "none",
        normalize.method = "vsn",
        normalize.param  = list(lts.quantile=0.5, niter=4),
        summary.method   = "avgdiff", verbose = TRUE)

  ## note: if you want to use another summary method, please make 
  ## sure that it does not take the logarithm!

  plot(exprs(es))
}

\keyword{robust}
