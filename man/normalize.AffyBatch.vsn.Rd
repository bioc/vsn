\name{normalize.AffyBatch.vsn}
\alias{normalize.AffyBatch.vsn}
\title{Wrapper for vsn to be used as a normalization method in the package affy}
\description{Wrapper for \code{\link{vsn}} to be used as a normalization method in the package affy}
\usage{normalize.AffyBatch.vsn(abatch, subsample=20000, niter = 4, ...)}
\arguments{
  \item{abatch}{An object of type \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{subsample}{The number of probes to be sampled for the fit of the transformation parameters.}
  \item{niter}{Parameter passed on to \code{\link{vsn}}.}
  \item{...}{Further parameters for \code{\link{vsn}}.}
}

\details{Please refer to the "details" and "references" sections of the man page for
  \code{\link{vsn}} for more details about this method.

  \bold{Important note}: after calling \code{\link{vsn}}, the function
  \code{normalize.AffyBatch.vsn} \bold{exponentiates} the data. This is done in
  order to make the behavior of this function similar to the other
  normalization methods in affy. There, it is assumed that in subsequent
  analysis steps (e.g. in \code{\link{medianpolish}}), the logarithm to
  base 2 needs to be taken.
}

\value{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
\author{Wolfgang Huber \url{mailto:w.huber@dkfz.de}}

\seealso{\code{\link{vsn}}}
\examples{
library(affy)
library(affydata)
data(Dilution)

## let affy know about vsn
normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods, "vsn")

es1 = expresso(Dilution[1:2], 
               bg.correct       = FALSE,   ## bg correction is done by vsn
               normalize.method = "vsn",
               pmcorrect.method = "pmonly", 
               summary.method   = "medianpolish")

es2 = expresso(Dilution[1:2], 
               bgcorrect.method = "rma",
               normalize.method = "quantiles", 
               pmcorrect.method = "pmonly",
               summary.method   = "medianpolish")

## graphics output
if(interactive()) x11()
oldpar = par(mfrow=c(2,2), pch=".")

## extract expression values
x1 = exprs(es1)
x2 = exprs(es2) 
 
## scatter plot
plot(x1, main="vsn: chip 3 vs 4")
plot(x2, main="rma: chip 3 vs 4")

## rank(mean) - difference plot
ylim = c(-0.7, 0.7)
plot(rank(rowSums(x1)), diff(t(x1)), ylim=ylim, main="rank(mean) vs differences")
abline(h=0, col="red")

plot(rank(rowSums(x2)), diff(t(x2)), ylim=ylim, main="rank(mean) vs differences")
abline(h=0, col="red")

## reset old plotting parameters 
par(oldpar)
}

\keyword{robust}
