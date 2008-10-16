\name{justvsn}
\alias{justvsn}
\alias{vsnrma}

\title{Wrapper functions for vsn}
\description{\code{justvsn} returns the vsn-normalised intensities in
  an object generally of the same class as its first
  argument (see the man page of \code{\link[vsn:vsn2trsf]{predict}} for
  details). It preserves the metadata. It is equivalent to calling
\preformatted{
  fit = vsn2(x, ...)
  nx = predict(fit, newdata=x, useDataInFit = TRUE)
}
 \code{vsnrma} is a wrapper around \code{\link{vsn2}} and \code{\link[affy:rma]{rma}}    
}

\usage{
justvsn(x, ...)
vsnrma(x, ...)
}

\arguments{
  \item{x}{For \code{justvsn}, any kind of object for which
    \code{\link{vsn2}} methods exist.
    For \code{vsnrma}, an
    \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{\dots}{Further arguments that get passed on to \code{\link{vsn2}}.}
}

\details{
  \code{vsnrma} does probe-wise
  background correction and between-array normalization by calling
  \code{\link{vsn2}} on the perfect match (PM) values only. Probeset
  summaries are calculated with the medianpolish algorithm of
  \code{\link[affy:rma]{rma}}. 
}

\seealso{
  \code{\link{vsn2}}
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
##--------------------------------------------------
## use "vsn2" to produce a "vsn" object
##--------------------------------------------------
data("kidney")
fit = vsn2(kidney)
nkid = predict(fit, newdata=kidney)

##--------------------------------------------------
## justvsn on ExpressionSet
##--------------------------------------------------
nkid2 = justvsn(kidney)
stopifnot(identical(exprs(nkid), exprs(nkid2)))

##--------------------------------------------------
## justvsn on RGList
##--------------------------------------------------
rg = new("RGList", list(R=exprs(kidney)[,1,drop=FALSE], G=exprs(kidney)[,2,drop=FALSE]))
erge = justvsn(rg)
}


