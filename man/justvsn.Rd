\name{justvsn}
\alias{justvsn}
\alias{justvsn-methods}
\alias{justvsn,ExpressionSet-method}
\alias{justvsn,AffyBatch-method}
\alias{justvsn,RGList-method}

\title{Normalization with vsn}
\description{\code{justvsn} fits the vsn model to data
  and returns an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
  It is a simple wrapper for
  \code{\link{vsn2}} and \code{\link[vsn:vsn2trsf]{predict}}.
}

\usage{
justvsn(x, reference, strata, ...)

\S4method{justvsn}{RGList}(x, reference, strata, backgroundsubtract=FALSE, ...)
}

\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted. Methods exist for \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package) and
    \code{\link[limma:rglist]{RGList}} (from the \code{limma} package).}
  \item{reference, strata, ...}{As in \code{\link{vsn2}}.}
  \item{backgroundsubtract}{Logical of length 1, it applies only to the
    method for \code{\link[limma:rglist]{RGList}} objects and indicates whether
    local background estimates in the \code{Rb} and \code{Gb} slots of
    \code{x} should be subtracted before fitting vsn.}
}

\details{
  See \code{\link{vsn2}}.
}

\value{
  An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
  
  If \code{x} is an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, then the
  dimensions and the metadata of the returned object will be the
  same, and its matrix of expression values are the normalized and
  transformed version of the input.

  If \code{x} is an \code{\link[affy:AffyBatch-class]{AffyBatch}}, probe-wise
  background correction and between-array normalization are done by
  \code{\link{vsn2}} on the perfect match (PM) values only. Probeset
  summaries are calculated with the medianpolish algorithm of
  \code{\link[affy:rma]{rma}}.
  
  If \code{x} is an \code{\link[limma:rglist]{RGList}}, then
  \code{\link{vsn2}} is applied to the matrix \code{cbind(x$R, x$G)}.
  If \code{backgroundsubtract=TRUE}, then \code{cbind(x$Rb, x$Gb)} is
  subtracted before. 
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
## justvsn on AffyBatch
##--------------------------------------------------
data("affybatch.example")
data("cdfenv.example")
eabe = justvsn(affybatch.example)

##--------------------------------------------------
## justvsn on RGList
##--------------------------------------------------
rg = new("RGList", list(R=exprs(kidney)[,1,drop=FALSE], G=exprs(kidney)[,2,drop=FALSE]))
erge = justvsn(rg)
}

\keyword{}
