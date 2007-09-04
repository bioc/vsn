\name{justvsn}
\alias{justvsn}
\alias{justvsn-methods}
\alias{justvsn,ExpressionSet-method}
\alias{justvsn,AffyBatch-method}
\alias{justvsn,RGList-method}
\alias{justvsn,NChannelSet-method}

\title{Normalization with vsn}
\description{\code{justvsn} accepts the same kind of arguments as 
  \code{\link{vsn2}}, but returns an appropriate 
  \code{\link[Biobase:class.eSet]{eSet}}-derived object
  (see Section \emph{Value} for details). This is in contrast
  to \code{\link{vsn2}}, which always returns
  an object of class \code{\linkS4class{vsn}}.
}

\usage{
justvsn(x, reference, strata, ...)
}

\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted. Methods exist for
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package) and
    \code{\link[limma:rglist]{RGList}} (from the \code{limma} package).}
  \item{reference, strata, ...}{Arguments that get passed on to \code{\link{vsn2}}.}
}

\details{
  The \code{justvsn} methods are implemented such that they call
  the corresponding \code{\link{vsn2}} method and then assemble the
  object that they are going to return from the normalized data and the original
  metadata in \code{x}. However, \code{justvsn} looses the information on the
  model fitting, i.e. what is stored in the \code{\linkS4class{vsn}}
  object produced by \code{link{vsn2}}.

  If \code{x} is an \code{\link[affy:AffyBatch-class]{AffyBatch}}, probe-wise
  background correction and between-array normalization are done by
  \code{\link{vsn2}} on the perfect match (PM) values only. Probeset
  summaries are calculated with the medianpolish algorithm of
  \code{\link[affy:rma]{rma}}. If you do not want to compute probeset
  summaries, please extract the intensities from the \code{AffyBatch} into a
  matrix and call \code{\link{vsn2}} on that.
  
  If \code{x} is an \code{\link[limma:rglist]{RGList}}
  or an \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, then
  \code{\link{vsn2}} is applied to the matrix that is obtained
  by horizontally concatenating the color channels.
  Optionally, available background estimates can be subtracted before.
}

\value{
  If \code{x} is an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} or an
  \code{\link[affy:AffyBatch-class]{AffyBatch}},
  the return value of this function is an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.

  If \code{x} is an
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}} or an
  \code{\link[limma:rglist]{RGList}}, the return value is an
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}.
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
require("limma")
rg = new("RGList", list(R=exprs(kidney)[,1,drop=FALSE], G=exprs(kidney)[,2,drop=FALSE]))
erge = justvsn(rg)
}

\keyword{}
