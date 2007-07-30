\name{justvsn}
\alias{justvsn}
\alias{justvsn-methods}
\alias{justvsn,ExpressionSet-method}
\alias{justvsn,AffyBatch-method}
\alias{justvsn,RGList-method}
\alias{justvsn,NChannelSet-method}

\title{Normalization with vsn}
\description{\code{justvsn} fits the vsn model to data
  and returns an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
  or an
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}.
  It is a simple wrapper for
  \code{\link{vsn2}} and \code{\link[vsn:vsn2trsf]{predict}}.
}

\usage{
justvsn(x, reference, strata, ...)

\S4method{justvsn}{ExpressionSet}(x, reference, strata, ...)
\S4method{justvsn}{AffyBatch}(x, reference, strata, ...)
\S4method{justvsn}{RGList}(x, reference, strata, backgroundsubtract=FALSE, ...)
\S4method{justvsn}{NChannelSet}(x, reference, strata, backgroundsubtract=FALSE,
       foreground=c("R","G"), background=c("Rb", "Gb"), ...)

}

\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted. Methods exist for
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
    \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
    \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package) and
    \code{\link[limma:rglist]{RGList}} (from the \code{limma} package).}
  \item{reference, strata, ...}{As in \code{\link{vsn2}}.}
  \item{backgroundsubtract}{Logical of length 1: should
    local background estimates be subtracted before fitting vsn?
    For \code{\link[limma:rglist]{RGList}}, the local background
    estimates are expected in the \code{Rb} and \code{Gb} slots of
    \code{x} and the foreground values in the \code{R} and \code{G}
    slots.}
  \item{foreground, background}{Character vectors of the same length,
    naming the channels of \code{x} that should be used
    as foreground and background values.
  }
}

\details{
  See \code{\link{vsn2}}.

  If \code{x} is an \code{\link[affy:AffyBatch-class]{AffyBatch}}, probe-wise
  background correction and between-array normalization are done by
  \code{\link{vsn2}} on the perfect match (PM) values only. Probeset
  summaries are calculated with the medianpolish algorithm of
  \code{\link[affy:rma]{rma}}.
  
  If \code{x} is an \code{\link[limma:rglist]{RGList}}
  or an \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, then
  \code{\link{vsn2}} is applied to the matrix that is obtained
  by horizontally concatenating the color channels.
  Optionally, available background estimates can be subtracted before. 
}

\value{
  If \code{x} is an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[affy:AffyBatch-class]{AffyBatch}} or
  \code{\link[limma:rglist]{RGList}}, the value of this function is an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.

  If \code{x} is an
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, the value is
  also an
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
