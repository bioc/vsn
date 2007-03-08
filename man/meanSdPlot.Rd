\name{meanSdPlot}
\docType{methods}
\alias{meanSdPlot}
\alias{meanSdPlot-methods}
\alias{meanSdPlot,matrix-method}
\alias{meanSdPlot,eSet-method}
\alias{meanSdPlot,exprSet-method}
\title{Plot row standard deviations
  versus row means}
\description{Methods for objects of classes
\code{eSet}, \code{exprSet} and \code{matrix}
  to plot row standard deviations versus row means}
\usage{
meanSdPlot(x, 
           ranks = TRUE,
           xlab  = ifelse(ranks, "rank(mean)", "mean"),
           ylab  = "sd",
           pch   = ".",
           col, ...)} 
\arguments{
\item{x}{An object of class \code{\link[base:matrix]{matrix}} or
  \code{\linkS4class{ExpressionSet}}.}
\item{ranks}{Logical, indicating whether the x-axis (means) should be plotted
  on the original scale (FALSE) or on the rank scale (TRUE). The latter
  distributes the data more evenly along the x-axis and allows a
  better visual assessment of the standard deviation as a function of 
  the mean.}
\item{xlab}{Character, label for the x-axis.}
\item{ylab}{Character, label for the y-axis.}
\item{pch}{Plot symbol.}
\item{col}{Color of plotted points. See details.}
\item{...}{Further arguments that get passed to plot.default.}
}

\section{Methods}{
\describe{

\item{meanSdPlot (signature \code{matrix}):}{method for
    objects of class \code{matrix}}
\item{meanSdPlot (signature \code{eSet}):}{method for
    objects of class \code{eSet}}
\item{meanSdPlot (signature \code{exprSet}):}{method for
    objects of class \code{exprSet}}
}}

\details{Standard deviation and mean are calculated row-wise from the
  expression matrix \code{x}. The scatterplot of these versus each other
  allows to visually verify whether there is
  a dependence of the standard deviation (or variance) on the mean.
  The red dots depict the running median estimator (window-width 10\%).
  If there is no variance-mean dependence, then the line formed by the
  red dots should be approximately horizontal.

  If the \code{preprocessing} slot of the \code{description} slot of
  \code{x} is a \code{\link{list}} and contains an element named
  \code{vsnTrimSelection}, then the coloring of the points reflects
  the trimming that was used in the least trimmed sum of squares (LTS)
  estimation (see \code{\link{vsn}}). If the condition does not apply,
  and \code{col} is \code{NULL}, the points are drawn in black. If
  \code{col} is not \code{NULL}, its value is used for the coloring of
  the points.
}

\value{
  The methods are called for their side effect, creating a plot on the
  active graphics device.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{vsn}}}
\examples{
  data(kidney)
  log.na = function(x) log(ifelse(x>0, x, NA))

  exprs(kidney) = log.na(exprs(kidney))
  meanSdPlot(kidney)

  ## ...try this out with non-logged data, the lymphoma data, your data...
}

\keyword{hplot}
\keyword{methods}
