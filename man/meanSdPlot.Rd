\name{meanSdPlot}
\docType{methods}
\alias{meanSdPlot}
\alias{meanSdPlot-methods}
\alias{meanSdPlot,matrix-method}
\alias{meanSdPlot,ExpressionSet-method}
\alias{meanSdPlot,vsn-method}
\alias{meanSdPlot,exprSet-method}
\title{Plot row standard deviations versus row means}

\description{Methods for objects of classes
  \code{\link{matrix}}, 
  \code{\linkS4class{ExpressionSet}} and
  \code{\linkS4class{exprSet}} (deprecated)
  to plot row standard deviations versus row means}
\usage{
meanSdPlot(x, 
           ranks = TRUE,
           xlab  = ifelse(ranks, "rank(mean)", "mean"),
           ylab  = "sd",
           pch   = ".", ...)} 
\arguments{
  \item{x}{An object of class
    \code{\link{matrix}}, 
    \code{\linkS4class{ExpressionSet}},
    \code{\linkS4class{vsn}} or
    \code{\linkS4class{exprSet}} (deprecated).}
\item{ranks}{Logical, indicating whether the x-axis (means) should be plotted
  on the original scale (\code{FALSE}) or on the rank scale (\code{TRUE}). The latter
  distributes the data more evenly along the x-axis and allows a
  better visual assessment of the standard deviation as a function of 
  the mean.}
\item{xlab}{Character, label for the x-axis.}
\item{ylab}{Character, label for the y-axis.}
\item{pch}{Plot symbol.}
\item{...}{Further arguments that get passed to plot.default.}
}

\details{Standard deviation and mean are calculated row-wise from the
  expression matrix (in) \code{x}. The scatterplot of these versus each other
  allows to visually verify whether there is
  a dependence of the standard deviation (or variance) on the mean.
  The red dots depict the running median estimator (window-width 10\%).
  If there is no variance-mean dependence, then the line formed by the
  red dots should be approximately horizontal.
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
