\name{meanSdPlot}
\docType{methods}
\alias{meanSdPlot}
\alias{meanSdPlot-methods}
\alias{meanSdPlot,matrix-method}
\alias{meanSdPlot,ExpressionSet-method}
\alias{meanSdPlot,vsn-method}
\alias{meanSdPlot,MAList-method}
\title{Plot row standard deviations versus row means}

\description{Methods for objects of classes
  \code{\link{matrix}}, 
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\linkS4class{vsn}} and
  \code{\link[limma:malist]{MAList}}
  to plot row standard deviations versus row means.}
\usage{
meanSdPlot(x, 
           ranks = TRUE,
           xlab  = ifelse(ranks, "rank(mean)", "mean"),
           ylab  = "sd",
           pch,
           plot  = TRUE,
	   bins  = 50,
           \dots)} 
\arguments{
  \item{x}{An object of class
    \code{\link{matrix}}, 
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, 
    \code{\linkS4class{vsn}} or
    \code{\link[limma:malist]{MAList}}.}
\item{ranks}{Logical, indicating whether the x-axis (means) should be plotted
  on the original scale (\code{FALSE}) or on the rank scale (\code{TRUE}). The latter
  distributes the data more evenly along the x-axis and allows a
  better visual assessment of the standard deviation as a function of 
  the mean.}
\item{xlab}{Character, label for the x-axis.}
\item{ylab}{Character, label for the y-axis.}
\item{pch}{Ignored - exists for backward compatibility.}
\item{plot}{Logical.  If \code{TRUE} (default), a plot is produced.
  Calling the function with \code{plot=FALSE} can be useful if only
  its return value is of interest.}
\item{bins}{Gets passed on to \code{\link[ggplot2:stat_binhex]{stat_binhex}}.}
\item{\dots}{Further arguments that get passed on to \code{\link[ggplot2:stat_binhex]{stat_binhex}}.}
}
\details{Standard deviation and mean are calculated row-wise from the
  expression matrix (in) \code{x}. The scatterplot of these versus each other
  allows you to visually verify whether there is a dependence of the standard
  deviation (or variance) on the mean.
  The red line depicts the running median estimator (window-width 10\%).
  If there is no variance-mean dependence, then the line should be approximately horizontal.
}
\value{
  A named list with five components: its elements \code{px} and
  \code{py} are the x- and y-coordinates of the individual data points
  in the plot; its first and second element are the x-coordinates and values of
  the running median estimator (the red line in the plot).
  Its element \code{gg} is the plot object.
  Depending on the value of \code{plot}, the method can also have a side effect,
  which is to print \code{gg} on the active graphics device.
}

\author{Wolfgang Huber}
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
