\name{meanSdPlot}
\alias{meanSdPlot}
\title{Plot row standard deviations versus row means}
\description{Plot row standard deviations versus row means}
\usage{
meanSdPlot(x, 
           ranks = TRUE,
           xlab  = ifelse(ranks, "rank(mean)", "mean"),
           ylab  = "sd",
           pch   = ".",
           col, ...)} 
\arguments{
\item{x}{An object of class \code{\link[base:matrix]{matrix}} or
  \code{\link[Biobase:exprSet-class]{exprSet}}}
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

\details{Standard deviation and mean are calculated row-wise from the
  matrix \code{exprs(x)}. The scatterplot of these versus each other
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
  The function is called for its side effect, creating a plot on the
  active graphics device.
}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{\code{\link{vsn}}}
\examples{
  data(kidney)
  exprs(kidney) = log.na(exprs(kidney))

  meanSdPlot(kidney)

  ## ...try this out with non-logged data, the lymphoma data, your data...
}

\keyword{hplot}
