\name{meanSdPlot}
\alias{meanSdPlot}
\title{Plot row standard deviations versus row means}
\description{Plot row standard deviations versus row means}
\usage{meanSdPlot(x, ranks=TRUE,
                     xlab = ifelse(ranks, "rank(mean)", "mean"),
                     ylab = "sd",
                     pch  = ".",
                     ...)} 
\arguments{
  \item{x}{An object of class \code{link{exprSet-class}}}
  \item{ranks}{Logical, indicating the x-axis (means) should be plotted
    in original scale (FALSE) or on the rank scale (TRUE). The latter
    distributes the data more evenly along the x-axis and allows a
    better visual of the standard deviation as a function of the mean.}
   \item{xlab}{Character, label for the x-axis.}
   \item{ylab}{Character, label for the y-axis.}
   \item{pch}{Plot symbol.}
   \item{...}{Further arguments that get passed to plot.default.}
}

\details{Standard deviation and mean are calculated row-wise from the
  matrix \code{exprs(x)}. This allows to verify whether there is
  dependence of the standard deviation (or variance) on the mean.
  The blue dots depict the running median estimator (window-width 10\%).
  If there is no variance-mean dependence, then the line formed by the
  blue dots should be approximately horizontal.

  If the \code{preprocessing} slot of the \code{description} slot of
  \code{x} is a \code{\link{list}} and contains an element named
  \code{vsnTrimSelection}, then the coloring of the points reflects
  the trimming that was used in the least trimmed sum of squares (LTS)
  estimation (see \code{\link{vsn}}).
     If the \code{plot} routine is called with the argument
   \code{what="offset"} or \code{what="factors"} or (see examples for
   \code{\link{vsn}}), then the iteration trajectory of the parameters is
   plotted. This helps in assessing whether the parameter estimation in
   \code{\link{vsn}} was sufficiently converged.
}

\value{
  The function is called for its side effect, creating a plot on the
  active graphics device.
}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\seealso{\code{\link{vsn}}}
\examples{
  data(lymphoma)

  if(interactive)
    x11()

  meanSdPlot(lymphoma)

  if(interactive)
    x11()

  meanSdPlot(log(lymphoma))
}

\keyword{hplot}
