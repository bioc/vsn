\name{vsn}
\alias{vsn}
\title{Variance stabilization and calibration for microarray data. }
\description{Robust estimation of variance-stabilizing and calibrating 
  transformations for microarray data. This is the main function of
  this package; see also the vignette vsn.pdf.}
\usage{
vsn(intensities,
    lts.quantile = 0.5,
    niter        = 10,
    verbose      = TRUE,
    pstart       = NULL,
    describe.preprocessing=TRUE)}
\arguments{
\item{intensities}{An object of class \code{\link{exprSet-class}}
  with raw intensity values from a microarray experiment. Alternatively,
  this may be a \code{\link{matrix}}, a \code{\link{data.frame}} with
  all numeric columns or an \code{\link{marrayRaw}} object (from package
  \code{marrayClasses}). The intensities are assumed to be the raw
  scanner data, summarized over the spots by an image analysis program,
  and possibly "background" subtracted.
  The intensities must not be logarithmically or otherwise transformed,
  and not thresholded or "floored". NAs are not accepted.
  See details.}
\item{lts.quantile}{The quantile that is used for the resistant least
  trimmed sum of squares regression. Allowed values are between
  0.5 and 1, corresponding to least median sum of squares regression,
  and to ordinary least sum of squares regression, respectively.}
\item{niter}{The number of iterations used in the least trimmed sum of
  squares regression.}
\item{verbose}{If TRUE, some messages are printed.}
\item{pstart}{Starting values for the model parameters in the iterative
  parameter estimation algorithm. If NULL, the function tries to determine
  reasonable starting values from the distribution of
  \code{intensities}.}
\item{describe.preprocessing}{If TRUE, calibration and transformation
  parameters, plus some other information are stored in the
  \code{preprocessing} slot of the returned object. See details.}
}

\details{The function calibrates for sample-to-sample variations through
  shifting and scaling, and transforms the intensities to a scale where
  the variance is approximately independent of the mean intensity.
  The variance stabilizing transformation is equivalent to the
  natural logarithm in the high-intensity range, and to a
  linear transformation in the low-intensity range. In an intermediate
  range, the \emph{arsinh} function interpolates smoothly between the
  two. The calibration consists of estimating an offset \code{offs[i]}
  and a scale factor \code{fac[i]} for each column \code{i} of the matrix
  \code{intensities}. Thus, the calibration mapping is:
  
  \code{intensities[k,i] <- intensities[k,i] * fac[i] + offs[i]}

  The parameters \code{offs[i]} and \code{fac[i]} are estimated through
  a robust variant of maximum likelihood. The model assumes that for
  the majority of genes the expression levels are not much different
  across the samples, i.e., that only a minority of genes (less than
  a fraction of \code{lts.quantile}) is differentially expressed.

  The format of the matrix of intensities is as follows:
  For the \bold{two-color printed array technology}, each row
  corresponds to one spot, and the columns to the different arrays
  and wave-lengths (usually red and green, but could be any number).
  For example, if there are 10 arrays, the matrix would have 20 columns,
  columns 1...10 containing the green intensities, and 11...20 the
  red ones. In fact, the ordering of the columns does not matter to
  \code{vsn}, but it is your responsibility to keep track of it for
  subsequent analyses.
  For \bold{one-color arrays}, each row corresponds to a probe, and each
  column to an array.

  \bold{Performance:} This function is slow. That is due to the nested
  iteration loops of the numerical optimization of the likelihood function
  and the heuristic that identifies the non-outlying data points in the
  least trimmed squares regression. For large arrays with many tens of
  thousands of probes, you may want to consider random subsetting: that is,
  only use a subset of the e.g. 10-20,000 rows of the data matrix
  \code{intensities} to fit the parameters, then apply the transformation
  to all the data, using \code{\link{vsnh}}. An example for this can be
  seen in the function \code{\link{normalize.AffyBatch.vsn}}, whose code
  you can inspect by typing \code{normalize.AffyBatch.vsn} on the R
  command line.

  \bold{Calibration and transformation parameters:} The parameters
  are stored in the \code{preprocessing} slot of the \code{description}
  slot of the \code{\link{exprSet-class}} object that is returned, in
  the form of a \code{\link{list}} with three elements
  \itemize{
    \item \code{vsnParams}: a length(2*d) numeric vector of parameters 
    \item \code{vsnParamsIter}: an (2*d) x niter numeric matrix that
    contains the parameter trajectory during the
    iterative fit process (see\code{link{vsnPlotPar}}).
    \item \code{vsnTrimSelection}: a length(n) logical vector that for
    each row of the intensities matrix reports whether it was below
    (TRUE) or above (FALSE) the trimming threshold.
    }

  If \code{intensities} has class \code{\link{exprSet-class}}, and its
  \code{description} slot has class \code{\link{MIAME-class}}, then this
  list is appended to any existing entries in the \code{preprocessing}
  slot. Otherwise, the \code{description} object and its
  \code{preprocessing} slot are created.
}

\value{An object of class \code{\link{exprSet-class}}. Differences
  between the columns of the transformed intensities may be interpreted
  as "regularized" or "shrunken" log-ratios. For the calibration and
  transformation parameters, see the \emph{Details} section.
}

\references{Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.

Parameter estimation for the calibration and variance stabilization 
of microarray data, 
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, 
Annemarie Poustka, and Martin Vingron;  
Statistical Applications in Genetics and Molecular Biology (2003)
Vol. 2 No. 1, Article 3.
http://www.bepress.com/sagmb/vol2/iss1/art3.}

\author{Wolfgang Huber \url{http://www.dkfz.de/mga/whuber}}

\seealso{\code{\link{exprSet-class}}, \code{\link{MIAME-class}},
  \code{\link{normalize.AffyBatch.vsn}}}

\examples{
data(kidney)

if(interactive()) {
  x11(width=9, height=4.5)
  par(mfrow=c(1,2))
}

vsnkid = vsn(kidney)   ## transform and calibrate

plot(exprs(kidney), pch=".", main="log-log", log="xy")
plot(exprs(vsnkid), pch=".", main="h-h")

if (interactive()) {
  x11(width=9, height=4)
  par(mfrow=c(1,3))
}

meanSdPlot(vsnkid)
vsnPlotPar(vsnkid, "factors")
vsnPlotPar(vsnkid, "offsets")

## this should always hold true
params = preproc(description(vsnkid))$vsnParams
stopifnot(all(vsnh(exprs(kidney), params) == exprs(vsnkid))) 
}

\keyword{robust}
