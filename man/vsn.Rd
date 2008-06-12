\name{vsn}
\alias{vsn}
\title{Variance stabilization and calibration for microarray data. }
\description{Robust estimation of variance-stabilizing and calibrating 
  transformations for microarray data.
  This function has been superseded by \code{\link{vsn2}}.
  The function \code{\link{vsn}} remains in the package for backward compatibility,
  but for new projects, please use \code{\link{vsn2}}.}
\usage{
vsn(intensities,
    lts.quantile = 0.5,
    verbose      = interactive(),
    niter        = 10,
    cvg.check    = NULL,
    describe.preprocessing = TRUE,
    subsample,
    pstart,
    strata)
}
\arguments{
\item{intensities}{An object that contains intensity values from
  a microarray experiment. See
  \code{\link[vsn:getIntensityMatrix]{getIntensityMatrix}} for details.
  The intensities are assumed to be the raw
  scanner data, summarized over the spots by an image analysis program,
  and possibly "background subtracted".
  The intensities must not be logarithmically or otherwise transformed,
  and not thresholded or "floored". NAs are not accepted.
  See details.}
\item{lts.quantile}{Numeric. The quantile that is used for the resistant
  least trimmed sum of squares regression. Allowed values are between
  0.5 and 1. A value of 1 corresponds to ordinary least sum of squares
  regression.}
\item{verbose}{Logical. If TRUE, some messages are printed.}
\item{niter}{Integer. The number of iterations to be used in the least
  trimmed sum of squares regression.}
\item{cvg.check}{List. If non-NULL, this allows finer control of the
  iterative least trimmed sum of squares regression. See details.}
\item{pstart}{Array. If not missing, user can specify start values
  for the iterative parameter estimation algorithm. See 
  \code{\link{vsnh}} for details.}
\item{describe.preprocessing}{Logical. If TRUE, calibration and
  transformation parameters, plus some other information are stored in
  the \code{preprocessing} slot of the returned object. See details.}
\item{subsample}{Integer. If specified, the model parameters are
  estimated from a subsample of the data only, the transformation is
  then applied to all data. This can be useful for performance reasons.}
\item{strata}{Integer vector. Its length must be the same as nrow(intensities).
  This parameter allows for the calibration and error model parameters
  to be stratified within each array, e.g to take into account probe 
  sequence properties, print-tip or plate effects.  
  If \code{strata} is not specified, one pair of parameters is fitted
  for every sample (i.e. for every column of \code{intensities}). If
  \code{strata} is specified, a pair of parameters is fitted for every 
  stratum within every sample. The strata are coded for by the different
  integer values. The integer vector \code{strata} can be obtained
  from a factor \code{fac} through \code{as.integer(fac)}, from
  a character vector \code{str} through \code{as.integer(factor(fac))}.}
}
\details{\bold{Overview:} 
  The function calibrates for sample-to-sample variations through
  shifting and scaling, and transforms the intensities to a scale where
  the variance is approximately independent of the mean intensity.
  The variance stabilizing transformation is equivalent to the
  natural logarithm in the high-intensity range, and to a
  linear transformation in the low-intensity range. In an intermediate
  range, the \emph{arsinh} function interpolates smoothly between the
  two. For details on the transformation, please see the help page for
  \code{\link{vsnh}}. The parameters are estimated through
  a robust variant of maximum likelihood. This assumes that for
  the majority of genes the expression levels are not much different
  across the samples, i.e., that only a minority of genes (less than
  a fraction \code{1-lts.quantile}) is differentially expressed.

  Even if most genes on an array are differentially expressed, it may still
  be possible to use the estimator: if a set of non-differentially expressed
  genes is known, e.g. because they are external controls or reliable
  'house-keeping genes', the transformation parameters can be fitted with
  \code{vsn} from the data of these genes, then the transformation can be
  applied to all data with \code{\link{vsnh}}.

  \bold{Format:} The format of the matrix of intensities is as follows:
  for the \bold{two-color printed array technology}, each row
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

  \bold{Iteration control:} 
  By default, if \code{cvg.check} is \code{NULL}, the function will run
  the fixed number \code{niter} of iterations in the least trimmed sum
  of squares regression. More fine-grained control can be obtained by
  passing a list with elements \code{eps} and \code{n}. If the maximum
  change between transformed data values is smaller than \code{eps} for
  \code{n} subsequent iterations, then the iteration terminates.

  \bold{Estimated transformation parameters:} 
  If \code{describe.preprocessing} is \code{TRUE}, the transformation
  parameters are returned in the \code{preprocessing} slot of the
  \code{experimentData} slot of the resulting 
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} object, in the form 
  of a \code{\link{list}} with three elements
  \itemize{
  \item \code{vsnParams}: the parameter array (see \code{\link{vsnh}} 
    for details) 
  \item \code{vsnParamsIter}: an array with dimensions 
    \code{c(dim(vsnParams, niter))} that contains the parameter 
    trajectory during the iterative fit process (see also 
    \code{\link{vsnPlotPar}}).
  \item \code{vsnTrimSelection}: a logical vector that for
    each row of the intensities matrix reports whether it was below
    (TRUE) or above (FALSE) the trimming threshold.
  }
  If \code{intensities} has class \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}},
  and its \code{experimentData} slot has class
  \code{\link[Biobase:MIAME-class]{MIAME}}, then this list is appended to any
  existing entries in the \code{preprocessing} slot. Otherwise, the
  \code{experimentData} object and its \code{preprocessing} slot are created.  
}

\value{An object of class \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
 Differences between the columns of the transformed intensities are 
 "generalized log-ratios", which are shrinkage estimators of the natural
 logarithm of the fold change. For the transformation parameters,
 please see the Details.
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

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\seealso{\code{\link{vsnh}}, \code{\link{vsnPlotPar}}, 
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}, 
  \code{\link[Biobase:MIAME-class]{MIAME-class}},
  \code{\link{normalize.AffyBatch.vsn}}}

\examples{
data(kidney)
log.na = function(x) log(ifelse(x>0, x, NA))

if(interactive()) {
  x11(width=9, height=4.5)
  par(mfrow=c(1,2))
}
plot(log.na(exprs(kidney)), pch=".", main="log-log")

vsnkid = vsn(kidney)   ## transform and calibrate
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
