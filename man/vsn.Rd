\name{vsn}
\alias{vsn}
\title{Variance stabilization and calibration for microarray data. }
\description{Robust estimation of variance-stabilizing and calibrating 
  transformations for microarray data. This is the main function of
  this package; see also the vignette vsn.pdf.}
\usage{vsn(intensities, lts.quantile=0.5, niter=10, verbose=TRUE, pstart=NULL)}
\arguments{
  \item{intensities}{Intensity values from a microarray experiment. 
    This may be a \code{\link{matrix}}, a \code{\link{data.frame}}
    with all numeric columns, an \code{\link{exprSet}} object
    (from package \code{Biobase}) or an \code{\link{marrayRaw}} object 
    (from package \code{marrayClasses}).
    For the two-color printed array technology,
    each row corresponds to one spot, and the columns to the different
    wave-lengths (usually red and green, but could be any number).
    Intensities may be "background" subtracted.
    For one-color arrays, each row corresponds to a probe, and each
    column to an array. The statistical model 
    that is assumed by this function expects the raw intensity data, 
    not logarithmically or otherwise transformed, and not thresholded
    or "floored". NAs are not accepted.}
  \item{lts.quantile}{The quantile that is used for the resistant least
    trimmed sum of squares regression. Allowed values are between
    0.5 and 1, corresponding to least median sum of squares regression,
    and to ordinary least sum of squares regression, respectively.}
  \item{niter}{The number of iterations used in the least trimmed sum of
    squares regression.}
  \item{verbose}{If TRUE, some mesages are printed.}
  \item{pstart}{Starting values for the model parameters in the iterative
    parameter estimation algorithm. If NULL, the function tries to determine
    reasonable starting values from the distribution of \code{intensities}.}
}

\details{The function performs (1) calibration for sample-to-sample
variations through shifting and scaling and (2) a variance stabilizing
transformation. The variance stabilizing transformation is equivalent to
a log-transformation (base e) in the high-intensity range, and to a
linear transformation in the low-intensity. In an intermediate range,
the "arsinh" function interpolates smoothly between the two. The
calibration consists of estimating an offset \code{a[i]} and a scale
factor \code{b[i]} for each column \code{i} of the matrix
\code{intensities}. Thus, the calibration mapping is:

      \code{intensities[k,i] <- intensities[k,i]*b[i] + a[i]}

The \code{a[i]} and \code{b[i]} are estimated through a robust weighted
regression method. The regression assumes that for the majority of genes
the expression levels are not much different across the samples, i.e.,
that only a minority of genes is differentially expressed.

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
}

\value{  
  An object of class \link{vsn.result}. The slot \code{h} contains the
  calibrated and transformed intensities (see examples). Differences
  between the transformed may be interpreted in much the same way as
  log-ratios. In the high intensity range, i.e., for intensities much
  above the background noise level, they are in fact identical to the
  convential log-ratios; In the low intensity range, i.e. for
  intensities close to the background noise level, they can be
  intepreted as "regularized" or "shrunken" log-ratios,
}

\references{Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.

Parameter estimation for the calibration and variance stabilization
of microarray data. Wolfgang Huber, Anja von Heydebreck, Holger
Sueltmann, Annemarie Poustka, Martin Vingron. To appear in Statistical
Applications in Genetics and Molecular Biology.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}

\seealso{\code{\link{vsn.result}}, \code{\link{normalize.AffyBatch.vsn}}}
\examples{
data(lymphoma)

if (interactive())
  x11(width=9, height=4.5)

y = exprs(lymphoma)[, 3:4]
colnames(y) = pData(lymphoma)$sample[3:4]

par(mfrow=c(1,2))
plot(y, log="xy", pch=".", main="log-log")

ny = vsn(y)

plot(ny@h, pch=".", main="h-h")
show(ny)

if (interactive())
  x11()

plot(ny)
plot(ny, what="offsets")
plot(ny, what="factors")

## this should always hold true
stopifnot(all(vsnh(y, params(ny)) == ny@h)) 
}

\keyword{robust}
