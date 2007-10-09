\name{vsn2}
\alias{vsn2}
\alias{vsn2-methods}
\alias{vsn2,matrix-method}
\alias{vsn2,numeric-method}
\alias{vsn2,ExpressionSet-method}
\alias{vsn2,AffyBatch-method}
\alias{vsn2,RGList-method}
\alias{vsn2,NChannelSet-method}
\alias{coerce,RGList,NChannelSet-method}
\alias{vsnMatrix}

\title{Fit the vsn model}
\description{\code{vsn2} fits the vsn model to the data 
  in \code{x} and returns a \code{\linkS4class{vsn}} object with
  the fit parameters and the transformed data matrix.
  The data are, typically, feature intensity readings from a
  microarray.
  To obtain an object of the same class as \code{x}, containing
  the normalised data and the same metdata as \code{x}, use
  \preformatted{
    fit = vsn2(x, ...)
    nx = predict(fit, newdata=x)
  }
  or the wrapper \code{\link{justvsn}}.
  Please see the vignette \emph{Introduction to vsn} for a description
  on how to use \code{vsn2} for different use cases.
}

\usage{
vsnMatrix(x,
          reference,
          strata,
          lts.quantile = 0.9,
          subsample    = 0L,
          verbose      = interactive(),
          returnData   = TRUE,
          pstart,
          minDataPointsPerStratum = 42L,
          optimpar     = list(),
          defaultpar   = list(factr=5e7, pgtol=2e-4, maxit=60000L,
                              trace=0L, cvg.niter=7L, cvg.eps=0))

\S4method{vsn2}{ExpressionSet}(x, reference, strata, ...)

\S4method{vsn2}{AffyBatch}(x, reference, strata, ...)

\S4method{vsn2}{matrix}(x, reference, strata, ...)

\S4method{vsn2}{NChannelSet}(x, reference, strata, backgroundsubtract=FALSE,
       foreground=c("R","G"), background=c("Rb", "Gb"), ...)

\S4method{vsn2}{RGList}(x, reference, strata, backgroundsubtract=FALSE, ...)
}
\arguments{
  \item{x}{An object containing the data to which the model is to be
    fitted.}
  \item{reference}{Optional, a \code{\linkS4class{vsn}} object from
    a previous fit. If this argument is specified, the data in \code{x}
    are normalized "towards" an existing set of reference arrays whose
    parameters are stored in the object \code{reference}. If this
    argument is not specified, then the data in \code{x} are normalized
    "among themselves". See Details for a more precise explanation.}
  \item{strata}{Optional, a \code{factor} or \code{integer}
    whose length is \code{nrow(x)}. Can
    be used for stratified normalization (i.e. separate offsets \code{a} and
    factors \code{b} for each level of \code{strata}). If missing, all
    rows of \code{x} are assumed to come from one stratum.
    If \code{strata} is an integer, its values must cover the range
    1...\emph{n}, where \emph{n} is the number of strata.}
  \item{lts.quantile}{Numeric of length 1. The quantile that is used for
    the resistant least trimmed sum of squares regression. Allowed
    values are between 0.5 and 1. A value of 1 corresponds to ordinary
    least sum of squares regression.}
  \item{subsample}{Integer of length 1. If specified, the model parameters are
    estimated from a subsample of the data of size \code{subsample}
    only, yet the fitted transformation is
    then applied to all data. For large datasets, this can substantially
    reduce the CPU time and memory consumption at a negligible loss of precision.}
  \item{backgroundsubtract}{Logical of length 1: should
    local background estimates be subtracted before fitting vsn?}
  \item{foreground, background}{Aligned character vectors of the same length,
    naming the channels of \code{x} that should be used
    as foreground and background values.}
  \item{verbose}{Logical. If TRUE, some messages are printed.}
  \item{returnData}{Logical. If TRUE, the transformed data are returned
    in a slot of the resulting \code{\linkS4class{vsn}} object.
    Setting this option to \code{FALSE} allows saving memory
    if the data are not needed.}
  \item{pstart}{Optional, a three-dimensional numeric array that
    specifies start values for the iterative parameter
    estimation algorithm.
    If not specified, the function tries to guess useful start values.
    The first dimension corresponds to the levels of \code{strata},
    the second dimension to the columns of \code{x} and the third dimension
    must be 2, corresponding to offsets and factors.}
  \item{minDataPointsPerStratum}{ The minimum number of data points per
    stratum. }
  \item{optimpar}{Optional, a list with parameters for the likelihood
    optimisation algorithm. Default parameters are taken from
    \code{defaultpar}. See details.}
  \item{defaultpar}{The default parameters  for the likelihood
    optimisation algorithm. Values in \code{optimpar} take precedence
    over those in \code{defaultpar}. The purpose of this argument is to
    expose the default values in this manual page - it is not
    intended to be changed, please use \code{optimpar} for that.}
  \item{...}{Arguments that get passed on to \code{vsnMatrix}.}
}

\value{An object of class \code{\linkS4class{vsn}}.}

\section{Note on overall scale and location of the glog
  transformation}{
  The data are returned on a glog scale to base 2. More precisely,
  the transformed data are subject to the transformation
  glog2(f(b)*x+a) + c, where
  glog2(u) = log2(u+sqrt(u*u+1)) = asinh(u)/log(2) is called the
  generalised logarithm, a and b are the fitted model parameters
  (see references), f is a parameter transformation [4],
  and the overall constant offset c is computed from b such that for
  large x the transformation approximately corresponds to the
  log2 function. The offset c is inconsequential for all differential
  expression calculations, but many users like to see the data in a
  range that they are familiar with.
}

\section{Specific behaviour of the different methods}{
  \code{vsn2} methods exist for
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
  \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package),
  \code{\link[limma:rglist]{RGList}} (from the \code{limma} package),
  \code{matrix} and \code{numeric}.
  If \code{x} is an \code{\link[Biobase:class.NChannelSet]{NChannelSet}}, then
  \code{vsn2} is applied to the matrix that is obtained
  by horizontally concatenating the color channels.
  Optionally, available background estimates can be subtracted before.
  If \code{x} is an \code{\link[limma:rglist]{RGList}}, it is
  converted into an \code{NChannelSet}
  using a copy of Martin Morgan's code for \code{RGList} to
  \code{NChannelSet} coercion, then the \code{NChannelSet} method is called.
}

\section{Standalone versus reference normalisation}{
  If the \code{reference} argument is \emph{not} specified, then the model
  parameters $\mu_k$ and $\sigma$ are fit from the data in \code{x}.
  This is the mode of operation described in [1]
  and that was the only option in versions 1.X of this package.
  If \code{reference} is specified, the model parameters
  $\mu_k$ and $\sigma$ are taken from it.
  This allows for 'incremental' normalization [4].}

\section{Convergence of the iterative likelihood optimisation}{
  \code{L-BFGS-B} uses three termination criteria:
  \enumerate{
    \item \code{(f_k - f_{k+1}) / max(|f_k|, |f_{k+1}|, 1) <= factr * epsmch}
      where \code{epsmch} is the machine precision.
    \item \code{|gradient| < pgtol}
    \item \code{iterations > maxit}
  }
  These are set by the elements \code{factr}, \code{pgtol} and
  \code{maxit} of \code{optimpar}. The remaining elements are
  \describe{
    \item{\code{trace}}{An integer between 0 and 6, indicating the
      verbosity level of \code{L-BFGS-B}, higher values
      create more output.}
    \item{\code{cvg.niter}}{The number of iterations to be used in the least
      trimmed sum of squares regression.}
    \item{\code{cvg.eps}}{Numeric. A convergence threshold for the least
      trimmed sum of squares regression.}
  }
}

\seealso{\code{\link{justvsn}}, \code{\link[vsn:vsn2trsf]{predict}}}

\references{
[1] Variance stabilization applied to microarray data
calibration and to the quantification of differential expression,
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104.

[2] Parameter estimation for the calibration and variance stabilization 
of microarray data, 
Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, 
Annemarie Poustka, and Martin Vingron;  
Statistical Applications in Genetics and Molecular Biology (2003)
Vol. 2 No. 1, Article 3.
http://www.bepress.com/sagmb/vol2/iss1/art3.

[3] L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
Optimization, C. Zhu, R.H. Byrd, P. Lu and J. Nocedal,
Technical Report, Northwestern University (1996).

[4] Package vignette: Likelihood Calculations for vsn
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
data("kidney")

fit = vsn2(kidney)                   ## fit
nkid = predict(fit, newdata=kidney)  ## apply fit

plot(exprs(nkid), pch=".")
abline(a=0, b=1, col="red")
}

\keyword{}
