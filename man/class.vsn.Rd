\name{vsn}
\docType{class}
\alias{class:vsn}
\alias{vsn-class}

\alias{[,vsn-method} 
\alias{dim,vsn-method} 
\alias{nrow,vsn-method} 
\alias{ncol,vsn-method} 
\alias{show,vsn-method} 
\alias{exprs,vsn-method} 
\alias{coef,vsn-method} 
\alias{coefficients,vsn-method} 

\title{
  Class to contain result of a vsn fit
}

\description{
  Class to contain result of a vsn fit
}

\section{Creating Objects}{
  \code{new("vsn")}
  \code{vsn2(x)} with \code{x} being an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}

\section{Slots}{
 \describe{
  \item{\code{coefficients}:}{A 3D array of size (number of strata) x (number of
    columns of the data matrix) x 2 containing fitted normalization parameters
    (see vignette).}
  \item{\code{strata}:}{A factor of length 0 or n. If its length is
    n, then its levels correspond to different normalization strata
    (see vignette).}
  \item{\code{mu}:}{A numeric vector of length n with the fitted
    parameters $\hat{\mu}_k$, $k=1,...,n$.}
  \item{\code{sigsq}:}{A numeric scalar, $\hat{\sigma}^2$.}
  \item{\code{hx}:}{A numeric matrix with 0 or n rows. If the number of
    rows is n, then \code{hx} contains the transformed data matrix.}
  \item{\code{lbfgsb}:}{An integer scalar containing the return code
    from the L-BFGS-B optimizer.}
  }
}


\section{Methods}{
  \describe{
    \item{\code{[}}{Subset}
    \item{\code{dim}}{Get dimensions of data matrix.}
    \item{\code{nrow}}{Get number of rows of data matrix.}
    \item{\code{ncol}}{Get number of columns of data matrix.}
    \item{\code{show}}{Print a summary of the object}
    \item{\code{exprs}}{Accessor to slot \code{hx}.}
    \item{\code{coef}, \code{coefficients}}{Accessors to slot \code{coefficients}.}
}}


\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\seealso{
  \code{\link{vsn2}}
}

\examples{
  data("kidney")
  v = vsn2(kidney)
  show(v)
  dim(v)
  v[1:10, ]
}

\keyword{classes}
