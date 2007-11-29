\name{vsnInput}
\docType{class}
\alias{class:vsnInput}
\alias{vsnInput}
\alias{vsnInput-class}

\alias{[,vsnInput-method} 
\alias{dim,vsnInput-method} 
\alias{nrow,vsnInput-method} 
\alias{ncol,vsnInput-method} 
\alias{show,vsnInput-method} 

\title{
  Class to contain input data and parameters for vsn functions
}

\description{
  Class to contain input data and parameters for vsn functions
}

\section{Creating Objects}{
  \code{new("vsnInput")}
}

\section{Slots}{
 \describe{
  \item{\code{x}:}{A numeric matrix with the input data.}
  \item{\code{reference}:}{An object of \code{\linkS4class{vsn}},
    typically this would have been
    obtained from a previous fit to a set of reference arrays (data).}
  \item{\code{strata}:}{A factor of length 0 or n. If its length is
    n, then its levels correspond to different normalization strata
    (see \code{\link{vsn2}}).}
  \item{\code{ordered}:}{Logical scalar; are the rows reordered so that
    the strata are contiguous.}
  \item{\code{lts.quantile}:}{Numeric scalar, see\code{\link{vsn2}}.}
  \item{\code{subsample}:}{Integer scalar, see\code{\link{vsn2}}.}
  \item{\code{verbose}:}{Logical scalar, see\code{\link{vsn2}}.}
  \item{\code{pstart}:}{A 3D array of size (number of strata) x (number of
    columns of the data matrix) x 2 containing start parameters.}
  \item{\code{cvg.niter}:}{Integer scalar, see\code{\link{vsn2}}.}
  \item{\code{cvg.eps}:}{Numeric scalar, see\code{\link{vsn2}}.}
  }
}


\section{Methods}{
  \describe{
    \item{\code{[}}{Subset}
    \item{\code{dim}}{Get dimensions of data matrix.}
    \item{\code{nrow}}{Get number of rows of data matrix.}
    \item{\code{ncol}}{Get number of columns of data matrix.}
    \item{\code{show}}{Print a summary of the object}
}}



\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\seealso{
  \code{\link{vsn2}}
}

\keyword{classes}
