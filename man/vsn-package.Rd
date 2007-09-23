\name{vsn-package}
\alias{vsn-package}
%% \alias{vsn}
\docType{package}
\title{vsn}
\description{vsn}
\details{
  The main function of the package is \code{\link{vsn2}}.
  Interesting for its applications are also
  \code{predict} and the wrapper function \code{\link{justvsn}}.
  
  \code{\link{vsn2}} can be applied to objects of class
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[Biobase:class.NChannelSet]{NChannelSet}},
  \code{\link[affy:AffyBatch-class]{AffyBatch}} (from the \code{affy} package) and
  \code{\link[limma:rglist]{RGList}} (from the \code{limma} package),
  \code{matrix} and \code{vector}. It returns an object of class
  \code{\linkS4class{vsn}}, which contains the results of fitting the
  \code{vsn} model to the data.

  The most common use case is that you will want to construct a new
  data object with the vsn-normalized data whose class is the same as
  that of the input data and which preserves the metadata. This can be
  achieved by

  \preformatted{
    fit = vsn2(x, ...)
    nx = predict(fit, newdata=x)
  }

  To simplify this, there exists also a simple wrapper
  \code{\link{justvsn}}.
}

\author{
   Wolfgang Huber
}
\keyword{package}
