\name{scalingFactorTransformation}
\alias{scalingFactorTransformation}

\title{The transformation that is applied to the scaling parameter of
  the vsn model}
\description{The transformation that is applied to the scaling parameter of
  the vsn model}
\usage{
scalingFactorTransformation(b)
}

\arguments{
  \item{b}{Real vector.}
}
\value{
  A real vector of same length as b, with transformation \code{f}
  applied (see
  vignette \emph{Likelihood Calculations for vsn}).
}

\author{Wolfgang Huber}

\examples{
b  = seq(-3, 2, length=20)
fb = scalingFactorTransformation(b)
if(interactive())
  plot(b, fb, type="b", pch=16)
}


