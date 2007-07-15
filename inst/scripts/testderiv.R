## Test gradient calculation in vsn C code
##     both for the likelihood with given mu_k, sigma and
##     for the profile likelihood.
##

library("vsn")
options(error=recover)

switch(1L,
       { data("kidney")
         x = exprs(kidney)
         nrpt = 25L ## number of points p0 from which to consider
         strata = cut(1:nrow(x), 3) },
       { data("lymphoma")
         x = exprs(lymphoma)[, 1:8]*1.0
         mode(x) = "double"
         nrpt  = 3L ## number of points p0 from which to consider
         strata = factor(rep(1L, nrow(x)))
       })

nrpar = 2L*ncol(x)*nlevels(strata)
eps   = 1e-4

norm = function(x) sqrt(sum(x*x))

v = new("vsnInput", x=x,
  pstart=array(as.numeric(NA), dim=c(nlevels(strata), ncol(x), 2L)),
  strata=strata, ordered=TRUE)

df = array(NA, dim=c(nrpar, nrpt, 2L))

doit = function(nm, fun) {
  cat("Wait for", nrpt, "points: ")
  for (ip in 1:nrpt) {
    cat(ip, "")
    p0 = runif(nrpar)+1.2
    df[,ip,1] = fun(p0)[-1]

    for(il in 1:nrpar) {
      dp = eps*(il==1:nrpar)
      fn = fun(cbind(p0-dp, p0+dp))[1, ]
      grn = diff(fn)/norm(2*dp)
      stopifnot(is.finite(grn))
      df[il,ip,2]  = grn
      }
  }
  cat("\n\n")

  x11(width=14, height=7)
  par(mfrow=c(2, ncol(x)*nlevels(strata)+1))
  lj = list(offset=1:(ncol(x)*nlevels(strata)),
            factor=ncol(x)*nlevels(strata)+(1:(ncol(x)*nlevels(strata))))
  for(j in seq(along=lj)){
    for(i in seq(along=lj[[j]])) {
      il = lj[[j]][i]
      plot(df[il,,1], df[il,,2], pch=16, xlab="analytic", ylab="numeric",
           main=sprintf("%s %s %d", nm, names(lj)[j], i))
      abline(a=0, b=1, col="blue")
    }
    hist(df[lj[[j]],,1]-df[lj[[j]],,2], col="orange", xlab="difference",
       main=paste(names(lj)[j], "s", sep=""), breaks=30)
    abline(v=0, col="blue")
  }
}

graphics.off()

doit("prof", function(p) logLik(v, p))
doit("ref",  function(p) logLik(v, p, mu = rowMeans(x), sigsq = mean(diff(t(x))^2)))
