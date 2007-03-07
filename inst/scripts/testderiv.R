## Test gradient calculation in vsn C code
##     both for the likelihood with given mu_k, sigma and
##     for the profile likelihood.
##

library("vsn")
options(error=recover)
data("kidney")
x = exprs(kidney)

nrpt  = 25 ## number of points p0 from which to consider
nrstr = 2
nrpar = 2*ncol(x)*nrstr
eps   = 1e-4

norm = function(x) sqrt(sum(x*x))

v = new("vsnInput", x=exprs(kidney),
  pstart=array(as.numeric(NA), dim=c(nrstr, ncol(kidney), 2)),
  strata=factor(seq(0, 1, length=nrow(kidney))<0.5), ordered=TRUE)

df = array(NA, dim=c(nrpar, nrpt, 2))

doit = function(fun) {
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
      ##browser()
    }
  }
  cat("\n\n")

  x11(width=14, height=7)
  par(mfrow=c(2, ncol(x)*nrstr+1))
  k = c(1, 2)
  lj = list(offset=1:(ncol(x)*nrstr), factor=ncol(x)*nrstr+(1:(ncol(x)*nrstr)))
  for(j in seq(along=lj)){
    for(il in lj[[j]]) {
      plot(df[il,,k[1]], df[il,,k[2]], pch=16, xlab=paste(k[1]), ylab=paste(k[2]),
           main=sprintf("%s %d", names(lj)[j], il))
      abline(a=0, b=1, col="blue")
    }
    hist(df[lj[[j]],,k[1]]-df[lj[[j]],,k[2]], col="orange", xlab="difference",
       main=paste(names(lj)[j], "s", sep=""), breaks=30)
    abline(v=0, col="blue")
  }
}

graphics.off()

doit(function(p) logLik(v, p))
doit(function(p) logLik(v, p, refh = rowMeans(x), refsigma = mean(diff(t(x)))))
