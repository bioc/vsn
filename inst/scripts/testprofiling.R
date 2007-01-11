## Test concordance of maxima of profile likelihood and normal likelihood

library("vsn")
data("kidney")

norm       = function(x) sqrt(sum(x*x))
almostequal= function(x,y) stopifnot(max(abs(x-y)/abs(x+y))<1e-10)

x = exprs(kidney)
nrstr = 1
nrpar = 2*ncol(x)*nrstr
istrat = as.integer(seq(0, 1, length=ncol(x)*nrstr+1)*ncol(x)*nrow(x))

if(!exists("fit"))
  fit = vsn2(x, lts.quantile=1)

v = new("vsnInput", x=exprs(kidney),
  pstart=array(as.numeric(NA), dim=c(1, ncol(kidney), 2)))

nplot = 41
par(mfrow=c(2, nrpar/2))

for(i in seq_len(nrpar))  {
  pars = matrix(fit@par, nrow=nrpar, ncol=nplot)
  if(i<=ncol(x)*nrstr) {
    pars[i,] = fit@par[i] +  (seq(-1, 1, length=nplot) * 0.1 *
          diff(quantile(fit@data[,i], probs=c(0.01, 0.99))))
    xlab = substitute(a[k], list(k=i))
    logarg  = ""
  } else {
    pars[i,] = fit@par[i] * exp(seq(-1, 1, length=nplot))    
    xlab = substitute(b[k], list(k=i-nrpar/2))
    logarg  = "x"
  }

  ## ll = vsnLikelihood(v, pars)
  ll = vsnLikelihood(v, pars, refh=fit@refh, refsigma=fit@refsigma)

  plot(pars[i,], ll[1, ], type="l",
       xlab = xlab, log=logarg,
       ylab = expression(-log(L)))
  abline(v=fit@par[i], col="red")
}



