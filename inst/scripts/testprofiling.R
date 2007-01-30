##
## Test concordance of maxima of profile likelihood and normal likelihood
##

library("vsn")

n = 10000
d = 2
nrstr = 1
nrpar = 2*d*nrstr
istrat = as.integer(seq(0, 1, length=d*nrstr+1)*d*n)

if(!exists("fit")) {
  dat = sagmbSimulateData(n=n, d=d, de=0, nrstrata=nrstr, miss=0, log2scale=TRUE)
  fit = vsn2(dat$y, lts.quantile=1)
  v = new("vsnInput", x=dat$y, pstart=array(as.numeric(NA), dim=c(nrstr, d, 2)))
}

refh = fit@refh
refsigma = fit@refsigma

nplot = 41
par(mfcol=c(2, nrpar))

for(i in seq_len(nrpar))  {
  pars = matrix(fit@par, nrow=nrpar, ncol=nplot)
  if(i<=d*nrstr) {
    pars[i,] = fit@par[i] +  (seq(-1, 1, length=nplot) * 0.1 *
          diff(quantile(fit@hx[,i], probs=c(0.01, 0.99))))
    xlab = substitute(a[k], list(k=i))
    logarg  = ""
  } else {
    pars[i,] = fit@par[i] * exp(seq(-1, 1, length=nplot))    
    xlab = substitute(b[k], list(k=i-nrpar/2))
    logarg  = "x"
  }
  
  for(what in 1:2){
    ll = switch(what,
      logLik(v, pars),                               ## without reference
      logLik(v, pars, refh=refh, refsigma=refsigma)) ## with reference

    plot(pars[i,], ll[1, ], type="l",
         xlab = xlab, log=logarg, main=round(ll[1, (ncol(ll)+1)/2], 1),
         ylab = expression(-log(L)))
    abline(v=fit@par[i], col="red")
  } ## for what
}



