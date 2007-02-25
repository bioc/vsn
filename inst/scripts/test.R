###################################################
### chunk number 1: setup
###################################################
library("vsn")
set.seed(0xabcd)
options(error=recover, warn=2)

sim = function(..., lts.quantile=0.9, nrrep=30) {
  callpar = list(...)
  ll      = listLen(callpar)
  stopifnot(ll[1]>=1, all(ll[-1]==1))
  res  = matrix(1, nrow=nrrep, ncol=ll[1])
 
  ## default parameters
  simpar = append(callpar, list(n=4096, d=2, de=0, up=0.5, nrstrata=1, miss=0, log2scale=TRUE)) 
  simpar = simpar[!duplicated(names(simpar))]
  for (i in 1:ll[1]) {
    simpar[[1]] = callpar[[1]][i]
    for (r in 1:nrrep) {
      sim = do.call("sagmbSimulateData", simpar)
      ny  = vsn2(sim$y, strata=factor(sim$strata), lts.quantile=lts.quantile, verbose=!TRUE)
      res[r, i] = sagmbAssess(ny@hx, sim)
      cat(i, r, res[r, i], "\n")
      if(res[r, i]>0.1){
        vi = new("vsnInput", x=sim$y, strata=factor(sim$strata),
             pstart=array(as.numeric(NA), dim=c(simpar$nrstrata, simpar$d, 2)))
        vi = vi[order(vi@strata), ]
        vi@ordered=TRUE
        plotVsnLogLik(vi, p=ny@par, whichp=1:2, expand=1000)
        stop("Bah")
      }
   } ## for r
  } ## for i
  return(res)
} ## sim

myPlot = function(n, res, log="xy", ...) {
  matplot(n, t(res), pch=20, log=log, ylab='r.m.s. error', col="orange", xlab=deparse(substitute(n)), ...)
  lines(n, colMeans(res), col="blue")
}


###################################################
### chunk number 3: fign2
###################################################
n = 1000*2^2

res = sim(n=n, nrstrata=8)
myPlot(n, res, main="n: 8 strata")

