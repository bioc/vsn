plotLogLik = function(x, strata, p0, wh, n=31) {

  stopifnot(is.matrix(wh), nrow(wh)==2, ncol(wh)==3)
  
  v = new("vsnInput",
    x = x,
    strata = strata,
    pstart = p0,
    lts.quantile = 1,
    cvg.niter  = as.integer(10),
    cvg.eps   = 0,
    subsample = as.integer(0),
    verbose   = FALSE,
    ordered   = FALSE)

  if((length(v@strata)>0)&&(nlevels(v@strata)>1)) {
    ord = order(v@strata)
    v   = v[ord,]
  }
  v@ordered = TRUE

  ll = matrix(as.numeric(NA), nrow=n, ncol=n)
  istrat = calcistrat(v)

  pwh = p0[wh]
  psteps = matrix(as.numeric(NA), ncol=nrow(wh), nrow=n)
  for(i in 1:ncol(psteps))
    psteps[,i] = switch(as.integer(wh[i, 3]),
            { delta = quantile(x[, wh[i, 2]], probs=0.05)
              pwh[i] + seq(-delta, +delta, length=n) },
            { pwh[i] * exp(seq(-log(sqrt(2)), log(sqrt(2)), length=n)) },
             stop("Zapperlot"))
 
  for(j in seq_len(n)) {
    p0[wh[1,,drop=FALSE]] = psteps[j,1]
    for(i in seq_len(n)) {
      p0[wh[2,,drop=FALSE]] = psteps[i,2]
      ll[i,j] = .Call("vsn2_ll", x, as.vector(p0), istrat, PACKAGE="vsn")[1]
    }
  }
  
  if(require("lattice"))
    print(levelplot(z~p1*p2,
      data=data.frame(z=as.vector(ll), p1=psteps[col(ll), 1], p2=psteps[row(ll), 2]),
      xlab=sprintf("par[%d,%d,%d]", wh[1,1], wh[1,2], wh[1,3]),
      ylab=sprintf("par[%d,%d,%d]", wh[2,1], wh[2,2], wh[2,3])))
  else
    warning("lattice package is not available, cannot call 'levelplot'.")
  
  return(invisible(ll))
}    
