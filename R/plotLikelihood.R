plotLogLik = function(x, strata, p0, which, delta, n=31) {

  stopifnot(length(which)==2, length(delta)==2)
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

  px = (2*(0:(n-1))/(n-1)-1)*delta[1] + p0[which[1]]
  py = (2*(0:(n-1))/(n-1)-1)*delta[2] + p0[which[2]]
  for(j in seq_len(n)) {
    p0[which[1]] = px[j]
    for(i in seq_len(n)) {
      p0[which[2]] = py[i]
      ll[i,j] = .Call("vsn2_c", x, p0, istrat, as.integer(1), PACKAGE="vsn")[1]
    }
  }
  
  if(require("lattice"))
    print(levelplot(z~x*y, data=data.frame(z=as.vector(ll), x=px[col(ll)], y=py[row(ll)])))
  else
    warning("lattice package is not available, cannot call 'levelplot'.")
  browser()
  
  return(invisible(ll))
}    
