plotVsnLogLik = function(object, p, whichp=1:2, expand=1, ngrid=31, ...) {

  stopifnot(length(whichp)==2)

  if(length(expand)==1)
    expand=rep(expand,2)

  d = ncol(object)
  stopifnot(2*d*nlevels(object@strata)==length(p))

  psteps = sapply(1:2, function(k) {
    z = ((whichp[k]-1) %/% nlevels(object@strata))
    i = (z %% d) + 1
    aorb = z %/% d
    stopifnot(aorb %in% c(0,1))
    if(aorb==0) {
      delta = quantile(object@x[,i], probs=0.045)*expand[1]
      p[whichp[k]] + seq(-delta, +delta, length=ngrid)
    } else {
      delta = log(sqrt(2))*expand[2]
      p[whichp[k]] * exp(seq(-delta, delta, length=ngrid))
    }
  })

  pgrid = expand.grid(psteps[,1], psteps[,2], KEEP.OUT.ATTRS = FALSE)
  psamp = matrix(p, nrow=length(p), ncol=nrow(pgrid))
  for(i in 1:2)
    psamp[whichp[i], ] = pgrid[, i]

  ll = logLik(object, psamp, ...)
  pgrid$logLik = ll[1, ]

  pkgs = c("lattice", "RColorBrewer")
  if(all(sapply(pkgs, function(p) do.call("require", args=list(p, character.only=TRUE)))))
    print(levelplot(logLik ~ Var1*Var2, data=pgrid,
                    col.regions=colorRampPalette(brewer.pal(9,"YlOrRd"))(128)))
  else
    warning("packages ", paste(pkgs, collapse=", "), " are not available, cannot call 'levelplot'.")
  
  return(invisible(pgrid))
}    
