plotVsnLogLik = function(object, p, whichp=1:2, expand=1, ngrid=31L, fun=logLik, ...) {

  stopifnot(length(whichp)==2L, ngrid>=2L)

  if(length(expand)==1)
    expand=rep(expand,2)

  d = ncol(object)
  stopifnot(2L*d*nlevels(object@strata)==length(p))

  psteps = sapply(1:2, function(k) {
    z = ((whichp[k]-1) %/% nlevels(object@strata))
    i = (z %% d) + 1L
    aorb = z %/% d
    stopifnot(aorb %in% c(0L,1L))
    if(aorb==0) {
      p[whichp[k]] + seq(-expand[1L], +expand[1L], length=ngrid)
    } else {
      p[whichp[k]] + seq(-expand[2L], +expand[2L], length=ngrid)
    }
  })

  pgrid = expand.grid(psteps[,1], psteps[,2], KEEP.OUT.ATTRS = FALSE)
  psamp = matrix(p, nrow=length(p), ncol=nrow(pgrid))
  for(i in 1:2)
    psamp[whichp[i], ] = pgrid[, i]

  ll = fun(object, psamp, ...)
  pgrid$logLik = ll[1, ]

  print(levelplot(logLik ~ Var1*Var2, data=pgrid,
       col.regions=colorRampPalette(c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026"))(128)))  ## colors are from RColorBrewer: brewer.pal(9,"YlOrRd")
  
  return(invisible(pgrid))
}    
