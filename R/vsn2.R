##----------------------------------------------------------------------
## Robust calibration and variance stabilization
##   revised version vsn2
## (C) Wolfgang Huber <huber@ebi.ac.uk> 2002-2007
## contributions from Markus Ruschhaupt, Dennis Kostka, David Kreil
##
## The computations are organised as a set of nested function calls.
## The innermost is vsnML, which calls the C code.
## The next one is vsnLTS, which does the robust modification of the ML.
## The next one is vsnStrata, which if necessary reorders the data.
## The next one is vsnSample, which allows for sampling.
## The outermost one is vsnMatrix
##----------------------------------------------------------------------

##----------------------------------------------------------------------------
## ML estimator. If the C optimization code fails to converge, retry with 
## different start parameters, for niter times.
##-----------------------------------------------------------------------------
vsnML = function(v, niter=4) {

  p = as.vector(v@pstart)
  istrat = calcistrat(v) ## pointers to the starts of strata
  
  for (iter in seq_len(niter)) {
    o = .Call("vsn2_optim", v@x, p, istrat, v@reference@refh, v@reference@refsigma, PACKAGE="vsn")
    
    conv = as.integer(o[length(o)])
    if (conv==0) 
      break
    
    if(v@verbose)
      cat("\nCONV=", conv, ":\n",
          "p=", paste(signif(p, 5), collapse=", "), "\n",
          "o=", paste(signif(o[-length(o)], 5), collapse=", "), "\n",
          "Restarting with new initial parameters.\n", sep="")

    ## The most frequent non-zero value that can occur for conv seems to be 52,
    ## which codes for "ABNORMAL_TERMINATION_IN_LNSRCH".
    ## This seems to indicate that a stepwidth to go along the gradient could 
    ## not be found, probably because the start point p was already too close 
    ## to the optimum (?). Hence, try again from a slightly different start point.
    r1 = exp(runif(length(p), min=-0.03, max=0.03))
    r2 = runif(length(p), min=0, max=0.02)
    p = p*r1 + r2 
    
  } 

  if(conv!=0)
    stop("\n", paste(strwrap(c("The likelihood optimization did not converge. A possible",
         "reason is that the normalization parameters are not identifiable",
         "from the provided data. For example, this will be the case if the",
         "columns of the data matrix are exactly co-linear or affine dependent.",
         "Please verify the data to make sure there were no mix-ups.")), collapse="\n"), "\n")
  
  par =  array(o[-length(o)], dim=dim(v@pstart))
  if (any(par[,,2]<0))
    stop("Likelihood optimization produced negative parameter estimates.\n",
         "Please contact the package maintainer.\n")

  return(par)
}

##------------------------------------------------------------
## istrat[j] is the starting position for stratum j
## i.e. istrat[j]...istrat[j+1] are the elements of ysel that
## belong to stratum j (using C indexing convention,
## i.e. starting at 0).
## Note: the counting over the different samples is folded into j,
## i.e., if there are 4 strata on the array and 3 samples, then
## j runs from 1:12
##------------------------------------------------------------
calcistrat = function(vp) {
  nrs = nlevels(vp@strata)
  one = as.integer(1)

  if(length(vp@strata)>0) {
    stopifnot(vp@ordered)
    istr = which(!duplicated(vp@strata))-one
  } else {
    istr = as.integer(0)
  }
  stopifnot(length(istr)==nrs)

  istrat = integer(ncol(vp@x)*nrs+one)
  for(i in 0:(ncol(vp@x)-one))
    istrat[i*nrs + seq_len(nrs)] = i*nrow(vp@x) + istr

  istrat[length(istrat)] = nrow(vp@x)*ncol(vp@x)  ## end point
  return(istrat)
}
         
##-------------------------------------------------------------------------
## LTS robust modification of the ML estimator
##-------------------------------------------------------------------------
vsnLTS = function(v) {

  ## a place to save the trajectory of estimated parameters along the iterations
  params = array(as.numeric(NA), dim=c(dim(v@pstart), v@cvg.niter))

  ## for calculating a convergence criterion: earlier result
  oldhy   = Inf
  
  ## the number of iterations that have already met the convergence criterion
  cvgcCnt = 0

  ## integer version of "v@strata"
  intStrat = if(length(v@strata)==0) rep(as.integer(1), nrow(v@x)) else as.integer(v@strata)

  for(iter in seq_len(v@cvg.niter)) {
    sv = if(iter==1) v else v[sel, ]
    par = vsnML(sv)

    ## apply to all data
    hy = vsn2trsf(v@x, par, intStrat)
    params[,,,iter] = v@pstart = par

    ## if LTS.quantile is 1, then the following stuff is not necessary
    if (abs(v@lts.quantile-1)<sqrt(.Machine$double.eps))
      break
       
    ## Calculate residuals
    hmean  = if(length(v@reference@refh)>0) {
      v@reference@refh  ## WITH reference
    } else {
      rowMeans(hy)      ## WITHOUT reference
    }
    
    sqres  = hy - hmean
    sqres  = rowSums(sqres*sqres) 

    ## select those data points whose sqres is within the quantile; do this separately
    ## within each stratum, and also within strata defined by hmean
    ## (see the SAGMB 2003 paper for details)
    nrslice = 5
    slice   = ceiling(rank(hmean)/length(hmean)*nrslice)
    slice   = factor((intStrat-1)*nrslice + slice)
    grmed   = tapply(sqres, slice, quantile, probs=v@lts.quantile)
    meds    = grmed[as.character(slice)]
    stopifnot(!any(is.na(meds)))
    sel     = (sqres <= meds)

    ## diagnostic plot (most useful for d=2)
    ## plot(hy, pch=".", col=2-sel, main=sprintf("iter=%d", iter))
    ## browser()
    
    ## Convergence check
    ## after a suggestion from David Kreil, kreil@ebi.ac.uk
    if(v@cvg.eps>0) {
      cvgc    = max(abs(hy - oldhy))
      cvgcCnt = if(cvgc<v@cvg.eps) (cvgcCnt+1) else 0 
      if (v@verbose)
        cat(sprintf("iter %2d: cvgc=%.5f%%, par=", iter, cvgc),
            sprintf("%9.3g",par), "\n")
      if(cvgcCnt>=3)
        break
      oldhy = hy
    }

    if(v@verbose)
      cat("\b\b\b\b\b\b\b\b\b",
       sprintf("%2d", as.integer(round((1-((iter-v@cvg.niter)/v@cvg.niter)^2)*100))),
          "% done.", sep="")

  } ## end of for-loop
  if(v@verbose) cat("\n")
  return(list(par=par, hy=hy, params=params))
}


##----------------------------------------------------------------------------------
## vsnStrata: if necessary,
## reorder the rows of the matrix so that each stratum sits in a contiguous block
##----------------------------------------------------------------------------------
vsnStrata = function(v, minDataPointsPerStratum=42) {

  ord = NULL
  if(length(v@strata)>0) {
    if(nlevels(v@strata)>1) {
      if(any(table(v@strata) < minDataPointsPerStratum))
        warning("*** There are less than ", minDataPointsPerStratum, " data points in ",
                "some of the strata.\n*** The fitted parameters may be unreliable.\n",
                "*** You could try to reduce the number of strata.\n")
    
      ## reorder the rows of the matrix so that each stratum sits in a contiguous block
      ord = order(v@strata)
      v   = v[ord,]
    }
    v@ordered = TRUE
  }

  res = vsnLTS(v)

  if(!is.null(ord)) {
    res$hy[ord, ] = res$hy
  }
  return(res)
} 

##----------------------------------------------------------------------------------
## vsnSample: allows for (balanced) subsetting
##----------------------------------------------------------------------------------
vsnSample = function(v) {

  if(v@subsample>0) {
    if((length(v@strata)>0) && (nlevels(v@strata)>1)) { 
      sp = split(seq(along=v@strata), v@strata)
      wh = unlist(lapply(sp, sample, size=v@subsample))
    } else {
      wh = sample(length(v@strata), size=v@subsample)
    }
    v = v[wh, ]
  } 
  vsnStrata(v)
}

##----------------------------------------------------------------------
## vsnMatrix
##----------------------------------------------------------------------
vsnMatrix = function(x,
  reference,
  strata,
  lts.quantile = 0.9,
  subsample    = as.integer(0),
  verbose      = interactive(),
  returnData   = TRUE,
  pstart,
  cvg.niter    = as.integer(7),
  cvg.eps      = 1e-3) {

  if(missing(pstart)) {
    pstart = array(0, dim=c(nlevels(strata), ncol(x), 2))
    pstart[,,2] = rep(1/apply(x, 2, IQR), each=dim(pstart)[1])
  }
  if(missing(reference)) {
    reference = new("vsn")
  } else {
    if(nrow(reference)!=nrow(x))
      stop("'nrow(reference)' must be equal to 'nrow(x)'.")
    if(nrow(reference)!=length(reference@refh))
      stop(sprintf("The slot 'reference@refh' has length 0, but expected is n=%d", nrow(reference)))
  }
  
  if(missing(strata)) {
    strata = factor(integer(0), levels="all")
  }
  
  v = new("vsnInput",
    x      = x,
    strata = strata,
    pstart = pstart,
    reference = reference,
    returnData = returnData,
    lts.quantile = lts.quantile,
    cvg.niter  = cvg.niter,
    cvg.eps   = cvg.eps,
    subsample = subsample,
    verbose   = verbose,
    ordered   = FALSE)
  
  ## Print welcome message
  if (verbose)
    cat("vsn: ", nrow(x), " x ", ncol(x), " matrix (", nlevels(strata), " strat",
        ifelse(nlevels(strata)==1, "um", "a"), ").  0% done.", sep="")

  res = vsnSample(v)
  
  return(res)
}


##---------------------------------------------------------------------
## The "arsinh" transformation
##
## Note: the constant -log2(2*facs[1]) is added to the transformed data
## in order to achieve h_1(y) \approx log2(y) for y\to\infty, that is,
## better comparability to the log2 transformation.
## It has no effect on the generalized log-ratios.
##--------------------------------------------------------------------
vsn2trsf = function(y, p, strata) {
  if (!is.matrix(y) || !is.numeric(y))
    stop("'y' must be a numeric matrix.\n")
  
  if (!is.array(p) || !is.numeric(p) || any(is.na(p)))
    stop("'p' must be an array with no NAs.\n")

  if(missing(strata)) {
    strata <- rep(as.integer(1), nrow(y))
  } else {
    if(!is.integer(strata) || !is.vector(strata) || 
       length(strata)!=nrow(y) || any(is.na(strata)))
      stop("'strata' must be an integer vector of length nrow(y) with no NAs.")
  }
  nrstrata <- max(strata)
  
  if(nrstrata==1 && length(dim(p))==2)
    dim(p) <- c(1, dim(p))
  
  if(length(dim(p))!=3 || dim(p)[1]!=nrstrata || dim(p)[2]!=ncol(y) || dim(p)[3]!=2)
    stop("'p' has wrong dimensions.")
  if (any(p[,,2]<=0))
    stop("'p' contains invalid values: factors must be non-negative.")

  hy = .Call("vsn2_trsf", y, as.vector(p), strata, PACKAGE="vsn")

  dimnames(hy) = dimnames(y)
  return(hy)
}

