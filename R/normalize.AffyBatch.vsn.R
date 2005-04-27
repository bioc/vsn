##------------------------------------------------------------
## a wrapper for vsn to be used as a normalization method in
## the package affy
##------------------------------------------------------------
normalize.AffyBatch.vsn = function (abatch, subsample = 20000, niter = 4, strata,...)  
{
  if (!("package:affy" %in% search())) 
    stop("Please load the package affy before calling normalize.AffyBatch.vsn")

  ## ind = the perfect match probes. If length(ind) is larger than the value in
  ## subsample, then only use a random sample of size subsample from these
  ind = unlist(indexProbes(abatch, "pm"))
  if (!is.na(subsample)) {
    if (!is.numeric(subsample)) 
      stop("'subsample' must be numeric.")
    if (length(ind) > subsample) 
      ind = sample(ind, subsample)
  }

  if(missing(strata))
    strind <- rep(as.integer(1), length(ind))
  else 
    strind <- strata[ind]
  
  ## call parameter estimation (on subset of data)
  vsnres = vsn(intensity(abatch)[ind, ], niter=niter, strata=strind, ...)

  ## add parameters to preprocessing slot
  pars = preproc(description(vsnres))$vsnParams
  parsIter = preproc(description(vsnres))$vsnParamsIter # J.: added 2005/04/26
  description(abatch)@preprocessing = c(description(abatch)@preprocessing, list(vsnParams=pars, vsnParamsIter=parsIter))

  ## apply the transformation (to all data)
  intensity(abatch) = exp(vsnh(intensity(abatch), pars, strata=strata))
  return(abatch)
}


