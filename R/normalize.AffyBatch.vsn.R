##------------------------------------------------------------
## a wrapper for vsn to be used as a normalization method in
## the package affy
##------------------------------------------------------------
normalize.AffyBatch.vsn = function (abatch, subsample = 20000, niter = 4, ...) 
{
  if (!("package:affy" %in% search())) 
    stop("Please load the package affy before calling normalize.AffyBatch.vsn")

  ## ind = the perfect match probes. If more than subsample, then only use
  ## a random sample of size subsample from these
  ind = unlist(indexProbes(abatch, "pm"))
  if (!is.na(subsample)) {
    if (!is.numeric(subsample)) 
      stop("Argument \"subsample\" must be numeric.")
    if (length(ind) > subsample) 
      ind = sample(ind, subsample)
  }
  
  ## call parameter estimation (on subset of data)
  vsnres = vsn(intensity(abatch)[ind, ], niter=niter, ...)

  ## apply the transformation (to all data)
  intensity(abatch) = exp(vsnh(intensity(abatch), preproc(description(vsnres))$vsnParams))
  return(abatch)
}

