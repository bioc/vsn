## 'manual' dispatch on the type of first argument
## this is quite clumsy, and obsolete now
getIntensityMatrix = function(intensities, verbose) {
  
  if(is.matrix(intensities)) {
    if (!is.numeric(intensities))
      stop("'intensities' was found to be a matrix, but is not numeric.")
    y = intensities
  } else if (is.data.frame(intensities)) {
    if (!all(sapply(intensities, is.numeric)))
      stop("'intensities' was found to be a data.frame, but contains non-numeric columns.")
    y = as.matrix(intensities)
  } else if (is(intensities, "ExpressionSet")) {
    y = exprs(intensities)
  } else if (is(intensities, "marrayRaw")) {
    nrslides = as.integer(ncol(intensities@maRf))
    nrspots  = as.integer(nrow(intensities@maRf))
    if (verbose)
      cat(sprintf("Converting marrayRaw (%d spots, %d slides) to %dx%d matrix.\n",
                  nrspots, nrslides, nrspots, as.integer(2*nrslides)),
          "Gf-Gb in odd columns, Rf-Rb in even columns.\n")
    y = matrix(NA, nrow=nrspots, ncol=2*nrslides)
    y[, (1:nrslides)*2-1 ] = intensities@maGf - intensities@maGb
    y[, (1:nrslides)*2   ] = intensities@maRf - intensities@maRb
  } else
    stop(sprintf("Don't know how to handle class %s of 'intensities'.", class(intensities)))
  ## end of cascaded ifs

  if (any(is.na(y)))
    stop(paste("'intensities' must not contain NA values.\n",
             "This could indicate that the input data has already undergone some\n",
             "thresholding or transformation (log?), and may not satisfy the\n",
             "requirements of the multiplicative-additive noise model.\n",
             "If you are sure that it is meaningful to proceed, please\n",
             "consider calling vsn on a subset of data where all values\n",
             "are defined, and then use vsnh on the full set of data.\n"))

  if (ncol(y)<=1) 
    stop(paste("'intensities' must be a matrix with at least two columns.\n",
               "Please read the documentation and the paper\n",
               "(Huber et al., Bioinformatics 18 (2002) S96-S104).\n"))

  storage.mode(y) = "double"
  return(y)
}
     
