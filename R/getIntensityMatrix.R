## Extract the intensity matrix from the argument "intensities"
## Probably this should, and may in a future version be done via S3- or S4
## method dispatching.
getIntensityMatrix = function(intensities, verbose) {
  y = switch(class(intensities),
    matrix     = { if (!is.numeric(intensities))
                     stop("'intensities' was found to be a matrix, but is not numeric.")
                   intensities
                 },
    data.frame = {  if (!all(sapply(intensities, is.numeric)))
                      stop("'intensities' was found to be a data.frame, but contains non-numeric columns.")
                    as.matrix(intensities)
                  },
    exprSet    = { exprs(intensities)
                 },
    marrayRaw  = { nrslides = as.integer(ncol(intensities@maRf))
                   nrspots  = as.integer(nrow(intensities@maRf))
                   if (verbose)
                     cat(sprintf("Converting marrayRaw (%d spots, %d slides) to %dx%d matrix.\n",
                                 nrspots, nrslides, nrspots, as.integer(2*nrslides)),
                         "Gf-Gb in odd columns, Rf-Rb in even columns.\n")
                   tmp = matrix(NA, nrow=nrspots, ncol=2*nrslides)
                   tmp[, (1:nrslides)*2-1 ] = intensities@maGf - intensities@maGb
                   tmp[, (1:nrslides)*2   ] = intensities@maRf - intensities@maRb
                   tmp
                 },
    stop(paste("'intensities' has class ", class(intensities),
         ". Permitted are: matrix, data.frame, exprSet, marrayRaw", sep=""))
  )  ## end of switch statement

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

  if(is.integer(y))
    y <- y+0.0 ## convert to double

  return(y)
}
     
