##---------------------------------------------------------
## Enumerate all the subsets of size k of the integers 1:n.
## The result is returned in a matrix with k rows and
## (n choose k) columns
## This function is not needed for VSN, but it is nice to
## have it for the examples / the vignette and I haven't
## yet found it anywhere else.
##---------------------------------------------------------
nchoosek = function(n, k) {
  if (!is.numeric(n)||!is.numeric(k)||is.na(n)||is.na(k)||length(n)!=1||length(k)!=1)
    stop("arguments must be non-NA numeric scalars.")
  if (k>n||k<0)
    stop("Arguments must satisfy 0 <= k <= n.")

  nck = choose(n, k)
  res = matrix(NA, nrow=k, ncol = nck)
  res[, 1] = 1:k
  j = 2
  repeat {
    res[, j] = res[, j-1]
    i = k
    repeat {
      res[i,j] = res[i,j]+1
      if(res[i,j] <= n-(k-i))
        break
      i = i-1
      stopifnot(i>=1)
    }
    if (i<k)
       res[(i+1):k,j] = res[i,j] + 1:(k-i)
    j = j+1
    if (j>nck) break
  }
  ## plausibility tests
  stopifnot(all(res[, nck]==(n-k+1):n))
  stopifnot(all(res<=n) && all(res>=1))
  return(res)
}

