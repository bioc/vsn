rowSds = function(x, ...) {
  sqr     = function(a) a*a  ## faster than a^2
  n       = rowSums(!is.na(x))
  n[n<=1] = NA
  return(sqrt(rowSums(sqr(x-rowMeans(x, ...)), ...)/(n-1)))
}
