rowSds = function(x, ...) {
  sqr = function(a) a*a  ## faster than a^2
  sqrt(rowSums(sqr(x-rowMeans(x, ...)), ...)/(ncol(x)-1))
}
