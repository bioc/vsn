rowSds = function(x, ...) {
  sqrt(rowSums(sqr(x-rowMeans(x, ...)), ...)/(ncol(x)-1))
}
