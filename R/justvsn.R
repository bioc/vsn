justvsn=function(x, ...) {
  fit = vsn2(x, ...)
  exprs(x) = fit@hx
  return(x)
}