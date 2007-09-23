justvsn = function(x,  ...) {
    fit = vsn2(x, ...)
    predict(fit, newdata=x, useDataInFit=TRUE)
  }

vsnrma = function(x, ...) {
  fit = vsn2(x, ...)
  exprs(x) = 2^exprs(fit)
  
  ## call RMA to do the probeset summarization,
  ## (with no background correction / normalization, that is already done)
  rv = rma(x, normalize=FALSE, background=FALSE)
  if(!(is(rv, "ExpressionSet") && validObject(rv)))
    stop("Failed to create a valid 'ExpressionSet' object.")
  return(rv)
}


