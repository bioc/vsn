require(Biobase)

## register vsn as a normalization method with the affy package, if that is loaded:
if(exists("normalize.AffyBatch.methods"))
  if(!"vsn" %in% normalize.AffyBatch.methods)
    normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods, "vsn")
