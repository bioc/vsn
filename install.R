require(Biobase)

## register vsn as a normalization method with the affy package, if that is loaded:
if ("package:affy" %in% search())
  if(!"vsn" %in% get("normalize.AffyBatch.methods", "package:affy"))
    assign("normalize.AffyBatch.methods",
           c(get("normalize.AffyBatch.methods", pos="package:affy"), "vsn"),
           pos="package:affy")
