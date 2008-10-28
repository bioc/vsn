.onAttach <- function(libname, pkgname) {
   addVigs2WinMenu("vsn")
}

.onLoad <- function(libname, pkgname) {
  ## register vsn as a normalization method with the affy package, if that is loaded:
  if ("package:affy" %in% search())
    assign("normalize.AffyBatch.methods", pos="package:affy", value=c("vsnrma", get("normalize.AffyBatch.methods", pos="package:affy")))
}
