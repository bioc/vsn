.onAttach <- function(libname, pkgname) {
   addVigs2WinMenu("vsn")
}

.onLoad <- function(libname, pkgname) {
  ## register vsn as a normalization method with the affy package, if that is loaded:
  if ("package:affy" %in% search())
    if(!"vsn" %in% affy::normalize.AffyBatch.methods())
      upDate.normalize.AffyBatch.methods(
             c(affy::normalize.AffyBatch.methods(), "vsn"))
}
