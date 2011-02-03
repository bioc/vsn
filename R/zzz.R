.onAttach <- function(libname, pkgname) {
   addVigs2WinMenu("vsn")
}

.onLoad <- function(libname, pkgname) {
  ## register vsn as a normalization method with the affy package
  if(!("vsn" %in% normalize.AffyBatch.methods()))
    upDate.normalize.AffyBatch.methods(c(normalize.AffyBatch.methods(), "vsn"))
}
