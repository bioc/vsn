##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib = function(lib, pkgname, where) {
  require(Biobase, quietly=TRUE) || stop("Cannot load without package \"Biobase\"")
  cat("vsn: The return type of the function vsn has changed since version 1.2.0, please see its help page.\n")
}

