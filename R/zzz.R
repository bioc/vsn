##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib = function(lib, pkgname, where) {
  require(Biobase, quietly=TRUE) || stop("Cannot load without package \"Biobase\"")
}

