##-----------------------------------------------------------------
## .First.lib: this function is called when the package is loaded
##-----------------------------------------------------------------
.First.lib <- function(lib, pkgname, where) {
  
  if(missing(where)) {
    where <- match(paste("package:", pkgname, sep=""), search())
    if(is.na(where)) {
      warning(paste("Not a package name: ",pkgname))
      return()
    }
    where <- pos.to.env(where)
  }
  .initvsn(where)
}

