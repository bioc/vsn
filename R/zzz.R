.First.lib <- function(lib, pkgname, where) {
  ## load the compiled code
  library.dynam(pkgname, pkgname, lib)
  if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("vsn")
    }
}
