.onLoad <- function(libname, pkgname) {
  # load helper functions into package namespace at load time
  ns <- asNamespace(pkgname)

  assign("PhD.m.est", safe_get("PhD.m.est", "iNEXT.3D"), envir = ns)
  assign("phyBranchAL_Abu", safe_get("phyBranchAL_Abu", "iNEXT.3D"), envir = ns)
  assign("coverage_to_size", safe_get("coverage_to_size", "iNEXT.beta3D"), envir = ns)
}
