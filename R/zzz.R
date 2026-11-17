#' @keywords internal
#' @import symcompR
"_PACKAGE"

# symcompR memoizes Add, Mul and Pow on the keys of their arguments. The
# variables of modelcompR carry their metadata as an attribute that is not
# part of the key, so a cached x^2 could hold an earlier declaration of x
# (another default, bounds or input flag). The cache is turned off, as
# symcompR itself does when two symbols share a key but not their
# assumptions.
.onLoad <- function(libname, pkgname) {
  .arith$native <- identical(getOption("modelcompR.arithmetic"), "native")
  st <- tryCatch(utils::getFromNamespace(".ctor_state", "symcompR"), error = function(e) NULL)
  if (is.environment(st)) assign("off", TRUE, envir = st)
  invisible()
}
