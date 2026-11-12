# The reference runs of most tests come from Julia, Python or C and were
# matched to the last bit on arm64 macOS, where R's long double is the double
# type and the C library's elementary functions round as the probes saw them.
# Other platforms differ in the last bits (x87 long double in sum() and the
# parser, other libm), and long iterations amplify that. There the tests
# compare results to a tolerance and leave out iterate trajectories and
# evaluation counts.
ref_platform <- function() {
  Sys.info()[["sysname"]] == "Darwin" && grepl("^(aarch64|arm64)", R.version$arch)
}

# the tolerance of a reference comparison: `tol` on the reference platform,
# `loose` elsewhere
ref_tol <- function(tol, loose = 1e-6) if (ref_platform()) tol else loose
