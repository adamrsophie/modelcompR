# No exported name may mask a function of the packages R attaches by default.
# Names shared with symcompR (flatten, compose, and Differential, extend and
# System in symcompR versions with differential geometry, unit systems and
# multibody mechanics) are allowed because the modelcompR functions hand
# foreign arguments to symcompR.
default_pkgs <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets")

test_that("exports do not mask base R functions", {
  ns <- getNamespaceExports("modelcompR")
  for (p in default_pkgs) {
    expect_identical(intersect(ns, getNamespaceExports(p)), character(0), info = p)
  }
})

test_that("exports do not mask symcompR functions except delegating ones", {
  ns <- getNamespaceExports("modelcompR")
  hit <- setdiff(intersect(ns, getNamespaceExports("symcompR")), c("flatten", "compose", "Differential", "extend", "System"))
  expect_identical(hit, character(0))
})

test_that("delegating functions keep symcompR behavior", {
  x <- symcompR::Symbol("x")
  expect_equal(format(compose(x^2, x + 1)), format(symcompR::compose(x^2, x + 1)))
  expect_equal(flatten(list(1, list(2, 3))), symcompR::flatten(list(1, list(2, 3))))
  if ("extend" %in% getNamespaceExports("symcompR")) {
    # a non-system goes to symcompR's extend() generic
    expect_error(extend(1, 2), "no applicable method for 'extend'")
  }
  if ("System" %in% getNamespaceExports("symcompR")) {
    # SymPy's multibody System: the constructor, $ methods, names and print
    # stay symcompR's
    ms <- System()
    expect_s3_class(ms, "physics_obj")
    q <- symcompR::dynamicsymbols("q")
    ms$add_coordinates(q)
    ref <- symcompR::System()
    ref$add_coordinates(q)
    expect_identical(capture.output(print(ms$q)), capture.output(print(ref$q)))
    expect_setequal(names(ms), ls(ms, all.names = TRUE))
    expect_output(print(ms))
    N <- symcompR::ReferenceFrame("N")
    expect_s3_class(System(N), "physics_obj")
    expect_s3_class(System(frame = N), "physics_obj")
  }
})

test_that("internal definitions do not shadow base R functions", {
  skip_if_not(file.exists(test_path("../../tools/check_masking.R")))
  out <- system2(file.path(R.home("bin"), "Rscript"), test_path("../../tools/check_masking.R"),
                 stdout = TRUE, stderr = TRUE, env = paste0("R_LIBS=", paste(.libPaths(), collapse = ":")))
  expect_match(paste(out, collapse = "\n"), "No function masks a base R function")
})
