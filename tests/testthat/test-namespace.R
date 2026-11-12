# No exported name may mask a function of the packages R attaches by default.
# Names shared with symcompR (flatten, compose, and Differential and extend in
# symcompR versions with differential geometry and unit systems) are allowed
# because the modelcompR functions hand foreign arguments to symcompR.
default_pkgs <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "datasets")

test_that("exports do not mask base R functions", {
  ns <- getNamespaceExports("modelcompR")
  for (p in default_pkgs) {
    expect_identical(intersect(ns, getNamespaceExports(p)), character(0), info = p)
  }
})

test_that("exports do not mask symcompR functions except delegating ones", {
  ns <- getNamespaceExports("modelcompR")
  hit <- setdiff(intersect(ns, getNamespaceExports("symcompR")), c("flatten", "compose", "Differential", "extend"))
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
})

test_that("internal definitions do not shadow base R functions", {
  skip_if_not(file.exists(test_path("../../tools/check_masking.R")))
  out <- system2(file.path(R.home("bin"), "Rscript"), test_path("../../tools/check_masking.R"),
                 stdout = TRUE, stderr = TRUE, env = paste0("R_LIBS=", paste(.libPaths(), collapse = ":")))
  expect_match(paste(out, collapse = "\n"), "No function masks a base R function")
})
