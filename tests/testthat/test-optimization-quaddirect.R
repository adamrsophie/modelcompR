# QuadDirect: QuadDIRECT.jl through OptimizationQuadDIRECT, ported to R.
# The method is deterministic; every evaluation of the runs below matches
# Julia bit for bit (numbers from tools/quaddirect-probe.jl).

qd_fun <- function(name) get(name, envir = asNamespace("modelcompR"))

qd_run <- function(f, splits, lower, upper, ...) {
  n <- 0L
  r <- qd_fun(".qd_analyze")(function(x) {
    n <<- n + 1L
    f(x)
  }, splits, lower, upper, ...)
  b <- qd_fun(".qd_minimum")(r$root)
  list(u = qd_fun(".qd_position")(b, r$x0), o = qd_fun(".qd_value")(b), n = n, root = r$root, x0 = r$x0)
}

test_that("QuadDirect runs match QuadDIRECT.jl", {
  # tools/quaddirect-probe.jl: the number of evaluations and the minimum
  canyon <- function(x) 0.1 * (x[1] + x[2])^2 + 10 * (x[1] - x[2])^2
  r <- qd_run(canyon, list(c(-11, -10, -9), c(-7, -6, -5)), c(-Inf, -Inf), c(Inf, Inf), atol = 0.01)
  expect_identical(r$n, 30L)
  expect_equal(r$u, c(-6.901146321069973e-13, -6.918909889463976e-13), tolerance = 1e-15)
  expect_equal(r$o, 1.910255080985268e-25, tolerance = 1e-15)
  oned <- function(x) (x[1] - 0.7)^2 * (x[1] + 0.3)^2 + 0.1 * x[1]
  r <- qd_run(oned, list(c(-1, 0, 1)), -2, 2)
  expect_identical(r$n, 236L)
  expect_equal(r$u, -0.34401695602142973, tolerance = 1e-15)
  expect_equal(r$o, -0.03228988427082661, tolerance = 1e-15)
  camel <- function(x) (4 - 2.1 * x[1]^2 + x[1]^2 * x[1]^2 / 3) * x[1]^2 + x[1] * x[2] + (-4 + 4 * x[2]^2) * x[2]^2
  r <- qd_run(camel, list(c(-2, 0, 2), c(-1, 0, 1)), c(-3, -2), c(3, 2))
  expect_identical(r$n, 430L)
  expect_equal(r$u, c(0.08984261066583982, -0.712657381571297), tolerance = 1e-15)
  expect_equal(r$o, -1.0316284534812283, tolerance = 1e-15)
})

test_that("QuadDirect stores the value of every evaluated point in its tree", {
  # tools/quaddirect-probe.jl, case edge: QuadDIRECT gives a box value to the
  # wrong child at its 79th evaluation; the port does not, and finds the
  # same minimum
  H <- matrix(c(2.2033, -0.8346, -0.1717, -0.8346, 1.6376, -0.1638, -0.1717, -0.1638, 1.393), 3)
  ctr <- c(0.5042, 1.2136, -0.5949)
  edge <- function(x) {
    d <- x - ctr
    s <- 0
    for (i in 1:3) for (j in 1:3) s <- s + d[i] * (H[i, j] * d[j])
    s
  }
  r <- qd_run(edge, rep(list(c(-1, 0, 1)), 3), rep(-1, 3), rep(1, 3), maxevals = 300)
  expect_equal(r$u, c(0.42052827743518195, 1, -0.63033008956524), tolerance = 1e-15)
  expect_equal(r$o, 0.05855961301783394, tolerance = 1e-15)
  nodes <- Filter(function(b) b$splitdim != 0L, qd_fun(".qd_preorder")(r$root))
  stored <- unlist(lapply(nodes, function(b) b$fvalues))
  actual <- unlist(lapply(nodes, function(b) vapply(b$children, function(ch) edge(qd_fun(".qd_position")(ch, r$x0)), 0)))
  expect_identical(stored, actual)
})

test_that("solve with QuadDirect follows OptimizationQuadDIRECT", {
  rosen <- function(x, p) (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
  prob <- OptimizationProblem(rosen, c(0, 0), c(1, 100), lb = c(-1, -1), ub = c(0.8, 0.8))
  s2 <- list(c(-0.5, 0, 0.5), c(-0.5, 0, 0.5))
  # the OptimizationQuadDIRECT test
  sol <- solve(prob, QuadDirect(), splits = s2)
  expect_lt(10 * sol$objective, rosen(c(0, 0), c(1, 100)))
  # tools/quaddirect-probe.jl, cases rosen2 and rosen_maxevals
  expect_equal(sol$u, c(0.7972222222222222, 0.6339364516762053), tolerance = 1e-15)
  expect_equal(sol$objective, 0.04138348146854611, tolerance = 1e-15)
  expect_identical(sol$stats$fevals, 1458L)
  expect_equal(sol$retcode, "Default")
  sol <- solve(prob, QuadDirect(), splits = s2, maxiters = 100, reltol = 1e-6)
  expect_equal(sol$u, c(0.6173653116315079, 0.3865517241379306), tolerance = 1e-15)
  expect_equal(sol$objective, 0.14933805848035295, tolerance = 1e-15)
  expect_identical(sol$stats$fevals, 112L)
  expect_true(is.environment(sol$original) && sol$original$splitdim == 1L)

  expect_error(solve(prob, QuadDirect()), "You must provide the initial locations at which to evaluate the function in `splits`")
  expect_warning(solve(prob, QuadDirect(), splits = s2, maxiters = 20, maxtime = 10),
                 "common maxtime is currently not used by QuadDirect()", fixed = TRUE)
  expect_error(solve(prob, QuadDirect(), splits = s2, maxiters = 0),
               "The number of maxiters has to be a non-negative and non-zero number.", fixed = TRUE)
  expect_error(solve(prob, QuadDirect(), splits = list(c(0.5, 0, -0.5), c(-0.5, 0, 0.5))),
               "splits must be in strictly increasing order and lie between lower and upper, got [0.5, 0.0, -0.5] along dimension 1",
               fixed = TRUE)
  expect_error(solve(prob, QuadDirect(), splits = list(c(-0.5, 0, 0.5), c(0, 0.5))),
               "3 values along each dimension required, got [0.0, 0.5] along dimension 2", fixed = TRUE)
  expect_error(solve(OptimizationProblem(rosen, c(0, 0), c(1, 100)), QuadDirect(), splits = s2),
               class = "IncompatibleOptimizerError")
  expect_true(requiresbounds(QuadDirect()))
  expect_false(allowscallback(QuadDirect()))
})
