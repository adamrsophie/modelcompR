# IPX, the interior point solver of HiGHS 1.12.0 with its crossover, as
# SciPy 1.18.1 linprog(method = "highs-ipm", presolve = FALSE) runs it.
# Reference values from tools/highs-ipx-probe.py and, for the analytic centre
# modes, tools/highs-ipx-centring-probe.py (highspy 1.12.0), both run on
# Linux x86-64, the reference build of the HiGHS ports. The random LPs of the
# probes are in fixtures/highs-ipx-probe-data.txt. Iteration counts and bit
# patterns are compared on the reference platform only.

ipx_case <- function(c, Aub = NULL, bub = NULL, Aeq = NULL, beq = NULL, lo = 0, up = Inf) {
  n <- length(c)
  Aub <- if (is.null(Aub)) matrix(0, 0, n) else rbind(Aub)
  Aeq <- if (is.null(Aeq)) matrix(0, 0, n) else rbind(Aeq)
  list(c = c, A = rbind(Aub, Aeq), row_lower = c(rep(-Inf, nrow(Aub)), beq),
       row_upper = c(bub, beq), col_lower = rep_len(lo, n), col_upper = rep_len(up, n))
}

ipx_probe_data <- function() {
  l <- readLines(test_path("fixtures", "highs-ipx-probe-data.txt"))
  l <- l[!startsWith(l, "#")]
  v <- lapply(sub("^[^ ]+ c\\((.*)\\)$", "\\1", l),
              function(s) as.numeric(strsplit(s, ", ")[[1]]))
  stats::setNames(v, sub(" .*", "", l))
}

ipx_cases <- function() {
  km <- matrix(0, 5, 5)
  for (i in 1:5) {
    if (i > 1) for (j in 1:(i - 1)) km[i, j] <- 2^(i - j + 1)
    km[i, i] <- 1
  }
  At <- matrix(0, 7, 12)
  for (i in 1:3) At[i, 4 * (i - 1) + 1:4] <- 1
  for (j in 1:4) At[3 + j, seq(j, 12, 4)] <- 1
  d <- ipx_probe_data()
  list(
    equality = ipx_case(c(-1, -2, 0, 0, 0), Aeq = matrix(c(1, 1, 1, 0, 0, 1, -1, 0, 1, 0, 0, 1, 0, 0, 1), 3, byrow = TRUE),
                        beq = c(4, 1, 3)),
    inequality = ipx_case(c(-1, -2), matrix(c(1, 1, -1, 1, 1, -2), 3, byrow = TRUE), c(4, 2, 1)),
    dualized = ipx_case(c(-1, -1), matrix(c(1, 2, 2, 1, 1, 0, 0, 1, 1, 1), 5, byrow = TRUE), c(8, 8, 3.5, 3.5, 5)),
    dualized_boxed_free = ipx_case(c(1, -1), matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 2, 1), 5, byrow = TRUE),
                                   c(3, 1, 2, 1, 5), lo = c(-Inf, -1), up = c(Inf, 2)),
    boxed = ipx_case(c(-3, -1, -2), matrix(c(1, 1, 3, 2, 2, 5, 4, 1, 2), 3, byrow = TRUE), c(30, 24, 36),
                     lo = c(0, 1, 0), up = c(5, 4, 6)),
    free = ipx_case(c(1, 1, 2), matrix(c(-1, 1, 0, 1, 1, 1), 2, byrow = TRUE), c(1, 6), matrix(c(1, -1, 1), 1), 1,
                    lo = c(-Inf, -Inf, 0)),
    upper_fixed = ipx_case(c(1, -1, 1), matrix(c(1, 1, 1, 1, -1, 0), 2, byrow = TRUE), c(4, 1), matrix(c(0, 1, 1), 1), 2,
                           lo = c(-Inf, 0, 0.5), up = c(3, Inf, 0.5)),
    scaled = ipx_case(c(-1, -1e-3, -50), matrix(c(1e3, 2, 1e-2, 3e-3, 4e2, 5, 1, 1, 1e4), 3, byrow = TRUE),
                      c(1e4, 2e2, 3e5)),
    redundant = ipx_case(c(1, 2, 3), Aeq = matrix(c(1, 1, 1, 2, 2, 2, 1, -1, 0), 3, byrow = TRUE), beq = c(3, 6, 0.5)),
    infeasible = ipx_case(c(1, 1), matrix(c(1, 1, -1, -1), 2, byrow = TRUE), c(1, -2)),
    infeasible_eq = ipx_case(c(1, 1, 1), Aeq = matrix(c(1, 1, 0, 0, 1, 1, 1, 2, 1), 3, byrow = TRUE), beq = c(1, 1, 3)),
    unbounded = ipx_case(c(-1, -1), matrix(c(1, -1), 1), 1),
    unbounded_free = ipx_case(c(-1, 0, 1), Aeq = matrix(c(1, 1, 1), 1), beq = 1, lo = c(-Inf, -Inf, 0)),
    klee_minty = ipx_case(-2^(5 - 1:5), km, 5^(1:5)),
    degenerate = ipx_case(c(-1, -1, 0), matrix(c(1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 2), 4, byrow = TRUE), c(1, 1, 2, 2)),
    dense = ipx_case(d$c1, matrix(d$A1, 8), d$b1, matrix(d$Ae1, 3), d$be1, 0, 2),
    transport = ipx_case(c(8, 6, 10, 9, 9, 12, 13, 7, 14, 9, 16, 5), At[1:3, ], c(20, 30, 25), At[4:7, ],
                         c(10, 25, 15, 25)),
    sparse40 = ipx_case(d$c2, matrix(d$A2, 40), d$b2, lo = 0, up = 3),
    sparse_eq25 = ipx_case(d$c3, Aeq = matrix(d$Ae3, 25), beq = d$be3)
  )
}

# status, IPM iterations, crossover updates and the optimal objective of the
# probe
ipx_expected <- list(
  equality = list("kOptimal", 5, 0, -0x1.c000000000000p+2),
  inequality = list("kOptimal", 4, 0, -0x1.c000000000000p+2),
  dualized = list("kOptimal", 5, 0, -0x1.4000000000000p+2),
  dualized_boxed_free = list("kOptimal", 5, 0, -0x1.0000000000000p+0),
  boxed = list("kOptimal", 6, 0, -0x1.5666666666666p+4),
  free = list("kUnboundedOrInfeasible", 4, 0, NULL),
  upper_fixed = list("kUnboundedOrInfeasible", 5, 0, NULL),
  scaled = list("kOptimal", 7, 0, -0x1.797cc7e8d78c9p+10),
  redundant = list("kOptimal", 4, 1, 0x1.1000000000000p+2),
  infeasible = list("kInfeasible", 4, 0, NULL),
  infeasible_eq = list("kInfeasible", 5, 0, NULL),
  unbounded = list("kUnboundedOrInfeasible", 4, 0, NULL),
  unbounded_free = list("kUnboundedOrInfeasible", 5, 0, NULL),
  klee_minty = list("kOptimal", 6, 0, -0x1.86a0000000000p+11),
  degenerate = list("kOptimal", 5, 0, -0x1.0000000000000p+1),
  dense = list("kOptimal", 7, 0, NA),
  transport = list("kOptimal", 6, 0, NA),
  sparse40 = list("kOptimal", 12, 0, NA),
  sparse_eq25 = list("kOptimal", 8, 0, NA)
)

test_that("IPX reproduces the status, iterations and solutions of the probe", {
  cases <- ipx_cases()
  for (nm in names(ipx_expected)) {
    e <- ipx_expected[[nm]]
    lp <- cases[[nm]]
    r <- .highs_ipx_solve(lp)
    expect_identical(r$status, e[[1]], info = nm)
    if (!is.na(e[[2]])) expect_identical(r$ipm_iterations, as.integer(e[[2]]), info = nm)
    expect_identical(r$crossover_iterations, as.integer(e[[3]]), info = nm)
    if (r$status != "kOptimal") next
    if (!is.na(e[[4]])) expect_equal(r$objective, e[[4]], tolerance = ref_tol(0, 1e-9), info = nm)
    # primal feasibility, dual feasibility and complementarity of the basic
    # solution
    expect_lt(max(abs(drop(lp$A %*% r$x) - r$row_activity)), 1e-9)
    expect_true(all(r$x >= lp$col_lower - 1e-9 & r$x <= lp$col_upper + 1e-9), info = nm)
    expect_true(all(r$row_activity >= lp$row_lower - 1e-9 & r$row_activity <= lp$row_upper + 1e-9),
                info = nm)
    expect_lt(max(abs(lp$c - drop(crossprod(lp$A, r$row_dual)) - r$col_dual)), 1e-9)
    expect_equal(sum(r$col_basis == "kBasic") + sum(r$row_basis == "kBasic"), nrow(lp$A), info = nm)
  }
})

test_that("IPX returns the solutions and bases of the probe", {
  tol <- ref_tol(0, 1e-9)
  cases <- ipx_cases()
  r <- .highs_ipx_solve(cases$equality)
  expect_equal(r$x, c(1, 3, 0, 3, 0), tolerance = tol)
  expect_identical(r$col_basis, c("kBasic", "kBasic", "kLower", "kBasic", "kLower"))
  expect_identical(r$row_basis, c("kUpper", "kLower", "kUpper"))
  expect_equal(r$row_dual, c(-1, 0, -1), tolerance = tol)
  r <- .highs_ipx_solve(cases$dualized)
  expect_equal(r$x, c(3, 2), tolerance = tol)
  expect_identical(r$row_basis, c("kBasic", "kUpper", "kBasic", "kBasic", "kUpper"))
  r <- .highs_ipx_solve(cases$boxed)
  expect_equal(r$x, c(5, 4, 0x1.3333333333333p+0), tolerance = tol)
  expect_identical(r$col_basis, c("kUpper", "kUpper", "kBasic"))
  r <- .highs_ipx_solve(cases$scaled)
  expect_equal(r$x, c(0x1.3ffd8ae0184e9p+3, 0, 0x1.dffbe77498c79p+4), tolerance = tol)
  r <- .highs_ipx_solve(cases$redundant)
  expect_equal(r$x, c(0x1.c000000000000p+0, 0x1.4000000000000p+0, 0), tolerance = tol)
  expect_identical(r$col_basis, c("kBasic", "kBasic", "kLower"))
  expect_identical(r$row_basis, c("kLower", "kBasic", "kUpper"))
  r <- .highs_ipx_solve(cases$transport)
  expect_equal(r$x, c(0, 20, 0, 0, 10, 0, 15, 5, 0, 5, 0, 20), tolerance = tol)
  r <- .highs_ipx_solve(cases$dense)
  expect_equal(r$x, c(0, 0x1.b20de3bae3bfbp+0, 0, 0x1.a50de51db82a0p-1, 0x1.4f8dd4d589216p+0, 0, 0,
                      0x1.ab3baa1231347p+0, 0, 0x1.7c04c5457dd60p+0, 0x1.5a8abbaa9f563p+0,
                      0x1.0df15c16f1b56p+0), tolerance = tol)
})

test_that("IPX computes the analytic centre as HiGHS does", {
  cases <- ipx_cases()
  zero <- function(lp) { lp$c[] <- 0; lp }
  tol <- ref_tol(0, 1e-9)
  # the analytic centre of the MIP solver: zero costs, 200 iterations at
  # most, no crossover
  ac <- list(ipm_iteration_limit = 200L, run_crossover = "off")
  r <- .highs_ipx_solve(zero(cases$boxed), ac)
  expect_identical(r$status, "kOptimal")
  expect_identical(r$ipm_iterations, 7L)
  expect_equal(r$x, c(0x1.39ca587e15032p+1, 0x1.63482fdf8efd4p+0, 0x1.6d6ca2fdef1b2p-1), tolerance = tol)
  expect_null(r$col_basis)
  r <- .highs_ipx_solve(zero(cases$dense), ac)
  expect_identical(r$ipm_iterations, 6L)
  expect_equal(r$x, c(0x1.3716167a57d80p-1, 0x1.7ba18423ff39ep-1, 0x1.75e10a7dd15cep-1, 0x1.1f18b6791f4cdp-2,
                      0x1.712f57d3daa4cp+0, 0x1.a2a65ac383d00p-1, 0x1.9dbcf5d20431fp-2, 0x1.8e076d6c4a3c4p-1,
                      0x1.b0802c03b6904p-1, 0x1.350c96d8d52e1p-1, 0x1.5a1573af32115p-1, 0x1.d9683891ece28p-1),
               tolerance = tol)
  # run_centring: centring steps after the IPM, crossover off
  r <- .highs_ipx_solve(zero(cases$boxed), list(run_centring = TRUE))
  expect_identical(r$status, "kOptimal")
  expect_identical(r$ipm_iterations, 8L)
  expect_equal(r$x, c(0x1.30c6a43cd584ep+1, 0x1.77a6b223f5152p+0, 0x1.8692b2f236ae6p-1), tolerance = tol)
  # centring that does not reach the ratio tolerance leaves the IPM status
  # imprecise (HiGHS's lpKktCheck then declares it optimal)
  r <- .highs_ipx_solve(cases$transport, list(run_centring = TRUE))
  expect_identical(r$status, "kUnknown")
  expect_identical(r$ipm_iterations, 8L)
  expect_equal(r$x, c(0x1.74f42c2f074aap-37, 0x1.3ffffffff9548p+4, 0x1.7f3e507fc62fdp-35, 0x1.756e0a961fc3ap-35,
                      0x1.3ffffffff8b89p+3, 0x1.4e83c234d1fddp-35, 0x1.dffffffff1c7fp+3, 0x1.4000000020668p+2,
                      0x1.7e27b0ee93613p-35, 0x1.400000000ff08p+2, 0x1.068348f08a8fdp-34, 0x1.3ffffffff4e97p+4),
               tolerance = ref_tol(0, 1e-6))
})

test_that("IPX stops at the iteration limit", {
  cases <- ipx_cases()
  for (k in c(0L, 1L, 3L)) {
    r <- .highs_ipx_solve(cases$dense, list(ipm_iteration_limit = k))
    expect_identical(r$status, "kIterationLimit")
    expect_identical(r$ipm_iterations, k)
  }
  r <- .highs_ipx_solve(cases$sparse40, list(ipm_iteration_limit = 5L))
  expect_identical(r$status, "kIterationLimit")
  expect_identical(r$ipm_iterations, 5L)
  r <- .highs_ipx_solve(cases$sparse40, list(ipm_iteration_limit = 12L))
  expect_identical(r$status, "kOptimal")
})

test_that("IPX prints the iteration log of the probe", {
  out <- capture.output(r <- .highs_ipx_solve(ipx_cases()$equality, list(log = TRUE)))
  it <- grep("^ +[0-9]+\\*? +-?[0-9]", out, value = TRUE)
  expect_length(it, 6L)
  toks <- strsplit(trimws(it[1]), " +")[[1]]
  expect_identical(toks, c("0", "-3.69090909e+00", "-3.70094664e+00", "7.68e-01", "1.09e+00", "2.14e-03",
                           "0.00", "0.00", "0", "4", "0", "0"))
  expect_identical(nrow(r$log), 6L)
})

test_that("the BASICLU factorization solves, updates and finds dependent columns", {
  tocols <- function(B) lapply(seq_len(ncol(B)), function(j) {
    nz <- which(B[, j] != 0)
    list(i = nz - 1L, x = B[nz, j])
  })
  set.seed(11)
  m <- 30
  B <- matrix(stats::rnorm(m * m) * (stats::runif(m * m) < 0.15), m) + diag(4, m)
  lu <- .ipx_lu_new(m, 0.0625)
  f <- lu$factorize(tocols(B), FALSE)
  expect_identical(f$flag, 0L)
  rhs <- stats::rnorm(m)
  expect_lt(max(abs(B %*% lu$solve_dense(rhs, "N") - rhs)), 1e-12)
  expect_lt(max(abs(crossprod(B, lu$solve_dense(rhs, "T")) - rhs)), 1e-12)
  # the factors of basiclu_get_factors: B[rowperm, colperm] = (I + L) U
  g <- lu$get_factors()
  Lm <- diag(m); Um <- matrix(0, m, m)
  for (k in seq_len(m)) {
    Lm[g$L[[k]]$i, k] <- g$L[[k]]$x
    Um[g$U[[k]]$i, k] <- g$U[[k]]$x
  }
  expect_lt(max(abs(Lm %*% Um - B[g$rowperm, g$colperm])), 1e-12)
  # Forrest-Tomlin updates
  for (k in 1:20) {
    p <- sample(m, 1)
    col <- stats::rnorm(m) * (stats::runif(m) < 0.2)
    col[p] <- 3
    nz <- which(col != 0)
    ft <- lu$ftran(nz, col[nz])
    expect_lt(max(abs(B %*% ft - col)), 1e-10)
    bt <- lu$btran(p)
    ep <- numeric(m); ep[p] <- 1
    expect_lt(max(abs(crossprod(B, bt) - ep)), 1e-10)
    expect_identical(lu$update(ft[p]), 0L)
    B[, p] <- col
    expect_lt(max(abs(B %*% lu$solve_dense(rhs, "N") - rhs)), 1e-10)
    if (lu$need_fresh()) lu$factorize(tocols(B), FALSE)
  }
  # a dependent column is replaced by a slack column
  B[, 3] <- B[, 5] + B[, 7]
  lu <- .ipx_lu_new(m, 0.0625)
  f <- lu$factorize(tocols(B), FALSE)
  expect_identical(bitwAnd(f$flag, 2L), 2L)
  expect_identical(nrow(f$dependent), 1L)
})
