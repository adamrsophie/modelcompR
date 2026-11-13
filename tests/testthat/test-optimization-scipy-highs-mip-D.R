# Group D of the HiGHS 1.12.0 MIP port: HighsPrimalHeuristics, the
# feasibility jump (mip/HighsFeasibilityJump.cpp, mip/feasibilityjump.hh) and
# the random number generators they draw from. References: the step traces
# of tools/fj-probe.py (g++ and libstdc++ on Linux x86-64) in
# fixtures/highs-fj-trace.txt and hand computation from the C++.

fj_fixture <- function() {
  readLines(test_path("fixtures", "highs-fj-trace.txt"))
}

# a function of another group the code under test calls: mocked in the
# namespace when it exists there, else defined in the global environment
# (the namespace's last resort) until heur_unmock removes it
heur_mock <- function(name, fn, env = parent.frame()) {
  ns <- asNamespace("modelcompR")
  if (exists(name, envir = ns, inherits = FALSE)) {
    if (is.function(fn)) do.call(local_mocked_bindings, c(stats::setNames(list(fn), name), list(.env = env)))
    return(invisible(NULL))
  }
  assign(name, fn, envir = globalenv())
  invisible(name)
}
heur_unmock <- function(names) {
  names <- unlist(names)
  if (length(names)) rm(list = names, envir = globalenv())
}
heur_src <- list(Branching = 0L, CentralRounding = 1L, FeasibilityPump = 2L, Heuristic = 3L, Shifting = 4L,
                 FeasibilityJump = 5L, SubMip = 6L, EmptyMip = 7L, RandomizedRounding = 8L, SolveLp = 9L,
                 EvaluateNode = 10L, Unbounded = 11L, UserSolution = 12L, HighsSolution = 13L, ZiRound = 14L)

test_that("std::mt19937 and the libstdc++ uniform_real_distribution match the C++", {
  g <- .hs_mt19937(5489)
  for (i in 1:9999) g()
  # the value the C++ standard mandates for the 10000th output
  expect_identical(g(), 4123659995)
  lines <- fj_fixture()
  mt <- as.numeric(strsplit(lines[startsWith(lines, "mt ")][1L], " ")[[1L]][-1L])
  g <- .hs_mt19937(0)
  expect_identical(vapply(1:24, function(i) g(), 0), mt)
  un <- as.numeric(strsplit(lines[startsWith(lines, "unif ")][1L], " ")[[1L]][-1L])
  g <- .hs_mt19937(0)
  expect_identical(vapply(1:24, function(i) .hs_uniform01(g), 0), un)
})

test_that("HighsRandom closedFraction, real and bit match the C++", {
  lines <- fj_fixture()
  ref <- function(key) strsplit(lines[startsWith(lines, paste0(key, " "))], " ")[[1L]][-1L]
  r <- .hs_rng(0)
  expect_identical(vapply(1:8, function(i) .hs_rng_closed_fraction(r), 0), as.numeric(ref("hrand_closed")))
  r <- .hs_rng(0)
  expect_identical(vapply(1:8, function(i) .hs_rng_real(r, 0.1, 0.9), 0), as.numeric(ref("hrand_real")))
  r <- .hs_rng(0)
  expect_identical(vapply(1:8, function(i) .hs_rng_bit(r), TRUE), ref("hrand_bit") == "1")
  r <- .hs_rng(0)
  mixed <- unlist(lapply(1:4, function(i) c(.hs_rng_real(r, 0.6, 0.6), .hs_rng_integer(r, 7))))
  expect_identical(mixed, as.numeric(ref("hrand_mixed")))
})

# parse one case of the trace fixture into the solver input
fj_case <- function(lines, name) {
  starts <- grep("^== fj ", lines)
  a <- starts[match(paste("== fj", name), lines[starts])]
  b <- c(starts, length(lines) + 1L)[match(a, starts) + 1L] - 1L
  blk <- lines[a:b]
  tok <- unlist(strsplit(sub("^in\\| ", "", blk[startsWith(blk, "in| ")]), " +"))
  p <- 0L
  nx <- function() { p <<- p + 1L; tok[p] }
  seed <- as.integer(nx()); eq <- as.numeric(nx()); vt <- as.numeric(nx()); nnz <- as.numeric(nx())
  n <- as.integer(nx())
  vars <- lapply(seq_len(n), function(j) list(int = nx() == "1", lb = as.numeric(nx()), ub = as.numeric(nx()),
                                             obj = as.numeric(nx())))
  m <- as.integer(nx())
  cons <- lapply(seq_len(m), function(i) {
    s <- as.integer(nx()); rhs <- as.numeric(nx()); k <- as.integer(nx())
    idx <- vapply(seq_len(k), function(t) as.integer(nx()), 0L)
    val <- vapply(seq_len(k), function(t) as.numeric(nx()), 0)
    list(sense = s, rhs = rhs, idx = idx, val = val)
  })
  init <- vapply(seq_len(n), function(j) as.numeric(nx()), 0)
  list(seed = seed, eq = eq, vt = vt, nnz = nnz, vars = vars, cons = cons, init = init,
       steps = blk[startsWith(blk, "s ")], end = blk[startsWith(blk, "end ")])
}

test_that("the feasibility jump solver follows the C++ step by step", {
  skip_on_cran()
  lines <- fj_fixture()
  for (name in c("t_milp1", "milp_mixed")) {
    cs <- fj_case(lines, name)
    fj <- .hmip_fj_new(cs$seed, cs$eq, cs$vt)
    for (v in cs$vars) fj$addVar(v$int, v$lb, v$ub, v$obj)
    for (c in cs$cons) fj$addConstraint(c$sense, c$rhs, c$idx, c$val, 0L)
    sol <- NULL; solobj <- NaN
    cb <- function(st) {
      if (!is.null(st$solution)) { sol <<- st$solution; solobj <<- st$solutionObjectiveValue }
      st$effortSinceLastImprovement > cs$nnz * 256 || st$totalEffort > cs$nnz * 1024
    }
    out <- character()
    tr <- function(step, v, val, eff, nv, ng, nb, obj) {
      if (step < 1000 || step %% 500 == 0)
        out[length(out) + 1L] <<- sprintf("s %d %d %a %.0f %d %d %.0f %a", step, v, val, eff, nv, ng, nb, obj)
    }
    fj$solve(cs$init, cb, tr)
    expect_identical(out, cs$steps, label = name)
    st <- fj$state()
    endref <- strsplit(cs$end, " ")[[1L]]
    expect_identical(st$totalEffort, as.numeric(endref[3L]))
    expect_identical(st$nBumps, as.numeric(endref[5L]))
    expect_identical(solobj, as.numeric(endref[9L]))
    expect_identical(sol, as.numeric(endref[-(1:10)]))
  }
})

# a HighsLp environment in the B1 layout from dense data
heur_lp <- function(cost, A, row_lower, row_upper, col_lower, col_upper, integrality) {
  lp <- new.env(parent = emptyenv())
  lp$num_col_ <- length(cost); lp$num_row_ <- nrow(A)
  lp$col_cost_ <- cost; lp$col_lower_ <- col_lower; lp$col_upper_ <- col_upper
  lp$row_lower_ <- row_lower; lp$row_upper_ <- row_upper
  lp$offset_ <- 0; lp$sense_ <- 1L; lp$integrality_ <- as.integer(integrality)
  start <- 0L; index <- integer(); value <- numeric()
  for (j in seq_len(ncol(A))) {
    nz <- which(A[, j] != 0)
    index <- c(index, nz - 1L); value <- c(value, A[nz, j]); start <- c(start, length(index))
  }
  lp$a_matrix_ <- list(format_ = "colwise", start_ = start, index_ = index, value_ = value)
  lp
}
heur_mipsolver <- function(lp) {
  ms <- new.env(parent = emptyenv())
  ms$model_ <- lp; ms$orig_model_ <- lp
  ms$options_mip_ <- list(random_seed = 0L, mip_root_presolve_only = FALSE)
  ms$submip <- FALSE
  md <- new.env(parent = emptyenv())
  md$mipsolver <- ms; md$feastol <- 1e-6; md$epsilon <- 1e-9; md$incumbent <- numeric()
  ms$mipdata_ <- md
  ms
}

test_that("feasibilityJump hands the solution of the C++ to trySolution", {
  skip_on_cran()
  # t_milp1 of tools/highs-probe.py with presolve off
  lp <- heur_lp(c(-5, -4, -3, -1), rbind(c(2, 3, 1, 1), c(4, 1, 2, 0.5), c(3, 4, 2, 1)), rep(-Inf, 3), c(5, 11, 8),
                rep(0, 4), c(10, 10, 10, 2.5), c(1, 1, 1, 0))
  ms <- heur_mipsolver(lp)
  got <- NULL
  mocked <- list(heur_mock(".hmip_HighsMipSolverData_trySolution", function(self, sol, src) {
    got <<- list(sol = sol, src = src); TRUE
  }), if (!exists(".hmip_src")) heur_mock(".hmip_src", heur_src))
  on.exit(heur_unmock(mocked), add = TRUE)
  expect_identical(.hmip_HighsMipSolverData_feasibilityJump(ms$mipdata_), "Notset")
  cs <- fj_case(fj_fixture(), "t_milp1")
  expect_identical(got$sol, as.numeric(strsplit(cs$end, " ")[[1L]][-(1:10)]))
  expect_identical(got$src, 5L)
  heur_unmock(mocked); mocked <- NULL
})

test_that("feasibilityJump reports integer columns with empty rounded bounds as infeasible", {
  # t_milp3: one integer column in [0.2, 0.8]
  lp <- heur_lp(1, matrix(0, 0, 1), numeric(), numeric(), 0.2, 0.8, 1)
  ms <- heur_mipsolver(lp)
  expect_identical(.hmip_HighsMipSolverData_feasibilityJump(ms$mipdata_), "Infeasible")
})

test_that("calcFixVal and determineTargetFixingRate follow the C++", {
  expect_identical(.hmip_calcFixVal(0.4, 2.3, 0), 3)
  expect_identical(.hmip_calcFixVal(-0.4, 2.3, 0), 2)
  expect_identical(.hmip_calcFixVal(0.1, 2.3, 1), 3)
  expect_identical(.hmip_calcFixVal(0.1, 2.3, -1), 2)
  expect_identical(.hmip_calcFixVal(0.1, 2.5, 0), 3)
  ms <- heur_mipsolver(heur_lp(1, matrix(0, 0, 1), numeric(), numeric(), 0, 1, 1))
  h <- .hmip_new_HighsPrimalHeuristics(ms)
  # no observations: real(0.6, 0.6), then the stream goes on as hrand_mixed
  expect_identical(.hmip_HighsPrimalHeuristics_determineTargetFixingRate(h), 0.6)
  expect_identical(.hs_rng_integer(h$randgen, 7), 3L)
  # infeasible at rate 0.5: low = high = 0.9 * 0.5
  h$numInfeasObservations <- 1L; h$infeasObservations <- 0.5
  expect_identical(.hmip_HighsPrimalHeuristics_determineTargetFixingRate(h), 0.9 * 0.5)
})

test_that("ziRound moves a fractional integer to the side its rows allow", {
  skip_if_not(exists(".cd_sub_cd") && exists(".hmip_fractionality"))
  # min -x  s.t.  x + y <= 2.5, x integer in [0, 3], y in [0, 1]; at the LP
  # point (1.5, 1) the row has no upper slack, so x can only move down by
  # 1.5 to the integer 0 (zi 0)
  lp <- heur_lp(c(-1, 0), matrix(c(1, 1), 1), -Inf, 2.5, c(0, 0), c(3, 1), c(1, 0))
  ms <- heur_mipsolver(lp)
  h <- .hmip_new_HighsPrimalHeuristics(ms)
  h$intcols <- 0L
  got <- NULL
  mocked <- list(heur_mock(".hmip_HighsMipSolverData_trySolution", function(self, sol, src) {
    got <<- list(sol = sol, src = src); TRUE
  }), if (!exists(".hmip_src")) heur_mock(".hmip_src", heur_src))
  on.exit(heur_unmock(mocked), add = TRUE)
  .hmip_HighsPrimalHeuristics_ziRound(h, c(1.5, 1))
  expect_identical(got$sol, c(0, 1))
  expect_identical(got$src, 14L)
  # an integral point is left alone
  got <- NULL
  .hmip_HighsPrimalHeuristics_ziRound(h, c(1, 1))
  expect_null(got)
  heur_unmock(mocked); mocked <- NULL
})

test_that("shifting rounds a fractional integer by its locks when no row is violated", {
  # x integer at 1.5 with one up lock and no down lock: the down rounding has
  # xi = uplocks = 1 > -1 and wins; the up rounding (xi 0) does not beat it
  lp <- heur_lp(c(1), matrix(1, 1), -Inf, 4, 0, 3, 1)
  ms <- heur_mipsolver(lp)
  md <- ms$mipdata_
  md$integer_cols <- 0L; md$uplocks <- 1; md$downlocks <- 0
  md$lp <- new.env(parent = emptyenv())
  md$lp$fractionalints <- list(first = 0L, second = 1.5)
  h <- .hmip_new_HighsPrimalHeuristics(ms)
  got <- NULL
  mocked <- list(
    heur_mock(".hmip_HighsMipSolverData_trySolution", function(self, sol, src) { got <<- list(sol = sol, src = src); TRUE }),
    heur_mock(".hmip_HighsMipSolverData_getInfeasibleRows", function(self, sol) list()),
    if (!exists(".hmip_HighsLpRelaxation_getFractionalIntegers"))
      heur_mock(".hmip_HighsLpRelaxation_getFractionalIntegers", function(lp) lp$fractionalints),
    if (!exists(".hmip_src")) heur_mock(".hmip_src", heur_src))
  on.exit(heur_unmock(mocked), add = TRUE)
  .hmip_HighsPrimalHeuristics_shifting(h, 1.5)
  expect_identical(got$sol, 1)
  expect_identical(got$src, 4L)
  heur_unmock(mocked); mocked <- NULL
})

test_that("flushStatistics moves the counters into mipdata", {
  ms <- heur_mipsolver(heur_lp(1, matrix(0, 0, 1), numeric(), numeric(), 0, 1, 1))
  md <- ms$mipdata_
  md$total_repair_lp <- 1; md$total_repair_lp_feasible <- 2; md$total_repair_lp_iterations <- 3
  md$heuristic_lp_iterations <- 10; md$total_lp_iterations <- 100
  h <- .hmip_new_HighsPrimalHeuristics(ms)
  h$total_repair_lp <- 4; h$total_repair_lp_feasible <- 5; h$total_repair_lp_iterations <- 6; h$lp_iterations <- 7
  .hmip_HighsPrimalHeuristics_flushStatistics(h)
  expect_identical(c(md$total_repair_lp, md$total_repair_lp_feasible, md$total_repair_lp_iterations,
                     md$heuristic_lp_iterations, md$total_lp_iterations), c(5, 7, 9, 17, 107))
  expect_identical(c(h$total_repair_lp, h$total_repair_lp_feasible, h$total_repair_lp_iterations, h$lp_iterations),
                   c(0, 0, 0, 0))
})

test_that("setupIntCols orders by lock score, clique score, hash and column", {
  skip_if_not(exists(".hmip_pdqsort") && exists(".hmip_hash"))
  ms <- heur_mipsolver(heur_lp(rep(1, 4), matrix(0, 0, 4), numeric(), numeric(), rep(0, 4), rep(1, 4), rep(1, 4)))
  md <- ms$mipdata_
  md$integer_cols <- 0:3
  md$uplocks <- c(1, 2, 1, 0); md$downlocks <- c(1, 2, 1, 5)
  md$cliquetable <- list(n1 = c(0, 0, 3, 0), n0 = c(0, 0, 1, 0))
  mocked <- list(heur_mock(".hmip_HighsCliqueTable_getNumImplications2", function(ct, col, val) {
    if (val == 1L) ct$n1[col + 1L] else ct$n0[col + 1L]
  }))
  on.exit(heur_unmock(mocked), add = TRUE)
  h <- .hmip_new_HighsPrimalHeuristics(ms)
  .hmip_HighsPrimalHeuristics_setupIntCols(h)
  # lock scores 1, 4, 1, 5e-6 (times feastol terms): col 1 first; cols 0 and
  # 2 tie on locks and col 2 has the larger clique score; col 3 last
  expect_identical(h$intcols, c(1L, 2L, 0L, 3L))
  heur_unmock(mocked); mocked <- NULL
})
