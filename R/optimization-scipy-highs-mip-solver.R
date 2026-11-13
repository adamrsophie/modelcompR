# The MIP solver driver of HiGHS 1.12.0: HighsMipSolver
# (mip/HighsMipSolver.cpp), HighsMipSolverData (mip/HighsMipSolverData.cpp),
# HighsObjectiveFunction (mip/HighsObjectiveFunction.cpp) and the MIP branch
# of Highs::run (Highs::callSolveMip and withoutSemiVariables of
# lp_data/Highs.cpp and lp_data/HighsLpUtils.cpp), as SciPy 1.18.1 vendors
# them. MIT license notice in R/optimization-scipy-highs-util.R.
#
# Every C++ class is an environment whose members keep their C++ names;
# stored index values (columns, rows, matrix starts) are 0-based as in C++
# and are read as v[i + 1L]. The other classes of the MIP solver (domain,
# clique table, LP relaxation, search, cuts, heuristics, symmetry, presolve)
# live in R/optimization-scipy-highs-mip-*.R and are called by the
# .hmip_<Class>_<method> names of their C++ methods.
#
# Deliberate differences: the solver is sequential, so the analytic center
# task that HiGHS spawns in startAnalyticCenterComputation runs when
# finishAnalyticCenterComputation synchronizes (its inputs do not change in
# between); callbacks are absent; time is read from a clock that stays at
# zero unless a finite time_limit is set, so logging and the primal dual
# integral behave as a fast run of the C++ does.

# HighsModelStatus values (without the k prefix) and solution sources
.hmip_src <- list(None = -1L, Branching = 0L, CentralRounding = 1L, FeasibilityPump = 2L,
                  Heuristic = 3L, Shifting = 4L, FeasibilityJump = 5L, SubMip = 6L,
                  EmptyMip = 7L, RandomizedRounding = 8L, SolveLp = 9L, EvaluateNode = 10L,
                  Unbounded = 11L, UserSolution = 12L, HighsSolution = 13L, ZiRound = 14L,
                  TrivialL = 15L, TrivialP = 16L, TrivialU = 17L, TrivialZ = 18L, Cleanup = 19L)
.hmip_src_code <- c("B", "C", "F", "H", "I", "J", "L", "P", "R", "S", "T", "U", "X", "Y",
                    "Z", "l", "p", "u", "z", " ")
.hmip_src_name <- c("Branching", "Central rounding", "Feasibility pump", "Heuristic",
                    "Shifting", "Feasibility jump", "Sub-MIP", "Empty MIP",
                    "Randomized rounding", "Solve LP", "Evaluate node", "Unbounded",
                    "User solution", "HiGHS solution", "ZI Round", "Trivial lower",
                    "Trivial point", "Trivial upper", "Trivial zero", "")
.hmip_kHighsIInf <- 2147483647

# HighsOptions defaults of 1.12.0 that the MIP solver reads
.hmip_default_options <- function() {
  list(mip_feasibility_tolerance = 1e-6, mip_rel_gap = 1e-4, mip_abs_gap = 1e-6,
       mip_max_nodes = .hmip_kHighsIInf, mip_max_stall_nodes = .hmip_kHighsIInf,
       mip_max_leaves = .hmip_kHighsIInf, mip_max_improving_sols = .hmip_kHighsIInf,
       mip_max_start_nodes = 500L, mip_lp_age_limit = 10L, mip_pool_age_limit = 30L,
       mip_pool_soft_limit = 10000L, mip_pscost_minreliable = 8L,
       mip_min_cliquetable_entries_for_parallelism = 100000L, mip_detect_symmetry = TRUE,
       mip_allow_restart = TRUE, mip_heuristic_effort = 0.05,
       mip_heuristic_run_feasibility_jump = TRUE, mip_heuristic_run_rins = TRUE,
       mip_heuristic_run_rens = TRUE, mip_heuristic_run_root_reduced_cost = TRUE,
       mip_heuristic_run_zi_round = FALSE, mip_heuristic_run_shifting = FALSE,
       mip_min_logging_interval = 5, mip_report_level = 1L, mip_root_presolve_only = FALSE,
       mip_lifting_for_probing = -1L, mip_improving_solution_save = FALSE,
       mip_ipm_solver = "choose", primal_feasibility_tolerance = 1e-7,
       dual_feasibility_tolerance = 1e-7, small_matrix_value = 1e-9, large_matrix_value = 1e15,
       time_limit = Inf, objective_bound = Inf, objective_target = -Inf, random_seed = 0L,
       presolve = "choose", presolve_reduction_limit = -1L,
       restart_presolve_reduction_limit = -1L, threads = 1L, timeless_log = FALSE, output_flag = FALSE,
       highs_debug_level = 0L)
}
.hmip_limit_set <- function(v) is.finite(v) && v != .hmip_kHighsIInf

# HighsLp as an environment (shared by reference like const HighsLp*)
.hmip_new_HighsLp <- function(col_cost, col_lower, col_upper, row_lower, row_upper,
                              start, index, value, integrality = integer(), offset = 0,
                              sense = 1L) {
  lp <- new.env(parent = emptyenv())
  lp$num_col_ <- length(col_cost); lp$num_row_ <- length(row_lower)
  lp$col_cost_ <- as.numeric(col_cost); lp$col_lower_ <- as.numeric(col_lower)
  lp$col_upper_ <- as.numeric(col_upper); lp$row_lower_ <- as.numeric(row_lower)
  lp$row_upper_ <- as.numeric(row_upper)
  lp$a_matrix_ <- list(format_ = "colwise", start_ = as.integer(start), index_ = as.integer(index),
                       value_ = as.numeric(value))
  lp$integrality_ <- as.integer(integrality)
  lp$offset_ <- offset; lp$sense_ <- as.integer(sense)
  lp
}
# HighsLp from a dense matrix (column-wise, zeros dropped)
.hmip_HighsLp_from_dense <- function(cost, A, row_lower, row_upper, col_lower, col_upper,
                                     integrality = integer(), offset = 0, sense = 1L) {
  n <- length(cost)
  start <- integer(n + 1L); index <- integer(); value <- numeric()
  for (j in seq_len(n)) {
    nz <- if (nrow(A)) which(A[, j] != 0) else integer()
    index <- c(index, nz - 1L); value <- c(value, A[nz, j])
    start[j + 1L] <- length(index)
  }
  .hmip_new_HighsLp(cost, col_lower, col_upper, row_lower, row_upper, start, index, value,
                    integrality, offset, sense)
}
.hmip_HighsLp_copy <- function(lp) {
  out <- new.env(parent = emptyenv())
  if (is.environment(lp)) lp <- as.list(lp, all.names = TRUE)
  for (k in names(lp)) assign(k, lp[[k]], envir = out)
  out
}
.hmip_HighsLp_dense <- function(lp) {
  A <- matrix(0, lp$num_row_, lp$num_col_)
  m <- lp$a_matrix_
  for (j in seq_len(lp$num_col_)) {
    k <- if (m$start_[j + 1L] > m$start_[j]) (m$start_[j] + 1L):m$start_[j + 1L] else integer()
    A[m$index_[k] + 1L, j] <- m$value_[k]
  }
  A
}
# calculateRowValuesQuad: HighsCDouble row sums in column order
.hmip_row_values_quad <- function(lp, x) {
  m <- lp$a_matrix_
  acc <- vector("list", lp$num_row_)
  for (i in seq_len(lp$num_row_)) acc[[i]] <- c(0, 0)
  for (j in seq_len(lp$num_col_)) {
    if (m$start_[j + 1L] > m$start_[j]) for (k in (m$start_[j] + 1L):m$start_[j + 1L]) {
      r <- m$index_[k] + 1L
      # row_value_quad[row] += col_value[col] * value: a rounded product
      acc[[r]] <- .cd_add(acc[[r]], x[j] * m$value_[k])
    }
  }
  vapply(acc, .cd_val, 0)
}
# highsSparseTranspose: the row-wise copy (0-based starts and indices)
.hmip_sparse_transpose <- function(nrow, ncol, start, index, value) {
  cnt <- integer(nrow)
  for (k in seq_along(index)) cnt[index[k] + 1L] <- cnt[index[k] + 1L] + 1L
  ARstart <- c(0L, cumsum(cnt))
  pos <- ARstart[seq_len(nrow)]
  ARindex <- integer(length(index)); ARvalue <- numeric(length(index))
  for (j in seq_len(ncol)) {
    if (start[j + 1L] > start[j]) for (k in (start[j] + 1L):start[j + 1L]) {
      r <- index[k] + 1L
      pos[r] <- pos[r] + 1L
      ARindex[pos[r]] <- j - 1L; ARvalue[pos[r]] <- value[k]
    }
  }
  list(ARstart_ = as.integer(ARstart), ARindex_ = ARindex, ARvalue_ = ARvalue)
}

# the solver clock (see the header): zero unless a time limit is set
.hmip_timer_read <- function(mipsolver) {
  if (is.finite(mipsolver$options_mip_$time_limit)) proc.time()[[3L]] - mipsolver$timer_start_ else 0
}
.hmip_log <- function(mipsolver, ...) {
  if (isTRUE(mipsolver$options_mip_$output_flag))
    mipsolver$log_lines <- c(mipsolver$log_lines, sprintf(...))
}

# ---------------------------------------------------------------------------
# HighsObjectiveFunction

.hmip_new_HighsObjectiveFunction <- function(mipsolver) {
  self <- new.env(parent = emptyenv())
  model <- mipsolver$model_
  self$model <- model
  n <- model$num_col_
  self$objectiveNonzeros <- which(model$col_cost_ != 0) - 1L
  self$colToPartition <- rep(-1L, n)
  self$cliquePartitionStart <- 0L
  self$objectiveVals <- numeric()
  if (!length(self$objectiveNonzeros)) {
    self$numBinary <- 0L; self$numIntegral <- 0L; self$objIntScale <- 1
    return(self)
  }
  vt <- model$integrality_
  p <- .hmip_std_partition(self$objectiveNonzeros, function(i) vt[i + 1L] != 0L)
  self$objectiveNonzeros <- p$v; self$numIntegral <- p$k
  if (self$numIntegral == 0L) {
    self$numBinary <- 0L
  } else {
    head <- self$objectiveNonzeros[seq_len(self$numIntegral)]
    q <- .hmip_std_partition(head, function(i)
      model$col_lower_[i + 1L] == 0 && model$col_upper_[i + 1L] == 1)
    self$objectiveNonzeros[seq_len(self$numIntegral)] <- q$v
    self$numBinary <- q$k
  }
  self$objectiveVals <- model$col_cost_[self$objectiveNonzeros + 1L]
  self$objIntScale <- 0
  self
}

.hmip_HighsObjectiveFunction_setupCliquePartition <- function(self, globaldom, cliqueTable) {
  if (self$numBinary <= 1L) return(invisible())
  model <- self$model
  bin <- self$objectiveNonzeros[seq_len(self$numBinary)]
  # CliqueVar(col, val) with val = cost < 0, encoded 2 * col + val
  clqvars <- 2L * bin + as.integer(model$col_cost_[bin + 1L] < 0)
  r <- .hmip_HighsCliqueTable_cliquePartitionObj(cliqueTable, model$col_cost_, clqvars)
  clqvars <- r$clqVars
  start <- r$partitionStart
  numPartitions <- length(start) - 1L
  if (numPartitions == self$numBinary) {
    self$cliquePartitionStart <- 0L
  } else {
    p <- 0L; k <- 0L
    for (i in seq_len(numPartitions) - 1L) {
      if (start[i + 2L] - start[i + 1L] == 1L) next
      start[p + 1L] <- k
      for (j in start[i + 1L]:(start[i + 2L] - 1L)) {
        col <- clqvars[j + 1L] %/% 2L
        self$colToPartition[col + 1L] <- k
        k <- k + 1L
      }
      p <- p + 1L
    }
    start[p + 1L] <- k
    self$cliquePartitionStart <- start[seq_len(p + 1L)]
    ctp <- self$colToPartition
    key <- function(i) {
      part <- ctp[i + 1L]
      list(part = if (part < 0L) part + 4294967296 else part, hash = .hmip_hash(i))
    }
    less <- function(a, b) {
      ka <- key(a); kb <- key(b)
      if (ka$part != kb$part) return(ka$part < kb$part)
      .hmip_u64_cmp(ka$hash, kb$hash) < 0L
    }
    self$objectiveNonzeros[seq_len(self$numBinary)] <-
      .hmip_pdqsort(self$objectiveNonzeros[seq_len(self$numBinary)], less)
    self$objectiveVals[seq_len(self$numBinary)] <-
      model$col_cost_[self$objectiveNonzeros[seq_len(self$numBinary)] + 1L]
  }
  invisible()
}

.hmip_HighsObjectiveFunction_checkIntegrality <- function(self, epsilon) {
  if (self$numIntegral == length(self$objectiveNonzeros)) {
    if (self$numIntegral) {
      self$objIntScale <- .hmip_integralScale(self$objectiveVals, epsilon, epsilon)
      if (self$objIntScale * 1e-14 > epsilon) self$objIntScale <- 0
    } else self$objIntScale <- 1
  }
  invisible()
}
.hmip_HighsObjectiveFunction_getObjectiveNonzeros <- function(self) self$objectiveNonzeros
.hmip_HighsObjectiveFunction_getObjectiveValuesPacked <- function(self) self$objectiveVals
.hmip_HighsObjectiveFunction_getNumBinariesInObjective <- function(self) self$numBinary
.hmip_HighsObjectiveFunction_getCliquePartitionStarts <- function(self) self$cliquePartitionStart
.hmip_HighsObjectiveFunction_getNumCliquePartitions <- function(self) length(self$cliquePartitionStart) - 1L
.hmip_HighsObjectiveFunction_getColCliquePartition <- function(self, col) self$colToPartition[col + 1L]
.hmip_HighsObjectiveFunction_integralScale <- function(self) self$objIntScale
.hmip_HighsObjectiveFunction_isIntegral <- function(self) self$objIntScale != 0
.hmip_HighsObjectiveFunction_isZero <- function(self) !length(self$objectiveNonzeros)

# ---------------------------------------------------------------------------
# HighsMipSolver

.hmip_new_HighsMipSolver <- function(options, lp, solution = NULL, submip = FALSE,
                                     submip_level = 0L) {
  ms <- new.env(parent = emptyenv())
  opts <- .hmip_default_options()
  opts[names(options)] <- options
  ms$callback_ <- NULL
  ms$options_mip_ <- opts
  ms$model_ <- lp
  ms$orig_model_ <- lp
  ms$modelstatus_ <- "Notset"
  ms$solution_ <- numeric()
  ms$solution_objective_ <- Inf
  ms$bound_violation_ <- 0; ms$integrality_violation_ <- 0; ms$row_violation_ <- 0
  ms$dual_bound_ <- 0; ms$primal_bound_ <- 0; ms$gap_ <- 0
  ms$node_count_ <- 0; ms$total_lp_iterations_ <- 0; ms$primal_dual_integral_ <- 0
  ms$saved_objective_and_solution_ <- list()
  ms$submip <- submip
  ms$submip_level <- submip_level
  ms$max_submip_level <- 0L
  ms$rootbasis <- NULL; ms$pscostinit <- NULL; ms$clqtableinit <- NULL; ms$implicinit <- NULL
  ms$mipdata_ <- NULL
  ms$termination_status_ <- "Notset"
  ms$terminator_ <- list(num_instance = 0L, my_instance = -1L, record = NULL)
  ms$timer_start_ <- proc.time()[[3L]]
  ms$log_lines <- character()
  if (!is.null(solution) && isTRUE(solution$value_valid)) {
    r <- .hmip_HighsMipSolver_solutionFeasible(ms, ms$orig_model_, solution$col_value,
                                               solution$row_value)
    ms$bound_violation_ <- r$bound_violation; ms$row_violation_ <- r$row_violation
    ms$integrality_violation_ <- r$integrality_violation
    ms$solution_objective_ <- .cd_val(r$obj)
    ms$solution_ <- solution$col_value
  }
  ms
}

.hmip_numCol <- function(mipsolver) mipsolver$model_$num_col_
.hmip_numRow <- function(mipsolver) mipsolver$model_$num_row_
.hmip_numNonzero <- function(mipsolver) length(mipsolver$model_$a_matrix_$index_)
.hmip_colCost <- function(mipsolver, col) mipsolver$model_$col_cost_[col + 1L]
.hmip_rowLower <- function(mipsolver, row) mipsolver$model_$row_lower_[row + 1L]
.hmip_rowUpper <- function(mipsolver, row) mipsolver$model_$row_upper_[row + 1L]
.hmip_variableType <- function(mipsolver, col) mipsolver$model_$integrality_[col + 1L]
.hmip_HighsMipSolver_terminate <- function(mipsolver) mipsolver$termination_status_ != "Notset"

# HighsMipSolver::initialiseTerminator(const HighsMipSolver&): a sub-MIP
# shares the terminator of its parent; without concurrent instances (always
# the case here) it is cleared
.hmip_HighsMipSolver_initialiseTerminator <- function(mipsolver, parent) {
  mipsolver$terminator_ <- list(num_instance = 0L, my_instance = -1L, record = NULL)
  if (!.hmip_HighsMipSolverData_terminatorActive(parent$mipdata_)) return(invisible())
  mipsolver$termination_status_ <- "Notset"
  mipsolver$terminator_ <- parent$terminator_
  invisible()
}

# HighsMipSolver::run
.hmip_HighsMipSolver_run <- function(mipsolver) {
  mipsolver$modelstatus_ <- "Notset"
  mipsolver$timer_start_ <- proc.time()[[3L]]
  opts <- mipsolver$options_mip_
  mipsolver$mipdata_ <- .hmip_new_HighsMipSolverData(mipsolver)
  md <- mipsolver$mipdata_
  .hmip_HighsMipSolverData_init(md)
  .hmip_HighsMipSolverData_runMipPresolve(md, opts$presolve_reduction_limit)
  if (mipsolver$modelstatus_ == "Notset" && .hmip_timer_read(mipsolver) >= opts$time_limit)
    mipsolver$modelstatus_ <- "TimeLimit"
  if (mipsolver$modelstatus_ != "Notset") {
    .hmip_log(mipsolver, "Presolve: %s", mipsolver$modelstatus_)
    if (mipsolver$modelstatus_ == "Optimal") {
      md$lower_bound <- 0
      md$upper_bound <- 0
      .hmip_HighsMipSolverData_transformNewIntegerFeasibleSolution(md, numeric())
      .hmip_HighsMipSolverData_saveReportMipSolution(md)
    }
    .hmip_HighsMipSolver_cleanupSolve(mipsolver)
    return(invisible())
  }
  .hmip_HighsMipSolverData_runSetup(md)
  repeat {  # restart:
    if (mipsolver$modelstatus_ == "Notset") {
      if (.hmip_HighsMipSolverData_checkLimits(md)) {
        .hmip_HighsMipSolver_cleanupSolve(mipsolver)
        return(invisible())
      }
      returned <- .hmip_HighsMipSolverData_trivialHeuristics(md)
      if (mipsolver$modelstatus_ == "Notset" && returned == "Infeasible") {
        mipsolver$modelstatus_ <- returned
        .hmip_HighsMipSolver_cleanupSolve(mipsolver)
        return(invisible())
      }
      if (isTRUE(opts$mip_heuristic_run_feasibility_jump)) {
        returned <- .hmip_HighsMipSolverData_feasibilityJump(md)
        if (mipsolver$modelstatus_ == "Notset" && returned == "Infeasible") {
          mipsolver$modelstatus_ <- returned
          .hmip_HighsMipSolver_cleanupSolve(mipsolver)
          return(invisible())
        }
      }
      .hmip_HighsMipSolverData_evaluateRootNode(md)
      if (.hmip_HighsMipSolver_terminate(mipsolver)) {
        mipsolver$modelstatus_ <- mipsolver$termination_status_
        .hmip_HighsMipSolver_cleanupSolve(mipsolver)
        return(invisible())
      }
      for (k in 1:5) .hmip_HighsCutPool_performAging(md$cutpool)
    }
    if (.hmip_HighsNodeQueue_empty(md$nodequeue) || .hmip_HighsMipSolverData_checkLimits(md)) {
      .hmip_HighsMipSolver_cleanupSolve(mipsolver)
      return(invisible())
    }
    if (!.hmip_search_loop(mipsolver)) break
  }
  .hmip_HighsMipSolver_cleanupSolve(mipsolver)
  invisible()
}

# lower_bound <- min(upper_bound, best queue bound) with the primal dual
# integral update that HighsMipSolver::run repeats after every change
.hmip_set_lower_bound <- function(md, value) {
  prev <- md$lower_bound
  md$lower_bound <- value
  if (!md$mipsolver$submip && md$lower_bound != prev)
    .hmip_HighsMipSolverData_updatePrimalDualIntegral(md, prev, md$lower_bound,
                                                      md$upper_bound, md$upper_bound)
}

# global bound changes found by propagation are pushed into the clique
# table, the implications and the local domain of the search
.hmip_apply_global_changes <- function(md, search) {
  if (!length(.hmip_HighsDomain_getChangedCols(md$domain))) return(invisible())
  .hmip_HighsCliqueTable_cleanupFixed(md$cliquetable, md$domain)
  for (col in .hmip_HighsDomain_getChangedCols(md$domain)) .hmip_HighsImplications_cleanupVarbounds(md$implications, col)
  .hmip_HighsDomain_setDomainChangeStack(md$domain, .hmip_domchg(numeric(), integer(), integer()))
  .hmip_HighsSearch_resetLocalDomain(search)
  .hmip_HighsDomain_clearChangedCols(md$domain)
  .hmip_HighsMipSolverData_removeFixedIndices(md)
}

# the branch and bound loop of HighsMipSolver::run after the root node;
# returns TRUE when the search asks for a restart (goto restart)
.hmip_search_loop <- function(mipsolver) {
  md <- mipsolver$mipdata_
  opts <- mipsolver$options_mip_
  search <- .hmip_new_HighsSearch(mipsolver, md$pseudocost)
  sepa <- .hmip_new_HighsSeparation(mipsolver)
  .hmip_HighsSearch_setLpRelaxation(search, md$lp)
  .hmip_HighsSeparation_setLpRelaxation(sepa, md$lp)
  .hmip_set_lower_bound(md, .hmip_HighsNodeQueue_getBestLowerBound(md$nodequeue))
  .hmip_HighsMipSolverData_printDisplayLine(md)
  .hmip_HighsSearch_installNode(search, .hmip_HighsNodeQueue_popBestBoundNode(md$nodequeue))
  numStallNodes <- 0; lastLbLeave <- 0; numQueueLeaves <- 0; numHugeTreeEstim <- 0L
  numNodesLastCheck <- md$num_nodes; nextCheck <- md$num_nodes
  treeweightLastCheck <- 0
  upperLimLastCheck <- md$upper_limit; lowerBoundLastCheck <- md$lower_bound
  while (.hmip_HighsSearch_hasNode(search)) {
    .hmip_HighsConflictPool_performAging(md$conflictPool)
    iterlimit <- 10 * max(.hmip_HighsLpRelaxation_getAvgSolveIters(md$lp), md$avgrootlpiters)
    iterlimit <- trunc(iterlimit)
    iterlimit <- max(10000, iterlimit, trunc((3 * md$firstrootlpiters) / 2))
    .hmip_HighsLpRelaxation_setIterationLimit(md$lp, iterlimit)
    plungestart <- md$num_nodes
    limit_reached <- FALSE
    considerHeuristics <- TRUE
    repeat {
      if (considerHeuristics && .hmip_HighsMipSolverData_moreHeuristicsAllowed(md)) {
        res <- .hmip_HighsSearch_evaluateNode(search)
        if (res == "SubOptimal") break
        if (.hmip_HighsSearch_currentNodePruned(search)) {
          md$num_leaves <- md$num_leaves + 1
          .hmip_HighsSearch_flushStatistics(search)
        } else {
          lpsol <- function() .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(md$lp))$col_value
          if (!length(md$incumbent))
            .hmip_HighsPrimalHeuristics_randomizedRounding(md$heuristics, lpsol())
          if (!length(md$incumbent)) {
            if (isTRUE(opts$mip_heuristic_run_rens)) .hmip_HighsPrimalHeuristics_RENS(md$heuristics, lpsol())
          } else {
            if (isTRUE(opts$mip_heuristic_run_rins)) .hmip_HighsPrimalHeuristics_RINS(md$heuristics, lpsol())
          }
          .hmip_HighsPrimalHeuristics_flushStatistics(md$heuristics)
          if (.hmip_HighsMipSolverData_terminatorTerminated(md)) return(FALSE)
        }
      }
      considerHeuristics <- FALSE
      if (.hmip_HighsDomain_infeasible(md$domain)) break
      if (!.hmip_HighsSearch_currentNodePruned(search)) {
        if (.hmip_HighsSearch_dive(search) == "SubOptimal") break
        md$num_leaves <- md$num_leaves + 1
        .hmip_HighsSearch_flushStatistics(search)
      }
      if (.hmip_HighsMipSolverData_checkLimits(md)) { limit_reached <- TRUE; break }
      if (md$num_nodes - plungestart >= 100) break
      if (!.hmip_HighsSearch_backtrackPlunge(search, md$nodequeue)) break
      if (.hmip_HighsConflictPool_getNumConflicts(md$conflictPool) > opts$mip_pool_soft_limit)
        .hmip_HighsConflictPool_performAging(md$conflictPool)
      .hmip_HighsSearch_flushStatistics(search)
      .hmip_HighsMipSolverData_printDisplayLine(md)
    }
    .hmip_HighsSearch_openNodesToQueue(search, md$nodequeue)
    .hmip_HighsSearch_flushStatistics(search)
    if (limit_reached) {
      .hmip_set_lower_bound(md, min(md$upper_bound, .hmip_HighsNodeQueue_getBestLowerBound(md$nodequeue)))
      .hmip_HighsMipSolverData_printDisplayLine(md)
      break
    }
    .hmip_HighsDomain_propagate(md$domain)
    md$pruned_treeweight <- .cd_add(md$pruned_treeweight,
      .hmip_HighsNodeQueue_pruneInfeasibleNodes(md$nodequeue, md$domain, md$feastol))
    if (.hmip_HighsDomain_infeasible(md$domain)) {
      .hmip_HighsNodeQueue_clear(md$nodequeue)
      md$pruned_treeweight <- c(1, 0)
      .hmip_set_lower_bound(md, min(Inf, md$upper_bound))
      .hmip_HighsMipSolverData_printDisplayLine(md)
      break
    }
    .hmip_set_lower_bound(md, min(md$upper_bound, .hmip_HighsNodeQueue_getBestLowerBound(md$nodequeue)))
    .hmip_HighsMipSolverData_printDisplayLine(md)
    if (.hmip_HighsNodeQueue_empty(md$nodequeue)) break
    .hmip_apply_global_changes(md, search)
    if (!mipsolver$submip && md$num_nodes >= nextCheck) {
      nTreeRestarts <- md$numRestarts - md$numRestartsRoot
      ptw <- .cd_val(md$pruned_treeweight)
      currNodeEstim <- numNodesLastCheck - md$num_nodes_before_run +
        (md$num_nodes - numNodesLastCheck) * .cd_val(.cd_add(-md$pruned_treeweight, 1)) /
        max(.cd_val(.cd_add(md$pruned_treeweight, -treeweightLastCheck)), md$epsilon)
      activeIntegerRatio <- 1 - .hmip_HighsMipSolverData_percentageInactiveIntegers(md) / 100
      activeIntegerRatio <- activeIntegerRatio * activeIntegerRatio
      gapReduction <- 1
      if (md$upper_limit != Inf) {
        oldGap <- upperLimLastCheck - lowerBoundLastCheck
        newGap <- md$upper_limit - md$lower_bound
        gapReduction <- oldGap / newGap
      }
      if (isTRUE(gapReduction < 1 + (0.05 / activeIntegerRatio)) &&
          isTRUE(currNodeEstim >= activeIntegerRatio * 20 * (md$num_nodes - md$num_nodes_before_run))) {
        nextCheck <- md$num_nodes + 100
        numHugeTreeEstim <- numHugeTreeEstim + 1L
      } else {
        numHugeTreeEstim <- 0L
        treeweightLastCheck <- ptw
        numNodesLastCheck <- md$num_nodes
        upperLimLastCheck <- md$upper_limit
        lowerBoundLastCheck <- md$lower_bound
      }
      doRestart <- FALSE
      if (isTRUE(opts$mip_allow_restart)) {
        minHugeTreeOffset <- (md$num_leaves - md$num_leaves_before_run) %/% 1000
        minHugeTreeEstim <- .hmip_nearestInteger(activeIntegerRatio * (10 + minHugeTreeOffset) *
                                                   1.5^nTreeRestarts)
        doRestart <- numHugeTreeEstim >= minHugeTreeEstim
      }
      if (doRestart) {
        .hmip_log(mipsolver, "\nRestarting search from the root node")
        .hmip_HighsMipSolverData_performRestart(md)
        return(TRUE)
      }
    }
    while (!.hmip_HighsNodeQueue_empty(md$nodequeue)) {
      if (numQueueLeaves - lastLbLeave >= 10) {
        .hmip_HighsSearch_installNode(search, .hmip_HighsNodeQueue_popBestBoundNode(md$nodequeue))
        lastLbLeave <- numQueueLeaves
      } else {
        bbStack <- .hmip_HighsNodeQueue_getBestBoundDomchgStackSize(md$nodequeue)
        bbLb <- .hmip_HighsNodeQueue_getBestLowerBound(md$nodequeue)
        nextNode <- .hmip_HighsNodeQueue_popBestNode(md$nodequeue)
        if (nextNode$lower_bound == bbLb && length(nextNode$domchgstack$column) == bbStack)
          lastLbLeave <- numQueueLeaves
        .hmip_HighsSearch_installNode(search, nextNode)
      }
      numQueueLeaves <- numQueueLeaves + 1
      if (.hmip_HighsSearch_getCurrentEstimate(search) >= md$upper_limit) {
        numStallNodes <- numStallNodes + 1
        if (.hmip_limit_set(opts$mip_max_stall_nodes) && numStallNodes >= opts$mip_max_stall_nodes) {
          limit_reached <- TRUE
          mipsolver$modelstatus_ <- "SolutionLimit"
          break
        }
      } else numStallNodes <- 0
      if (.hmip_HighsSearch_evaluateNode(search) == "SubOptimal")
        .hmip_HighsSearch_currentNodeToQueue(search, md$nodequeue)
      if (.hmip_HighsSearch_currentNodePruned(search)) {
        .hmip_HighsSearch_backtrack(search)
        md$num_leaves <- md$num_leaves + 1
        md$num_nodes <- md$num_nodes + 1
        .hmip_HighsSearch_flushStatistics(search)
        .hmip_HighsDomain_propagate(md$domain)
        md$pruned_treeweight <- .cd_add(md$pruned_treeweight,
          .hmip_HighsNodeQueue_pruneInfeasibleNodes(md$nodequeue, md$domain, md$feastol))
        if (.hmip_HighsDomain_infeasible(md$domain)) {
          .hmip_HighsNodeQueue_clear(md$nodequeue)
          md$pruned_treeweight <- c(1, 0)
          .hmip_set_lower_bound(md, min(Inf, md$upper_bound))
          break
        }
        if (.hmip_HighsMipSolverData_checkLimits(md)) { limit_reached <- TRUE; break }
        .hmip_set_lower_bound(md, min(md$upper_bound, .hmip_HighsNodeQueue_getBestLowerBound(md$nodequeue)))
        .hmip_HighsMipSolverData_printDisplayLine(md)
        .hmip_apply_global_changes(md, search)
        next
      }
      .hmip_HighsSeparation_separate(sepa, .hmip_HighsSearch_getLocalDomain(search))
      if (.hmip_HighsDomain_infeasible(md$domain)) {
        .hmip_HighsSearch_cutoffNode(search)
        .hmip_HighsSearch_openNodesToQueue(search, md$nodequeue)
        .hmip_HighsNodeQueue_clear(md$nodequeue)
        md$pruned_treeweight <- c(1, 0)
        .hmip_set_lower_bound(md, min(Inf, md$upper_bound))
        break
      }
      st <- .hmip_HighsLpRelaxation_getStatus(md$lp)
      if (st != "Error" && st != "NotSet") .hmip_HighsLpRelaxation_storeBasis(md$lp)
      basis <- .hmip_HighsLpRelaxation_getStoredBasis(md$lp)
      if (is.null(basis) || !.hmip_isBasisConsistent(.hmip_HighsLpRelaxation_getLp(md$lp), basis)) {
        b <- md$firstrootbasis
        nr <- .hmip_HighsLpRelaxation_numRows(md$lp)
        b$row_status <- c(b$row_status, rep("basic", max(0L, nr - length(b$row_status))))[seq_len(nr)]
        .hmip_HighsLpRelaxation_setStoredBasis(md$lp, b)
      }
      break
    }
    if (limit_reached) break
  }
  FALSE
}

# isBasisConsistent: sizes match and the number of basic variables is the
# number of rows
.hmip_isBasisConsistent <- function(lp, basis) {
  length(basis$col_status) == lp$num_col_ && length(basis$row_status) == lp$num_row_ &&
    sum(basis$col_status == "basic") + sum(basis$row_status == "basic") == lp$num_row_
}

# HighsMipSolver::cleanupSolve
.hmip_HighsMipSolver_cleanupSolve <- function(mipsolver) {
  md <- mipsolver$mipdata_
  opts <- mipsolver$options_mip_
  .hmip_HighsMipSolverData_printDisplayLine(md, .hmip_src$Cleanup)
  .hmip_HighsMipSolverData_updatePrimalDualIntegral(md, md$lower_bound, md$lower_bound,
                                                    md$upper_bound, md$upper_bound, FALSE)
  havesolution <- mipsolver$solution_objective_ != Inf
  tol <- opts$mip_feasibility_tolerance
  feasible <- havesolution && mipsolver$bound_violation_ <= tol &&
    mipsolver$integrality_violation_ <= tol && mipsolver$row_violation_ <= tol
  mipsolver$dual_bound_ <- md$lower_bound
  if (.hmip_HighsObjectiveFunction_isIntegral(md$objectiveFunction)) {
    s <- .hmip_HighsObjectiveFunction_integralScale(md$objectiveFunction)
    rounded <- ceiling(md$lower_bound * s - md$feastol) / s
    mipsolver$dual_bound_ <- max(mipsolver$dual_bound_, rounded)
  }
  mipsolver$dual_bound_ <- mipsolver$dual_bound_ + mipsolver$model_$offset_
  mipsolver$primal_bound_ <- md$upper_bound + mipsolver$model_$offset_
  mipsolver$node_count_ <- md$num_nodes
  mipsolver$total_lp_iterations_ <- md$total_lp_iterations
  mipsolver$dual_bound_ <- min(mipsolver$dual_bound_, mipsolver$primal_bound_)
  mipsolver$primal_dual_integral_ <- md$primal_dual_integral$value
  if (mipsolver$orig_model_$sense_ == -1L) {
    mipsolver$dual_bound_ <- -mipsolver$dual_bound_
    mipsolver$primal_bound_ <- -mipsolver$primal_bound_
  }
  if (mipsolver$modelstatus_ %in% c("Notset", "Infeasible"))
    mipsolver$modelstatus_ <- if (feasible && havesolution) "Optimal" else "Infeasible"
  pb <- mipsolver$primal_bound_; db <- mipsolver$dual_bound_
  mipsolver$gap_ <- if (pb == 0) { if (db == 0) 0 else Inf }
                    else if (pb != Inf) abs(pb - db) / abs(pb) else Inf
  .hmip_log(mipsolver, "\nSolving report")
  .hmip_log(mipsolver, "  Status            %s", mipsolver$modelstatus_)
  .hmip_log(mipsolver, "  Primal bound      %.12g", pb)
  .hmip_log(mipsolver, "  Dual bound        %.12g", db)
  .hmip_log(mipsolver, "  Nodes             %.0f", md$num_nodes)
  .hmip_log(mipsolver, "  LP iterations     %.0f", md$total_lp_iterations)
  invisible()
}

# HighsMipSolver::solutionFeasible
.hmip_HighsMipSolver_solutionFeasible <- function(mipsolver, lp, col_value, pass_row_value = NULL) {
  tol <- mipsolver$options_mip_$mip_feasibility_tolerance
  bv <- 0; rv <- 0; iv <- 0
  obj <- c(lp$offset_, 0)
  for (i in seq_len(lp$num_col_)) {
    value <- col_value[i]
    obj <- .cd_add(obj, lp$col_cost_[i] * value)
    if (length(lp$integrality_) && lp$integrality_[i] == 1L)
      iv <- max(.hmip_fractionality(value), iv)
    lower <- lp$col_lower_[i]; upper <- lp$col_upper_[i]
    if (value < lower - tol) bv <- max(bv, lower - value)
    else if (value > upper + tol) bv <- max(bv, value - upper)
  }
  if (lp$num_row_ > 0L) {
    row_value <- if (is.null(pass_row_value)) .hmip_row_values_quad(lp, col_value) else pass_row_value
    for (i in seq_len(lp$num_row_)) {
      value <- row_value[i]; lower <- lp$row_lower_[i]; upper <- lp$row_upper_[i]
      if (value < lower - tol) rv <- max(rv, lower - value)
      else if (value > upper + tol) rv <- max(rv, value - upper)
    }
  }
  list(feasible = bv <= tol && iv <= tol && rv <= tol, bound_violation = bv,
       row_violation = rv, integrality_violation = iv, obj = obj)
}

# HighsMipSolver::runMipPresolve (Highs::runPresolve for a MIP)
.hmip_HighsMipSolver_runMipPresolve <- function(mipsolver, presolve_reduction_limit) {
  mipsolver$mipdata_ <- .hmip_new_HighsMipSolverData(mipsolver)
  .hmip_HighsMipSolverData_init(mipsolver$mipdata_)
  .hmip_HighsMipSolverData_runMipPresolve(mipsolver$mipdata_, presolve_reduction_limit)
}

# ---------------------------------------------------------------------------
# HighsMipSolverData

.hmip_new_HighsMipSolverData <- function(mipsolver) {
  md <- new.env(parent = emptyenv())
  opts <- mipsolver$options_mip_
  ncol <- .hmip_numCol(mipsolver)
  md$mipsolver <- mipsolver
  # members in the C++ declaration (initialization) order
  md$cutpool <- .hmip_new_HighsCutPool(ncol, opts$mip_pool_age_limit, opts$mip_pool_soft_limit)
  md$conflictPool <- .hmip_new_HighsConflictPool(5L * opts$mip_pool_age_limit, opts$mip_pool_soft_limit)
  md$domain <- .hmip_new_HighsDomain(mipsolver)
  md$lp <- .hmip_new_HighsLpRelaxation(mipsolver)
  md$pseudocost <- .hmip_new_HighsPseudocost(NULL)
  md$cliquetable <- .hmip_new_HighsCliqueTable(ncol)
  md$implications <- .hmip_new_HighsImplications(mipsolver)
  md$heuristics <- .hmip_new_HighsPrimalHeuristics(mipsolver)
  md$redcostfixing <- .hmip_new_HighsRedcostFixing()
  md$objectiveFunction <- .hmip_new_HighsObjectiveFunction(mipsolver)
  md$postSolveStack <- .hmip_new_HighsPostsolveStack()
  md$presolve_status <- "NotSet"
  md$presolvedModel <- NULL
  md$cliquesExtracted <- FALSE; md$rowMatrixSet <- FALSE
  md$analyticCenterComputed <- FALSE; md$analyticCenterStatus <- "Notset"
  md$analyticCenterPending <- FALSE
  md$detectSymmetries <- FALSE
  md$numRestarts <- 0L; md$numRestartsRoot <- 0L
  md$numCliqueEntriesAfterPresolve <- 0L; md$numCliqueEntriesAfterFirstPresolve <- 0L
  md$ARstart_ <- integer(); md$ARindex_ <- integer(); md$ARvalue_ <- numeric()
  md$maxAbsRowCoef <- numeric(); md$rowintegral <- logical()
  md$uplocks <- integer(); md$downlocks <- integer()
  md$integer_cols <- integer(); md$implint_cols <- integer()
  md$integral_cols <- integer(); md$continuous_cols <- integer()
  md$symmetries <- NULL; md$globalOrbits <- NULL
  md$feastol <- 0; md$epsilon <- 0; md$heuristic_effort <- 0; md$dispfreq <- 0
  md$analyticCenter <- numeric(); md$firstlpsol <- numeric(); md$rootlpsol <- numeric()
  md$firstlpsolobj <- -Inf
  md$firstrootbasis <- list(valid = FALSE, alien = FALSE, useful = FALSE,
                            col_status = character(), row_status = character())
  md$rootlpsolobj <- -Inf
  md$numintegercols <- 0L; md$maxTreeSizeLog2 <- 0L
  md$pruned_treeweight <- c(0, 0)
  md$avgrootlpiters <- 0; md$disptime <- 0; md$last_disptime <- 0
  md$firstrootlpiters <- 0
  for (k in c("num_nodes", "num_leaves", "num_leaves_before_run", "num_nodes_before_run",
              "total_repair_lp", "total_repair_lp_feasible", "total_repair_lp_iterations",
              "total_lp_iterations", "heuristic_lp_iterations", "sepa_lp_iterations",
              "sb_lp_iterations", "total_lp_iterations_before_run",
              "heuristic_lp_iterations_before_run", "sepa_lp_iterations_before_run",
              "sb_lp_iterations_before_run", "num_disp_lines")) md[[k]] <- 0
  md$numImprovingSols <- 0L
  md$lower_bound <- -Inf; md$upper_bound <- Inf; md$upper_limit <- Inf; md$optimality_limit <- Inf
  md$incumbent <- numeric()
  md$nodequeue <- .hmip_new_HighsNodeQueue()
  md$primal_dual_integral <- list(value = -Inf, prev_lb = 0, prev_ub = 0, prev_gap = 0, prev_time = 0)
  .hmip_HighsDomain_addCutpool(md$domain, md$cutpool)
  .hmip_HighsDomain_addConflictPool(md$domain, md$conflictPool)
  md
}

.hmip_HighsMipSolverData_getRow <- function(self, row) {
  start <- self$ARstart_[row + 1L]; end <- self$ARstart_[row + 2L]
  k <- if (end > start) (start + 1L):end else integer()
  list(rowlen = end - start, rowinds = self$ARindex_[k], rowvals = self$ARvalue_[k])
}

.hmip_HighsMipSolverData_solutionSourceToString <- function(self, source, code = TRUE) {
  if (source == .hmip_src$None) return(if (code) " " else "None")
  if (code) .hmip_src_code[source + 1L] else .hmip_src_name[source + 1L]
}

# row activity in plain double arithmetic, as checkSolution and trySolution
# compute it
.hmip_row_activity <- function(self, row, solution) {
  start <- self$ARstart_[row + 1L]; end <- self$ARstart_[row + 2L]
  a <- 0
  if (end > start) for (j in (start + 1L):end)
    a <- a + solution[self$ARindex_[j] + 1L] * self$ARvalue_[j]
  a
}

.hmip_HighsMipSolverData_checkSolution <- function(self, solution) {
  ms <- self$mipsolver; model <- ms$model_
  for (i in seq_len(model$num_col_)) {
    if (solution[i] < model$col_lower_[i] - self$feastol) return(FALSE)
    if (solution[i] > model$col_upper_[i] + self$feastol) return(FALSE)
    if (model$integrality_[i] == 1L && .hmip_fractionality(solution[i]) > self$feastol) return(FALSE)
  }
  for (i in seq_len(model$num_row_) - 1L) {
    a <- .hmip_row_activity(self, i, solution)
    if (a > model$row_upper_[i + 1L] + self$feastol) return(FALSE)
    if (a < model$row_lower_[i + 1L] - self$feastol) return(FALSE)
  }
  TRUE
}

.hmip_HighsMipSolverData_getInfeasibleRows <- function(self, solution) {
  model <- self$mipsolver$model_
  out <- list()
  for (i in seq_len(model$num_row_) - 1L) {
    start <- self$ARstart_[i + 1L]; end <- self$ARstart_[i + 2L]
    q <- c(0, 0)
    if (end > start) for (j in (start + 1L):end)
      q <- .cd_iadd(q, .cd_two_product(solution[self$ARindex_[j] + 1L], self$ARvalue_[j]))
    a <- .cd_val(q)
    if (a > model$row_upper_[i + 1L] + self$feastol)
      out[[length(out) + 1L]] <- list(i, 1L, abs(a - model$row_upper_[i + 1L]))
    if (a < model$row_lower_[i + 1L] - self$feastol)
      out[[length(out) + 1L]] <- list(i, -1L, abs(model$row_lower_[i + 1L] - a))
  }
  out
}

.hmip_HighsMipSolverData_trySolution <- function(self, solution, solution_source = -1L) {
  model <- self$mipsolver$model_
  if (length(solution) != model$num_col_) return(FALSE)
  obj <- c(0, 0)
  for (i in seq_len(model$num_col_)) {
    if (solution[i] < model$col_lower_[i] - self$feastol) return(FALSE)
    if (solution[i] > model$col_upper_[i] + self$feastol) return(FALSE)
    if (model$integrality_[i] == 1L && .hmip_fractionality(solution[i]) > self$feastol) return(FALSE)
    obj <- .cd_add(obj, model$col_cost_[i] * solution[i])
  }
  for (i in seq_len(model$num_row_) - 1L) {
    a <- .hmip_row_activity(self, i, solution)
    if (a > model$row_upper_[i + 1L] + self$feastol) return(FALSE)
    if (a < model$row_lower_[i + 1L] - self$feastol) return(FALSE)
  }
  .hmip_HighsMipSolverData_addIncumbent(self, solution, .cd_val(obj), solution_source)
}

.hmip_HighsMipSolverData_solutionRowFeasible <- function(self, solution) {
  model <- self$mipsolver$model_
  for (i in seq_len(model$num_row_) - 1L) {
    start <- self$ARstart_[i + 1L]; end <- self$ARstart_[i + 2L]
    q <- c(0, 0)
    if (end > start) for (j in (start + 1L):end)
      q <- .cd_add(q, solution[self$ARindex_[j] + 1L] * self$ARvalue_[j])
    a <- .cd_val(q)
    if (a > model$row_upper_[i + 1L] + self$feastol) return(FALSE)
    if (a < model$row_lower_[i + 1L] - self$feastol) return(FALSE)
  }
  TRUE
}

.hmip_HighsMipSolverData_trivialHeuristics <- function(self) {
  if (length(self$continuous_cols) > 0L) return("Notset")
  ms <- self$mipsolver; model <- ms$model_
  sources <- c(.hmip_src$TrivialZ, .hmip_src$TrivialL, .hmip_src$TrivialU, .hmip_src$TrivialP)
  col_lower <- model$col_lower_; col_upper <- model$col_upper_
  lower_non_positive <- TRUE; lower_zero <- TRUE; lower_finite <- TRUE; upper_finite <- TRUE
  for (k in seq_len(self$numintegercols)) {
    c1 <- self$integer_cols[k] + 1L
    col_lower[c1] <- ceiling(col_lower[c1]); col_upper[c1] <- floor(col_upper[c1])
    legal <- col_lower[c1] <= col_upper[c1] && col_lower[c1] < Inf && col_upper[c1] > -Inf &&
      !is.nan(col_lower[c1]) && !is.nan(col_upper[c1])
    if (!legal) {
      .hmip_log(ms, "HighsMipSolverData::trivialHeuristics() has detected infeasible/illegal bounds [%g, %g] for column %d: MIP is infeasible",
                col_lower[c1], col_upper[c1], c1 - 1L)
      return("Infeasible")
    }
    if (col_lower[c1] > 0) lower_non_positive <- FALSE
    if (col_lower[c1] != 0) lower_zero <- FALSE
    if (col_lower[c1] <= -Inf) lower_finite <- FALSE
    if (col_upper[c1] >= Inf) upper_finite <- FALSE
    if (!(lower_non_positive || lower_zero || upper_finite)) break
  }
  boxed <- lower_finite && upper_finite
  tol <- ms$options_mip_$mip_feasibility_tolerance
  solution <- numeric(model$num_col_)
  m <- model$a_matrix_
  for (h in 0:3) {
    if (h == 0L) {
      if (!lower_non_positive) next
      if (any(model$row_lower_ > tol | model$row_upper_ < -tol)) next
      solution <- numeric(model$num_col_)
    } else if (h == 1L) {
      if (lower_zero) next
      if (!.hmip_HighsMipSolverData_solutionRowFeasible(self, col_lower)) next
      solution <- col_lower
    } else if (h == 2L) {
      if (!upper_finite) next
      if (!.hmip_HighsMipSolverData_solutionRowFeasible(self, col_upper)) next
      solution <- col_upper
    } else {
      if (!boxed) next
      for (k in seq_len(self$numintegercols)) {
        c1 <- self$integer_cols[k] + 1L
        els <- if (m$start_[c1 + 1L] > m$start_[c1]) (m$start_[c1] + 1L):m$start_[c1 + 1L] else integer()
        npos <- sum(m$value_[els] > 0); nneg <- length(els) - npos
        solution[c1] <- if (npos > nneg) col_lower[c1] else col_upper[c1]
      }
      if (!.hmip_HighsMipSolverData_solutionRowFeasible(self, solution)) next
    }
    obj <- c(0, 0)
    for (j in seq_len(model$num_col_)) obj <- .cd_add(obj, model$col_cost_[j] * solution[j])
    .hmip_HighsMipSolverData_addIncumbent(self, solution, .cd_val(obj), sources[h + 1L])
  }
  "Notset"
}

# startAnalyticCenterComputation: the task runs at the synchronization
.hmip_HighsMipSolverData_startAnalyticCenterComputation <- function(self) {
  self$analyticCenterPending <- TRUE
  invisible()
}

# the task body: Highs ipm (.hs_highs_run_ipm) with presolve off, IPX, ipm_iteration_limit 200,
# run_crossover off, on the model with zero costs and no integrality
.hmip_analytic_center_task <- function(self) {
  model <- self$mipsolver$model_
  n <- model$num_col_
  r <- .hs_highs_run_ipm(numeric(n), .hmip_HighsLp_dense(model), model$row_lower_, model$row_upper_,
                         model$col_lower_, model$col_upper_,
                         list(presolve = FALSE, run_crossover = "off", iteration_limit = 200L))
  sol <- r$col_value
  if (length(sol) != n) return(invisible())
  self$analyticCenterStatus <- r$status
  self$analyticCenter <- sol
  invisible()
}

.hmip_HighsMipSolverData_finishAnalyticCenterComputation <- function(self) {
  if (isTRUE(self$analyticCenterPending)) {
    self$analyticCenterPending <- FALSE
    .hmip_analytic_center_task(self)
  }
  self$analyticCenterComputed <- TRUE
  ms <- self$mipsolver; model <- ms$model_; dom <- self$domain
  if (self$analyticCenterStatus == "Optimal") {
    nfixed <- 0L; nintfixed <- 0L
    for (i in seq_len(model$num_col_)) {
      boundRange <- dom$col_upper_[i] - dom$col_lower_[i]
      if (boundRange == 0) next
      tolerance <- self$feastol * min(boundRange, 1)
      if (self$analyticCenter[i] <= model$col_lower_[i] + tolerance) {
        .hmip_HighsDomain_changeBound3(dom, .hmip_kUpper, i - 1L, model$col_lower_[i],
                                       .hmip_reason_unspecified())
        if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
        nfixed <- nfixed + 1L
        if (model$integrality_[i] == 1L) nintfixed <- nintfixed + 1L
      } else if (self$analyticCenter[i] >= model$col_upper_[i] - tolerance) {
        .hmip_HighsDomain_changeBound3(dom, .hmip_kLower, i - 1L, model$col_upper_[i],
                                       .hmip_reason_unspecified())
        if (.hmip_HighsDomain_infeasible(dom)) return(invisible())
        nfixed <- nfixed + 1L
        if (model$integrality_[i] == 1L) nintfixed <- nintfixed + 1L
      }
    }
    .hmip_HighsDomain_propagate(dom)
  }
  invisible()
}

# symmetry detection runs inline (group F: .hmip_symmetry_start/finish)
.hmip_HighsMipSolverData_startSymmetryDetection <- function(self) {
  self$symData <- .hmip_symmetry_start(self)
  invisible()
}

.hmip_HighsMipSolverData_finishSymmetryDetection <- function(self) {
  msg <- .hmip_symmetry_finish(self, self$symData)
  .hmip_log(self$mipsolver, "\nSymmetry detection completed in%s\n",
            if (isTRUE(self$mipsolver$options_mip_$timeless_log)) "" else " 0.0s")
  .hmip_log(self$mipsolver, "%s\n", msg)
  self$symData <- NULL
  invisible()
}

.hmip_HighsMipSolverData_limitsToGap <- function(self, use_lower_bound, use_upper_bound) {
  offset <- self$mipsolver$model_$offset_
  lb <- use_lower_bound + offset
  if (abs(lb) <= self$epsilon) lb <- 0
  ub <- Inf; gap <- Inf
  if (use_upper_bound != Inf) {
    ub <- use_upper_bound + offset
    if (abs(ub) <= self$epsilon) ub <- 0
    lb <- min(ub, lb)
    gap <- if (ub == 0) { if (lb == 0) 0 else Inf } else (ub - lb) / abs(ub)
  }
  list(gap = gap, lb = lb, ub = ub)
}

.hmip_HighsMipSolverData_computeNewUpperLimit <- function(self, ub, mip_abs_gap, mip_rel_gap) {
  of <- self$objectiveFunction
  offset <- self$mipsolver$model_$offset_
  if (.hmip_HighsObjectiveFunction_isIntegral(of)) {
    s <- .hmip_HighsObjectiveFunction_integralScale(of)
    nul <- floor(s * ub - 0.5) / s
    if (mip_rel_gap != 0)
      nul <- min(nul, ub - ceiling(mip_rel_gap * abs(ub + offset) * s - self$epsilon) / s)
    if (mip_abs_gap != 0)
      nul <- min(nul, ub - ceiling(mip_abs_gap * s - self$epsilon) / s)
    nul <- nul + self$feastol
  } else {
    nul <- min(ub - self$feastol, .hmip_nextafter_down(ub))
    if (mip_rel_gap != 0) nul <- min(nul, ub - mip_rel_gap * abs(ub + offset))
    if (mip_abs_gap != 0) nul <- min(nul, ub - mip_abs_gap)
  }
  nul
}

.hmip_HighsMipSolverData_moreHeuristicsAllowed <- function(self) {
  ptw <- .cd_val(self$pruned_treeweight)
  if (self$mipsolver$submip) {
    return(self$heuristic_lp_iterations < self$total_lp_iterations * self$heuristic_effort)
  } else if (ptw < 1e-3 && self$num_leaves - self$num_leaves_before_run < 10 &&
             self$num_nodes - self$num_nodes_before_run < 1000) {
    if (self$heuristic_lp_iterations < self$total_lp_iterations * self$heuristic_effort + 10000)
      return(TRUE)
  } else if (self$heuristic_lp_iterations <
             100000 + ((self$total_lp_iterations - self$heuristic_lp_iterations -
                        self$sb_lp_iterations) %/% 2)) {
    heur_curr <- self$heuristic_lp_iterations - self$heuristic_lp_iterations_before_run
    sb_curr <- self$sb_lp_iterations - self$sb_lp_iterations_before_run
    node_curr <- self$total_lp_iterations - self$total_lp_iterations_before_run - heur_curr - sb_curr
    estim <- self$heuristic_lp_iterations /
      ((self$total_lp_iterations - node_curr) + node_curr / max(0.01, ptw))
    if (estim < max(0.3 / 0.8, min(ptw, 0.8) / 0.8) * self$heuristic_effort) return(TRUE)
  }
  FALSE
}

.hmip_HighsMipSolverData_removeFixedIndices <- function(self) {
  keep <- function(cols) cols[!vapply(cols, function(c) .hmip_HighsDomain_isFixed(self$domain, c), TRUE)]
  self$integral_cols <- keep(self$integral_cols)
  self$integer_cols <- keep(self$integer_cols)
  self$implint_cols <- keep(self$implint_cols)
  self$continuous_cols <- keep(self$continuous_cols)
  invisible()
}

.hmip_HighsMipSolverData_init <- function(self) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  .hmip_HighsPostsolveStack_initializeIndexMaps(self$postSolveStack, ms$model_$num_row_, ms$model_$num_col_)
  ms$orig_model_ <- ms$model_
  self$feastol <- opts$mip_feasibility_tolerance
  self$epsilon <- opts$small_matrix_value
  if (!is.null(ms$clqtableinit)) .hmip_HighsCliqueTable_buildFrom(self$cliquetable, ms$orig_model_, ms$clqtableinit)
  .hmip_HighsCliqueTable_setMinEntriesForParallelism(self$cliquetable, .hmip_kHighsIInf)
  if (!is.null(ms$implicinit)) .hmip_HighsImplications_buildFrom(self$implications, ms$implicinit)
  self$heuristic_effort <- opts$mip_heuristic_effort
  self$detectSymmetries <- isTRUE(opts$mip_detect_symmetry)
  self$firstlpsolobj <- -Inf; self$rootlpsolobj <- -Inf
  self$analyticCenterComputed <- FALSE; self$analyticCenterStatus <- "Notset"
  self$maxTreeSizeLog2 <- 0L; self$numRestarts <- 0L; self$numRestartsRoot <- 0L
  self$numImprovingSols <- 0L; self$pruned_treeweight <- c(0, 0); self$avgrootlpiters <- 0
  for (k in c("num_nodes", "num_nodes_before_run", "num_leaves", "num_leaves_before_run",
              "total_repair_lp", "total_repair_lp_feasible", "total_repair_lp_iterations",
              "total_lp_iterations", "heuristic_lp_iterations", "sepa_lp_iterations",
              "sb_lp_iterations", "total_lp_iterations_before_run",
              "heuristic_lp_iterations_before_run", "sepa_lp_iterations_before_run",
              "sb_lp_iterations_before_run", "num_disp_lines")) self[[k]] <- 0
  self$numCliqueEntriesAfterPresolve <- 0L; self$numCliqueEntriesAfterFirstPresolve <- 0L
  self$cliquesExtracted <- FALSE; self$rowMatrixSet <- FALSE
  self$lower_bound <- -Inf; self$upper_bound <- Inf
  self$upper_limit <- opts$objective_bound; self$optimality_limit <- opts$objective_bound
  self$primal_dual_integral$value <- -Inf
  self$dispfreq <- if (opts$mip_report_level == 0) 0 else if (opts$mip_report_level == 1) 2000 else 100
  invisible()
}

# runMipPresolve: HPresolve on the model (group E); on return
# mipsolver$model_ is the presolved model held in presolvedModel
.hmip_HighsMipSolverData_runMipPresolve <- function(self, presolve_reduction_limit) {
  ms <- self$mipsolver
  presolve <- .hmip_new_HPresolve()
  if (!.hmip_HPresolve_okSetInput(presolve, ms, presolve_reduction_limit)) {
    ms$modelstatus_ <- "MemoryLimit"
    self$presolve_status <- "OutOfMemory"
  } else {
    ms$modelstatus_ <- .hmip_HPresolve_run(presolve, self$postSolveStack)
    self$presolve_status <- .hmip_HPresolve_getPresolveStatus(presolve)
  }
  invisible()
}

.hmip_HighsMipSolverData_runSetup <- function(self) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  model <- ms$model_
  .hmip_HighsLpRelaxation_setSolvedFirstLp(self$lp, FALSE)
  self$last_disptime <- -Inf
  self$disptime <- 0
  self$upper_limit <- self$upper_limit - model$offset_
  self$optimality_limit <- self$optimality_limit - model$offset_
  self$lower_bound <- self$lower_bound - model$offset_
  self$upper_bound <- self$upper_bound - model$offset_
  if (ms$solution_objective_ != Inf) {
    self$incumbent <- .hmip_HighsPostsolveStack_getReducedPrimalSolution(self$postSolveStack, ms$solution_)
    solobj <- ms$solution_objective_ * ms$orig_model_$sense_ - model$offset_
    tol <- opts$mip_feasibility_tolerance
    feasible <- ms$bound_violation_ <= tol && ms$integrality_violation_ <= tol && ms$row_violation_ <= tol
    if (self$numRestarts == 0L)
      .hmip_log(ms, "\nMIP start solution is %s, objective value is %.12g",
                if (feasible) "feasible" else "infeasible", ms$solution_objective_)
    if (feasible && solobj < self$upper_bound) {
      prev <- self$upper_bound
      self$upper_bound <- solobj
      if (!ms$submip && self$upper_bound != prev)
        .hmip_HighsMipSolverData_updatePrimalDualIntegral(self, self$lower_bound, self$lower_bound,
                                                          prev, self$upper_bound)
      nul <- .hmip_HighsMipSolverData_computeNewUpperLimit(self, solobj, 0, 0)
      .hmip_HighsMipSolverData_saveReportMipSolution(self, nul)
      if (nul < self$upper_limit) {
        self$upper_limit <- nul
        self$optimality_limit <- .hmip_HighsMipSolverData_computeNewUpperLimit(
          self, solobj, opts$mip_abs_gap, opts$mip_rel_gap)
        .hmip_HighsNodeQueue_setOptimalityLimit(self$nodequeue, self$optimality_limit)
      }
    }
  }
  if (model$num_col_ == 0L)
    .hmip_HighsMipSolverData_addIncumbent(self, numeric(), 0, .hmip_src$EmptyMip)
  self$redcostfixing <- .hmip_new_HighsRedcostFixing()
  .hmip_HighsPseudocost_init(self$pseudocost, ms)
  .hmip_HighsNodeQueue_setNumCol(self$nodequeue, model$num_col_)
  .hmip_HighsNodeQueue_setOptimalityLimit(self$nodequeue, self$optimality_limit)
  self$continuous_cols <- integer(); self$integer_cols <- integer()
  self$implint_cols <- integer(); self$integral_cols <- integer()
  m <- model$a_matrix_
  tr <- .hmip_sparse_transpose(model$num_row_, model$num_col_, m$start_, m$index_, m$value_)
  self$ARstart_ <- tr$ARstart_; self$ARindex_ <- tr$ARindex_; self$ARvalue_ <- tr$ARvalue_
  self$rowMatrixSet <- TRUE
  up <- integer(model$num_col_); down <- integer(model$num_col_)
  for (i in seq_len(model$num_col_)) {
    if (m$start_[i + 1L] > m$start_[i]) for (j in (m$start_[i] + 1L):m$start_[i + 1L]) {
      row <- m$index_[j] + 1L
      if (model$row_lower_[row] != -Inf) {
        if (m$value_[j] < 0) up[i] <- up[i] + 1L else down[i] <- down[i] + 1L
      }
      if (model$row_upper_[row] != Inf) {
        if (m$value_[j] < 0) down[i] <- down[i] + 1L else up[i] <- up[i] + 1L
      }
    }
  }
  self$uplocks <- up; self$downlocks <- down
  self$rowintegral <- logical(model$num_row_)
  self$maxAbsRowCoef <- numeric(model$num_row_)
  pm <- self$presolvedModel
  for (i in seq_len(model$num_row_)) {
    maxabsval <- 0
    integral <- TRUE
    start <- self$ARstart_[i]; end <- self$ARstart_[i + 1L]
    if (end > start) for (j in (start + 1L):end) {
      integral <- integral && model$integrality_[self$ARindex_[j] + 1L] != 0L &&
        .hmip_fractionality(self$ARvalue_[j]) <= self$epsilon
      maxabsval <- max(maxabsval, abs(self$ARvalue_[j]))
    }
    if (integral && !is.null(pm)) {
      if (pm$row_lower_[i] != -Inf) pm$row_lower_[i] <- ceiling(pm$row_lower_[i] - self$feastol)
      if (pm$row_upper_[i] != Inf) pm$row_upper_[i] <- floor(pm$row_upper_[i] + self$feastol)
    }
    self$rowintegral[i] <- integral
    self$maxAbsRowCoef[i] <- maxabsval
  }
  .hmip_HighsObjectiveFunction_setupCliquePartition(self$objectiveFunction, self$domain, self$cliquetable)
  .hmip_HighsDomain_setupObjectivePropagation(self$domain)
  .hmip_HighsDomain_computeRowActivities(self$domain)
  .hmip_HighsDomain_propagate(self$domain)
  if (.hmip_HighsDomain_infeasible(self$domain)) {
    ms$modelstatus_ <- "Infeasible"
    prev <- self$lower_bound
    self$lower_bound <- Inf
    if (!ms$submip && self$lower_bound != prev)
      .hmip_HighsMipSolverData_updatePrimalDualIntegral(self, prev, self$lower_bound,
                                                        self$upper_bound, self$upper_bound)
    self$pruned_treeweight <- c(1, 0)
    return(invisible())
  }
  if (model$num_col_ == 0L) {
    ms$modelstatus_ <- "Optimal"
    return(invisible())
  }
  if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
  for (col in .hmip_HighsDomain_getChangedCols(self$domain))
    .hmip_HighsImplications_cleanupVarbounds(self$implications, col)
  .hmip_HighsDomain_clearChangedCols(self$domain)
  lps <- .hmip_HighsLpRelaxation_getLpSolver(self$lp)
  .hmip_Highs_setOptionValue(lps, "presolve", "off")
  .hmip_Highs_setOptionValue(lps, "simplex_initial_condition_check", FALSE)
  .hmip_HighsMipSolverData_checkObjIntegrality(self)
  self$rootlpsol <- numeric(); self$firstlpsol <- numeric()
  num_binary <- 0L; num_domain_fixed <- 0L
  self$maxTreeSizeLog2 <- 0L
  dom <- self$domain
  for (i in seq_len(model$num_col_)) {
    vt <- model$integrality_[i]
    if (vt == 0L) {
      if (.hmip_HighsDomain_isFixed(dom, i - 1L)) { num_domain_fixed <- num_domain_fixed + 1L; next }
      self$continuous_cols <- c(self$continuous_cols, i - 1L)
    } else if (vt == 4L) {
      if (.hmip_HighsDomain_isFixed(dom, i - 1L)) { num_domain_fixed <- num_domain_fixed + 1L; next }
      self$implint_cols <- c(self$implint_cols, i - 1L)
      self$integral_cols <- c(self$integral_cols, i - 1L)
    } else if (vt == 1L) {
      if (.hmip_HighsDomain_isFixed(dom, i - 1L)) {
        num_domain_fixed <- num_domain_fixed + 1L
        if (.hmip_fractionality(dom$col_lower_[i]) > self$feastol) {
          ms$modelstatus_ <- "Infeasible"
          prev <- self$lower_bound
          self$lower_bound <- Inf
          if (!ms$submip && self$lower_bound != prev)
            .hmip_HighsMipSolverData_updatePrimalDualIntegral(self, prev, self$lower_bound,
                                                              self$upper_bound, self$upper_bound)
          self$pruned_treeweight <- c(1, 0)
          return(invisible())
        }
        next
      }
      self$integer_cols <- c(self$integer_cols, i - 1L)
      self$integral_cols <- c(self$integral_cols, i - 1L)
      self$maxTreeSizeLog2 <- self$maxTreeSizeLog2 +
        as.integer(ceiling(log2(min(1024, 1 + model$col_upper_[i] - model$col_lower_[i]))))
      num_binary <- num_binary + as.integer(model$col_lower_[i] == 0 && model$col_upper_[i] == 1)
    } else {
      stop("Unexpected variable type", call. = FALSE)
    }
  }
  .hmip_HighsMipSolverData_basisTransfer(self)
  self$numintegercols <- length(self$integer_cols)
  self$detectSymmetries <- self$detectSymmetries && num_binary > 0L
  self$numCliqueEntriesAfterPresolve <- .hmip_HighsCliqueTable_getNumEntries(self$cliquetable)
  nnz <- .hmip_numNonzero(ms)
  pl <- function(k) if (k == 1) "" else "s"
  if (self$numRestarts == 0L) {
    self$numCliqueEntriesAfterFirstPresolve <- self$numCliqueEntriesAfterPresolve
    .hmip_log(ms, "\nSolving MIP model with:")
    .hmip_log(ms, "   %d row%s", model$num_row_, pl(model$num_row_))
    .hmip_log(ms, "   %d col%s (%d binary, %d integer, %d implied int., %d continuous, %d domain fixed)",
              model$num_col_, pl(model$num_col_), num_binary, self$numintegercols - num_binary,
              length(self$implint_cols), length(self$continuous_cols), num_domain_fixed)
    .hmip_log(ms, "   %d nonzero%s", nnz, pl(nnz))
  } else {
    .hmip_log(ms, "Model after restart has %d row%s, %d col%s (%d bin., %d int., %d impl., %d cont., %d dom.fix.), and %d nonzero%s",
              model$num_row_, pl(model$num_row_), model$num_col_, pl(model$num_col_), num_binary,
              self$numintegercols - num_binary, length(self$implint_cols), length(self$continuous_cols),
              num_domain_fixed, nnz, pl(nnz))
  }
  .hmip_HighsPrimalHeuristics_setupIntCols(self$heuristics)
  if (self$upper_limit == Inf) self$analyticCenterComputed <- FALSE
  self$analyticCenterStatus <- "Notset"
  self$analyticCenter <- numeric()
  self$symmetries <- NULL
  invisible()
}

# transformNewIntegerFeasibleSolution: primal postsolve into the original
# space, with the repair LP (Highs on the model with integers fixed) when
# the postsolved point violates the original model
.hmip_HighsMipSolverData_transformNewIntegerFeasibleSolution <- function(self, sol,
                                                                       possibly_store_as_new_incumbent = TRUE) {
  ms <- self$mipsolver; opts <- ms$options_mip_; orig <- ms$orig_model_
  col_value <- .hmip_HighsPostsolveStack_undoPrimal(self$postSolveStack, opts, sol)
  row_value <- .hmip_row_values_quad(orig, col_value)
  allow_try_again <- TRUE
  repeat {
    r <- .hmip_HighsMipSolver_solutionFeasible(ms, orig, col_value, row_value)
    obj_value <- .cd_val(r$obj)
    if (!r$feasible && allow_try_again) {
      fl <- orig$col_lower_; fu <- orig$col_upper_
      for (i in seq_len(orig$num_col_)) if (orig$integrality_[i] == 1L) {
        solval <- .hmip_round(col_value[i])
        fl[i] <- max(fl[i], solval); fu[i] <- min(fu[i], solval)
      }
      self$total_repair_lp <- self$total_repair_lp + 1
      rep <- .hs_highs_run(orig$col_cost_, .hmip_HighsLp_dense(orig), orig$row_lower_, orig$row_upper_,
                           fl, fu, list(ptol = opts$mip_feasibility_tolerance,
                                        presolve = !isTRUE(opts$mip_root_presolve_only)))
      self$total_repair_lp_iterations <- rep$nit
      if (isTRUE(rep$primal_feasible) && !is.null(rep$col_value)) {
        self$total_repair_lp_feasible <- self$total_repair_lp_feasible + 1
        col_value <- rep$col_value
        row_value <- rep$row_value
        allow_try_again <- FALSE
        next
      }
    }
    break
  }
  transformed <- .cd_val(.cd_add(.cd_mul(r$obj, orig$sense_), -ms$model_$offset_))
  if (transformed >= self$upper_bound && length(sol)) return(transformed)
  if (possibly_store_as_new_incumbent) {
    store <- function() {
      ms$row_violation_ <- r$row_violation; ms$bound_violation_ <- r$bound_violation
      ms$integrality_violation_ <- r$integrality_violation
      ms$solution_ <- col_value; ms$solution_objective_ <- obj_value
    }
    if (r$feasible) {
      store()
    } else {
      tol <- opts$mip_feasibility_tolerance
      currentFeasible <- ms$solution_objective_ != Inf && ms$bound_violation_ <= tol &&
        ms$integrality_violation_ <= tol && ms$row_violation_ <= tol
      .hmip_log(ms, "Solution with objective %g has untransformed violations: bound = %.4g; integrality = %.4g; row = %.4g",
                obj_value, r$bound_violation, r$integrality_violation, r$row_violation)
      if (!currentFeasible) store()
      return(Inf)
    }
  }
  transformed
}

.hmip_HighsMipSolverData_percentageInactiveIntegers <- function(self) {
  nsub <- length(.hmip_HighsCliqueTable_getSubstitutions(self$cliquetable)$substcol)
  100 * (1 - (length(self$integer_cols) - nsub) / self$numintegercols)
}

.hmip_HighsMipSolverData_performRestart <- function(self) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  root_basis <- NULL
  ms$pscostinit <- .hmip_new_HighsPseudocostInitialization(self$pseudocost, opts$mip_pscost_minreliable,
                                                           self$postSolveStack)
  self$numRestarts <- self$numRestarts + 1L
  self$num_leaves_before_run <- self$num_leaves
  self$num_nodes_before_run <- self$num_nodes
  self$total_lp_iterations_before_run <- self$total_lp_iterations
  self$heuristic_lp_iterations_before_run <- self$heuristic_lp_iterations
  self$sepa_lp_iterations_before_run <- self$sepa_lp_iterations
  self$sb_lp_iterations_before_run <- self$sb_lp_iterations
  lplp <- .hmip_HighsLpRelaxation_getLp(self$lp)
  numCuts <- lplp$num_row_ - .hmip_numRow(ms)
  if (numCuts > 0L) .hmip_HighsPostsolveStack_appendCutsToModel(self$postSolveStack, numCuts)
  integrality <- self$presolvedModel$integrality_
  offset <- self$presolvedModel$offset_
  self$presolvedModel <- .hmip_HighsLp_copy(lplp)
  self$presolvedModel$offset_ <- offset
  self$presolvedModel$integrality_ <- integrality
  ms$model_ <- self$presolvedModel
  basis <- self$firstrootbasis
  if (isTRUE(basis$valid)) {
    ps <- self$postSolveStack
    cs <- rep("nonbasic", .hmip_HighsPostsolveStack_getOrigNumCol(ps))
    rs <- rep("basic", .hmip_HighsPostsolveStack_getOrigNumRow(ps))
    for (i in seq_len(ms$model_$num_col_))
      cs[.hmip_HighsPostsolveStack_getOrigColIndex(ps, i - 1L) + 1L] <- basis$col_status[i]
    for (i in seq_along(basis$row_status))
      rs[.hmip_HighsPostsolveStack_getOrigRowIndex(ps, i - 1L) + 1L] <- basis$row_status[i]
    root_basis <- list(valid = TRUE, useful = TRUE, col_status = cs, row_status = rs)
    ms$rootbasis <- root_basis
  }
  offs <- ms$model_$offset_
  self$upper_limit <- self$upper_limit + offs
  self$optimality_limit <- self$optimality_limit + offs
  self$upper_bound <- self$upper_bound + offs
  self$lower_bound <- self$lower_bound + offs
  self$incumbent <- numeric()
  self$pruned_treeweight <- c(0, 0)
  .hmip_HighsNodeQueue_clear(self$nodequeue)
  self$globalOrbits <- NULL
  num_reductions <- .hmip_HighsPostsolveStack_numReductions(self$postSolveStack)
  rlim <- opts$restart_presolve_reduction_limit
  further <- if (rlim >= 0) num_reductions + rlim else -1L
  .hmip_HighsMipSolverData_runMipPresolve(self, further)
  if (ms$modelstatus_ != "Notset") {
    self$upper_limit <- self$upper_limit - ms$model_$offset_
    self$optimality_limit <- self$optimality_limit - ms$model_$offset_
    if (ms$modelstatus_ == "Optimal") {
      self$upper_bound <- 0
      .hmip_HighsMipSolverData_transformNewIntegerFeasibleSolution(self, numeric())
    } else {
      self$upper_bound <- self$upper_bound - ms$model_$offset_
    }
    prev <- self$lower_bound - ms$model_$offset_
    self$lower_bound <- self$upper_bound
    if (!ms$submip && self$lower_bound != prev)
      .hmip_HighsMipSolverData_updatePrimalDualIntegral(self, prev, self$lower_bound,
                                                        self$upper_bound, self$upper_bound)
    if (ms$solution_objective_ != Inf && ms$modelstatus_ == "Infeasible") ms$modelstatus_ <- "Optimal"
    return(invisible())
  }
  .hmip_HighsMipSolverData_runSetup(self)
  if (.hmip_HighsMipSolver_terminate(ms)) return(invisible())
  .hmip_HighsPostsolveStack_removeCutsFromModel(self$postSolveStack, numCuts)
  if (!is.null(root_basis)) ms$rootbasis <- NULL
  ms$pscostinit <- NULL
  invisible()
}

.hmip_HighsMipSolverData_basisTransfer <- function(self) {
  ms <- self$mipsolver
  if (is.null(ms$rootbasis)) return(invisible())
  ps <- self$postSolveStack
  nr <- .hmip_numRow(ms); nc <- .hmip_numCol(ms)
  rs <- character(nr); cs <- character(nc)
  for (i in seq_len(nr)) rs[i] <- ms$rootbasis$row_status[.hmip_HighsPostsolveStack_getOrigRowIndex(ps, i - 1L) + 1L]
  for (i in seq_len(nc)) cs[i] <- ms$rootbasis$col_status[.hmip_HighsPostsolveStack_getOrigColIndex(ps, i - 1L) + 1L]
  self$firstrootbasis <- list(valid = TRUE, alien = TRUE, useful = TRUE, col_status = cs, row_status = rs)
  invisible()
}

.hmip_HighsMipSolverData_getSolution <- function(self) self$incumbent

.hmip_HighsMipSolverData_addIncumbent <- function(self, sol, solobj, solution_source,
                                                  print_display_line = TRUE, is_user_solution = FALSE) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  possibly_store <- solobj < self$upper_bound
  transformed <- if (possibly_store)
    .hmip_HighsMipSolverData_transformNewIntegerFeasibleSolution(self, sol, possibly_store) else 0
  if (possibly_store) {
    solobj <- transformed
    if (solobj >= self$upper_bound) return(FALSE)
    prev <- self$upper_bound
    self$upper_bound <- solobj
    if (!ms$submip && self$upper_bound != prev)
      .hmip_HighsMipSolverData_updatePrimalDualIntegral(self, self$lower_bound, self$lower_bound,
                                                        prev, self$upper_bound)
    self$incumbent <- sol
    nul <- .hmip_HighsMipSolverData_computeNewUpperLimit(self, solobj, 0, 0)
    if (!is_user_solution && !ms$submip) .hmip_HighsMipSolverData_saveReportMipSolution(self, nul)
    if (nul < self$upper_limit) {
      self$numImprovingSols <- self$numImprovingSols + 1L
      self$upper_limit <- nul
      self$optimality_limit <- .hmip_HighsMipSolverData_computeNewUpperLimit(
        self, solobj, opts$mip_abs_gap, opts$mip_rel_gap)
      .hmip_HighsNodeQueue_setOptimalityLimit(self$nodequeue, self$optimality_limit)
      .hmip_HighsDomain_propagate(self$domain)
      if (!.hmip_HighsDomain_infeasible(self$domain))
        .hmip_HighsRedcostFixing_propagateRootRedcost(self$redcostfixing, ms)
      if (.hmip_HighsDomain_infeasible(self$domain)) {
        self$pruned_treeweight <- c(1, 0)
        .hmip_HighsNodeQueue_clear(self$nodequeue)
        if (print_display_line) .hmip_HighsMipSolverData_printDisplayLine(self, solution_source)
        return(TRUE)
      }
      .hmip_HighsCliqueTable_extractObjCliques(self$cliquetable, ms)
      if (.hmip_HighsDomain_infeasible(self$domain)) {
        self$pruned_treeweight <- c(1, 0)
        .hmip_HighsNodeQueue_clear(self$nodequeue)
        if (print_display_line) .hmip_HighsMipSolverData_printDisplayLine(self, solution_source)
        return(TRUE)
      }
      self$pruned_treeweight <- .cd_add(self$pruned_treeweight,
        .hmip_HighsNodeQueue_performBounding(self$nodequeue, self$upper_limit))
      .hmip_HighsMipSolverData_printDisplayLine(self, solution_source)
    }
  } else if (!length(self$incumbent)) {
    self$incumbent <- sol
  }
  TRUE
}

# convertToPrintString for int64 and double values
.hmip_print_int <- function(val) {
  l <- floor(log10(max(1, val)))
  if (l <= 5) sprintf("%.0f", val) else if (l <= 8) sprintf("%.0fk", val %/% 1000)
  else sprintf("%.0fm", val %/% 1e6)
}
.hmip_print_double <- function(val, trailing = "") {
  l <- if (abs(val) == Inf) 0 else log10(max(1e-6, abs(val)))
  l <- trunc(l)
  fmt <- if (l >= 0 && l <= 3) "%.10g" else if (l == 4) "%.11g" else if (l == 5) "%.12g"
         else if (l >= 6 && l <= 10) "%.13g" else "%.9g"
  s <- if (abs(val) == Inf) (if (val > 0) "inf" else "-inf") else sprintf(fmt, val)
  paste0(s, trailing)
}

.hmip_HighsMipSolverData_printDisplayLine <- function(self, solution_source = -1L) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  if (!isTRUE(opts$output_flag)) return(invisible())
  timeless <- isTRUE(opts$timeless_log)
  self$disptime <- if (timeless) self$disptime + 1 else .hmip_timer_read(ms)
  if (solution_source == .hmip_src$None &&
      self$disptime - self$last_disptime < opts$mip_min_logging_interval) return(invisible())
  self$last_disptime <- self$disptime
  time_string <- if (timeless) "" else sprintf(" %7.1fs", self$disptime)
  if (self$num_disp_lines %% 20 == 0) {
    if (self$num_disp_lines == 0) {
      .hmip_log(ms, "\nSrc: B => Branching; C => Central rounding; F => Feasibility pump; H => Heuristic;")
      .hmip_log(ms, "     I => Shifting; J => Feasibility jump; L => Sub-MIP; P => Empty MIP; R => Randomized rounding;")
      .hmip_log(ms, "     S => Solve LP; T => Evaluate node; U => Unbounded; X => User solution; Y => HiGHS solution;")
      .hmip_log(ms, "     Z => ZI Round; l => Trivial lower; p => Trivial point; u => Trivial upper; z => Trivial zero")
    }
    w0 <- if (timeless) "   Work" else "      Work      "
    w1 <- if (timeless) "LpIters" else "LpIters     Time"
    .hmip_log(ms, "\n        Nodes      |    B&B Tree     |            Objective Bounds              |  Dynamic Constraints | %s\nSrc  Proc. InQueue |  Leaves   Expl. | BestBound       BestSol              Gap |   Cuts   InLp Confl. | %s\n",
              w0, w1)
  }
  self$num_disp_lines <- self$num_disp_lines + 1
  explored <- 100 * .cd_val(self$pruned_treeweight)
  g <- .hmip_HighsMipSolverData_limitsToGap(self, self$lower_bound, self$upper_bound)
  gap <- g$gap * 1e2; lb <- g$lb; ub <- g$ub
  if (opts$objective_bound < ub) ub <- opts$objective_bound
  # as in C++ the test below runs after ub was lowered, so it never marks ub
  star <- if (opts$objective_bound < ub) "*" else ""
  sense <- ms$orig_model_$sense_
  nrows <- .hmip_HighsLpRelaxation_numRows(self$lp)
  dyn <- if (nrows > 0) nrows - .hmip_HighsLpRelaxation_getNumModelRows(self$lp) else 0
  fmt_f2 <- function(v, w) if (is.infinite(v)) formatC(if (v > 0) "inf" else "-inf", width = w)
                           else sprintf(paste0("%", w, ".2f"), v)
  gap_string <- if (self$upper_bound != Inf) {
    sprintf("%8s", if (gap >= 9999) "Large" else sprintf("%.2f%%", gap))
  } else fmt_f2(gap, 8L)
  .hmip_log(ms, " %s %7s %7s   %7s %6.2f%%   %-15s %-15s %s   %6d %6d %6d   %7s%s",
            .hmip_HighsMipSolverData_solutionSourceToString(self, solution_source),
            .hmip_print_int(self$num_nodes), .hmip_print_int(.hmip_HighsNodeQueue_numActiveNodes(self$nodequeue)),
            .hmip_print_int(self$num_leaves - self$num_leaves_before_run), explored,
            .hmip_print_double(sense * lb), .hmip_print_double(sense * ub, star), gap_string,
            as.integer(.hmip_HighsCutPool_getNumCuts(self$cutpool)), as.integer(dyn),
            as.integer(.hmip_HighsConflictPool_getNumConflicts(self$conflictPool)),
            .hmip_print_int(self$total_lp_iterations), time_string)
  invisible()
}

.hmip_HighsMipSolverData_rootSeparationRound <- function(self, sepa, status) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  tmp <- -.hmip_HighsLpRelaxation_getNumLpIterations(self$lp)
  sr <- .hmip_HighsSeparation_separationRound(sepa, self$domain, status)
  ncuts <- sr$ncuts; status <- sr$status
  tmp <- tmp + .hmip_HighsLpRelaxation_getNumLpIterations(self$lp)
  self$avgrootlpiters <- .hmip_HighsLpRelaxation_getAvgSolveIters(self$lp)
  self$total_lp_iterations <- self$total_lp_iterations + tmp
  self$sepa_lp_iterations <- self$sepa_lp_iterations + tmp
  status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
  if (status == "Infeasible") return(list(stop = TRUE, ncuts = ncuts, status = status))
  solvals <- .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(self$lp))$col_value
  if (ms$submip || !length(self$incumbent)) {
    .hmip_HighsPrimalHeuristics_randomizedRounding(self$heuristics, solvals)
    if (isTRUE(opts$mip_heuristic_run_shifting)) .hmip_HighsPrimalHeuristics_shifting(self$heuristics, solvals)
    .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
    status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
    if (status == "Infeasible") return(list(stop = TRUE, ncuts = ncuts, status = status))
  }
  list(stop = FALSE, ncuts = ncuts, status = status)
}

# the "count this root node as a pruned leaf" bookkeeping of evaluateRootLp
.hmip_root_pruned <- function(self) {
  self$pruned_treeweight <- c(1, 0)
  self$num_nodes <- self$num_nodes + 1
  self$num_leaves <- self$num_leaves + 1
}

.hmip_HighsMipSolverData_evaluateRootLp <- function(self) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  lp <- self$lp; dom <- self$domain
  repeat {
    .hmip_HighsDomain_propagate(dom)
    if (!is.null(self$globalOrbits) && !.hmip_HighsDomain_infeasible(dom))
      .hmip_StabilizerOrbits_orbitalFixing(self$globalOrbits, dom)
    if (.hmip_HighsDomain_infeasible(dom)) {
      .hmip_set_lower_bound(self, min(Inf, self$upper_bound))
      .hmip_root_pruned(self)
      return("Infeasible")
    }
    lpBoundsChanged <- FALSE
    if (length(.hmip_HighsDomain_getChangedCols(dom))) {
      lpBoundsChanged <- TRUE
      .hmip_HighsMipSolverData_removeFixedIndices(self)
      .hmip_HighsLpRelaxation_flushDomain(lp, dom)
    }
    lpWasSolved <- FALSE
    lps <- .hmip_HighsLpRelaxation_getLpSolver(lp)
    if (lpBoundsChanged || .hmip_Highs_getModelStatus(lps) == "Notset") {
      it <- -.hmip_HighsLpRelaxation_getNumLpIterations(lp)
      status <- .hmip_HighsLpRelaxation_resolveLp(lp, dom)
      it <- it + .hmip_HighsLpRelaxation_getNumLpIterations(lp)
      self$total_lp_iterations <- self$total_lp_iterations + it
      self$avgrootlpiters <- .hmip_HighsLpRelaxation_getAvgSolveIters(lp)
      lpWasSolved <- TRUE
      if (status == "Unbounded") {
        ms$modelstatus_ <- if (!length(ms$solution_)) "UnboundedOrInfeasible" else "Unbounded"
        .hmip_root_pruned(self)
        return(status)
      }
      if (status == "Optimal" && !length(.hmip_HighsLpRelaxation_getFractionalIntegers(lp)) &&
          .hmip_HighsMipSolverData_addIncumbent(self, .hmip_Highs_getSolution(lps)$col_value,
                                                .hmip_HighsLpRelaxation_getObjective(lp),
                                                .hmip_src$EvaluateNode)) {
        ms$modelstatus_ <- "Optimal"
        .hmip_set_lower_bound(self, self$upper_bound)
        .hmip_root_pruned(self)
        return("Infeasible")
      }
      if (status == "Optimal" && isTRUE(opts$mip_heuristic_run_zi_round))
        .hmip_HighsPrimalHeuristics_ziRound(self$heuristics, .hmip_Highs_getSolution(lps)$col_value)
    } else {
      status <- .hmip_HighsLpRelaxation_getStatus(lp)
    }
    if (status == "Infeasible") {
      .hmip_set_lower_bound(self, min(Inf, self$upper_bound))
      .hmip_root_pruned(self)
      return(status)
    }
    if (.hmip_HighsLpRelaxation_unscaledDualFeasible(.hmip_HighsLpRelaxation_getStatus(lp))) {
      .hmip_set_lower_bound(self, max(.hmip_HighsLpRelaxation_getObjective(lp), self$lower_bound))
      if (lpWasSolved) {
        .hmip_HighsRedcostFixing_addRootRedcost(self$redcostfixing, ms, .hmip_Highs_getSolution(lps)$col_dual,
                                                .hmip_HighsLpRelaxation_getObjective(lp))
        if (self$upper_limit != Inf) .hmip_HighsRedcostFixing_propagateRootRedcost(self$redcostfixing, ms)
      }
    }
    if (self$lower_bound > self$optimality_limit) {
      .hmip_root_pruned(self)
      return("Infeasible")
    }
    if (!length(.hmip_HighsDomain_getChangedCols(dom))) return(status)
  }
}

.hmip_HighsMipSolverData_evaluateRootNode <- function(self) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  lp <- self$lp
  maxSepaRounds <- if (ms$submip) 5L else .hmip_kHighsIInf
  if (self$numRestarts == 0L) maxSepaRounds <- min(trunc(2 * sqrt(self$maxTreeSizeLog2)), maxSepaRounds)
  repeat {  # restart:
    self$symData <- NULL
    if (self$detectSymmetries) .hmip_HighsMipSolverData_startSymmetryDetection(self)
    if (!self$analyticCenterComputed) .hmip_HighsMipSolverData_startAnalyticCenterComputation(self)
    .hmip_HighsLpRelaxation_setIterationLimit(lp)
    .hmip_HighsLpRelaxation_loadModel(lp)
    .hmip_HighsDomain_clearChangedCols(self$domain)
    .hmip_HighsLpRelaxation_setObjectiveLimit(lp, self$upper_limit)
    .hmip_set_lower_bound(self, max(self$lower_bound, .hmip_HighsDomain_getObjectiveLowerBound(self$domain)))
    .hmip_HighsMipSolverData_printDisplayLine(self)
    lps <- .hmip_HighsLpRelaxation_getLpSolver(lp)
    if (isTRUE(self$firstrootbasis$valid)) .hmip_Highs_setBasis(lps, self$firstrootbasis)
    else if (isTRUE(opts$mip_root_presolve_only)) .hmip_Highs_setOptionValue(lps, "presolve", "off")
    else .hmip_Highs_setOptionValue(lps, "presolve", "on")
    status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
    if (self$numRestarts == 0L) self$firstrootlpiters <- self$total_lp_iterations
    .hmip_Highs_setOptionValue(lps, "output_flag", FALSE)
    .hmip_Highs_setOptionValue(lps, "presolve", "off")
    .hmip_Highs_setOptionValue(lps, "parallel", "off")
    if (status %in% c("Infeasible", "Unbounded")) return(invisible())
    self$firstlpsol <- .hmip_HighsLpRelaxation_getSolution(lp)$col_value
    self$firstlpsolobj <- .hmip_HighsLpRelaxation_getObjective(lp)
    self$rootlpsolobj <- self$firstlpsolobj
    b <- .hmip_Highs_getBasis(lps)
    if (isTRUE(b$valid) && .hmip_HighsLpRelaxation_numRows(lp) == .hmip_numRow(ms)) {
      self$firstrootbasis <- b
    } else {
      self$firstrootbasis <- list(valid = TRUE, alien = FALSE, useful = TRUE,
                                  col_status = rep("nonbasic", .hmip_numCol(ms)),
                                  row_status = rep("basic", .hmip_numRow(ms)))
    }
    if (.hmip_HighsCutPool_getNumCuts(self$cutpool) != 0) {
      cutset <- .hmip_new_HighsCutSet()
      .hmip_HighsCutPool_separateLpCutsAfterRestart(self$cutpool, cutset)
      .hmip_HighsLpRelaxation_addCuts(lp, cutset)
      status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
      .hmip_HighsLpRelaxation_removeObsoleteRows(lp)
      if (status == "Infeasible") return(invisible())
    }
    .hmip_HighsLpRelaxation_setIterationLimit(lp, max(10000, trunc(10 * self$avgrootlpiters)))
    self$last_disptime <- -Inf
    self$disptime <- 0
    if (isTRUE(opts$mip_heuristic_run_zi_round)) .hmip_HighsPrimalHeuristics_ziRound(self$heuristics, self$firstlpsol)
    .hmip_HighsPrimalHeuristics_randomizedRounding(self$heuristics, self$firstlpsol)
    if (isTRUE(opts$mip_heuristic_run_shifting)) .hmip_HighsPrimalHeuristics_shifting(self$heuristics, self$firstlpsol)
    .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
    status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
    if (status == "Infeasible") return(invisible())
    self$rootlpsolobj <- self$firstlpsolobj
    .hmip_HighsMipSolverData_removeFixedIndices(self)
    if (isTRUE(opts$mip_allow_restart) && opts$presolve != "off") {
      fixingRate <- .hmip_HighsMipSolverData_percentageInactiveIntegers(self)
      if (isTRUE(fixingRate >= 10)) {
        self$analyticCenterPending <- FALSE  # tg.cancel()
        .hmip_log(ms, "\n%.1f%% inactive integer columns, restarting", fixingRate)
        .hmip_HighsMipSolverData_performRestart(self)
        self$numRestartsRoot <- self$numRestartsRoot + 1L
        if (ms$modelstatus_ == "Notset") next
        return(invisible())
      }
    }
    avgdirection <- numeric(.hmip_numCol(ms)); curdirection <- numeric(.hmip_numCol(ms))
    stall <- 0L; smoothprogress <- 0; nseparounds <- 0L
    sepa <- .hmip_new_HighsSeparation(ms)
    .hmip_HighsSeparation_setLpRelaxation(sepa, lp)
    nfrac <- function() length(.hmip_HighsLpRelaxation_getFractionalIntegers(lp))
    while (.hmip_HighsLpRelaxation_scaledOptimal(status) && nfrac() > 0L && stall < 3L) {
      .hmip_HighsMipSolverData_printDisplayLine(self)
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      if (nseparounds == maxSepaRounds) break
      .hmip_HighsMipSolverData_removeFixedIndices(self)
      if (!ms$submip && opts$presolve != "off") {
        if (isTRUE(.hmip_HighsMipSolverData_percentageInactiveIntegers(self) >= 10)) { stall <- -1L; break }
      }
      nseparounds <- nseparounds + 1L
      rs <- .hmip_HighsMipSolverData_rootSeparationRound(self, sepa, status)
      ncuts <- rs$ncuts; status <- rs$status
      if (rs$stop) return(invisible())
      if (nseparounds >= 5L && !ms$submip && !self$analyticCenterComputed) {
        if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
        .hmip_HighsMipSolverData_finishAnalyticCenterComputation(self)
        .hmip_HighsPrimalHeuristics_centralRounding(self$heuristics)
        .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
        if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
        status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
        if (status == "Infeasible") return(invisible())
      }
      solvals <- .hmip_HighsLpRelaxation_getSolution(lp)$col_value
      sqrnorm <- c(0, 0)
      for (i in seq_along(curdirection)) {
        curdirection[i] <- self$firstlpsol[i] - solvals[i]
        sqrnorm <- .cd_add(sqrnorm, curdirection[i] * curdirection[i])
      }
      scale <- .cd_val(.cd_div_cd(c(1, 0), .cd_sqrt(sqrnorm)))
      sqrnorm <- c(0, 0); dotproduct <- c(0, 0)
      for (i in seq_along(curdirection)) {
        avgdirection[i] <- (scale * curdirection[i] - avgdirection[i]) / nseparounds
        sqrnorm <- .cd_add(sqrnorm, avgdirection[i] * avgdirection[i])
        dotproduct <- .cd_add(dotproduct, avgdirection[i] * curdirection[i])
      }
      progress <- .cd_val(.cd_div_cd(dotproduct, .cd_sqrt(sqrnorm)))
      if (nseparounds == 1L) {
        smoothprogress <- progress
      } else {
        alpha <- 1 / 3
        nextprogress <- (1 - alpha) * smoothprogress + alpha * progress
        if (nextprogress < smoothprogress * 1.01 &&
            (.hmip_HighsLpRelaxation_getObjective(lp) - self$firstlpsolobj) <=
            (self$rootlpsolobj - self$firstlpsolobj) * 1.001) stall <- stall + 1L
        else stall <- 0L
        smoothprogress <- nextprogress
      }
      self$rootlpsolobj <- .hmip_HighsLpRelaxation_getObjective(lp)
      .hmip_HighsLpRelaxation_setIterationLimit(lp, max(10000, trunc(10 * self$avgrootlpiters)))
      if (ncuts == 0L) break
    }
    .hmip_HighsLpRelaxation_setIterationLimit(lp)
    status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
    if (status == "Infeasible") return(invisible())
    self$rootlpsol <- .hmip_Highs_getSolution(lps)$col_value
    self$rootlpsolobj <- .hmip_HighsLpRelaxation_getObjective(lp)
    .hmip_HighsLpRelaxation_setIterationLimit(lp, max(10000, trunc(10 * self$avgrootlpiters)))
    if (isTRUE(opts$mip_heuristic_run_zi_round)) {
      .hmip_HighsPrimalHeuristics_ziRound(self$heuristics, self$firstlpsol)
      .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
    }
    if (isTRUE(opts$mip_heuristic_run_shifting)) {
      .hmip_HighsPrimalHeuristics_shifting(self$heuristics, self$rootlpsol)
      .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
    }
    # one more separation round when evaluateRootLp sees new global bounds
    sep_once <- function() {
      separate <- length(.hmip_HighsDomain_getChangedCols(self$domain)) > 0L
      status <<- .hmip_HighsMipSolverData_evaluateRootLp(self)
      if (status == "Infeasible") return(TRUE)
      if (separate && .hmip_HighsLpRelaxation_scaledOptimal(status)) {
        rs <- .hmip_HighsMipSolverData_rootSeparationRound(self, sepa, status)
        status <<- rs$status
        if (rs$stop) return(TRUE)
        nseparounds <<- nseparounds + 1L
        .hmip_HighsMipSolverData_printDisplayLine(self)
      }
      FALSE
    }
    if (!self$analyticCenterComputed) {
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      .hmip_HighsMipSolverData_finishAnalyticCenterComputation(self)
      .hmip_HighsPrimalHeuristics_centralRounding(self$heuristics)
      .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      if (sep_once()) return(invisible())
    }
    .hmip_HighsMipSolverData_printDisplayLine(self)
    if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
    repeat {  # do { ... } while (false)
      if (!length(self$rootlpsol)) break
      if (self$upper_limit != Inf && !.hmip_HighsMipSolverData_moreHeuristicsAllowed(self)) break
      if (isTRUE(opts$mip_heuristic_run_root_reduced_cost)) {
        .hmip_HighsPrimalHeuristics_rootReducedCost(self$heuristics)
        .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
      }
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      if (sep_once()) return(invisible())
      if (self$upper_limit != Inf && !.hmip_HighsMipSolverData_moreHeuristicsAllowed(self)) break
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      if (isTRUE(opts$mip_heuristic_run_rens)) {
        .hmip_HighsPrimalHeuristics_RENS(self$heuristics, self$rootlpsol)
        .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
      }
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      if (sep_once()) return(invisible())
      if (self$upper_limit != Inf || ms$submip) break
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      .hmip_HighsPrimalHeuristics_feasibilityPump(self$heuristics)
      .hmip_HighsPrimalHeuristics_flushStatistics(self$heuristics)
      if (.hmip_HighsMipSolverData_checkLimits(self)) return(invisible())
      status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
      if (status == "Infeasible") return(invisible())
      break
    }
    if (self$lower_bound > self$upper_limit) {
      ms$modelstatus_ <- "Optimal"
      .hmip_root_pruned(self)
      return(invisible())
    }
    if (sep_once()) return(invisible())
    .hmip_HighsMipSolverData_removeFixedIndices(self)
    if (isTRUE(.hmip_Highs_getBasis(lps)$valid)) .hmip_HighsLpRelaxation_removeObsoleteRows(lp)
    self$rootlpsolobj <- .hmip_HighsLpRelaxation_getObjective(lp)
    .hmip_HighsMipSolverData_printDisplayLine(self)
    if (self$lower_bound <= self$upper_limit) {
      if (!ms$submip && isTRUE(opts$mip_allow_restart) && opts$presolve != "off") {
        if (!self$analyticCenterComputed) .hmip_HighsMipSolverData_finishAnalyticCenterComputation(self)
        fixingRate <- .hmip_HighsMipSolverData_percentageInactiveIntegers(self)
        if (isTRUE(fixingRate >= 2.5 + 7.5 * ms$submip) || (!ms$submip && isTRUE(fixingRate > 0) && self$numRestarts == 0L)) {
          .hmip_log(ms, "\n%.1f%% inactive integer columns, restarting", fixingRate)
          if (stall != -1L) maxSepaRounds <- min(maxSepaRounds, nseparounds)
          .hmip_HighsMipSolverData_performRestart(self)
          if (.hmip_HighsMipSolver_terminate(ms)) return(invisible())
          self$numRestartsRoot <- self$numRestartsRoot + 1L
          if (ms$modelstatus_ == "Notset") next
          return(invisible())
        }
      }
      if (self$detectSymmetries) {
        .hmip_HighsMipSolverData_finishSymmetryDetection(self)
        status <- .hmip_HighsMipSolverData_evaluateRootLp(self)
        if (status == "Infeasible") return(invisible())
      }
      .hmip_HighsNodeQueue_emplaceNode(self$nodequeue, .hmip_domchg(numeric(), integer(), integer()),
                                       integer(), self$lower_bound,
                                       .hmip_HighsLpRelaxation_computeBestEstimate(lp, self$pseudocost), 1L)
    }
    return(invisible())
  }
}

.hmip_HighsMipSolverData_checkLimits <- function(self, nodeOffset = 0) {
  ms <- self$mipsolver; opts <- ms$options_mip_
  if (.hmip_HighsMipSolverData_terminatorActive(self) && .hmip_HighsMipSolverData_terminatorTerminated(self))
    return(TRUE)
  set <- function(st) if (ms$modelstatus_ == "Notset") ms$modelstatus_ <- st
  if (!ms$submip && ms$solution_objective_ < Inf && opts$objective_target > -Inf) {
    s <- ms$orig_model_$sense_
    if (s * ms$solution_objective_ < s * opts$objective_target) { set("ObjectiveTarget"); return(TRUE) }
  }
  if (.hmip_limit_set(opts$mip_max_nodes) && self$num_nodes + nodeOffset >= opts$mip_max_nodes) {
    set("SolutionLimit"); return(TRUE)
  }
  if (.hmip_limit_set(opts$mip_max_leaves) && self$num_leaves >= opts$mip_max_leaves) {
    set("SolutionLimit"); return(TRUE)
  }
  if (.hmip_limit_set(opts$mip_max_improving_sols) && self$numImprovingSols >= opts$mip_max_improving_sols) {
    set("SolutionLimit"); return(TRUE)
  }
  if (opts$time_limit < Inf && .hmip_timer_read(ms) >= opts$time_limit) {
    set("TimeLimit"); return(TRUE)
  }
  FALSE
}

.hmip_HighsMipSolverData_checkObjIntegrality <- function(self) {
  .hmip_HighsObjectiveFunction_checkIntegrality(self$objectiveFunction, self$epsilon)
  if (.hmip_HighsObjectiveFunction_isIntegral(self$objectiveFunction) && self$numRestarts == 0L)
    .hmip_log(self$mipsolver, "Objective function is integral with scale %g",
              .hmip_HighsObjectiveFunction_integralScale(self$objectiveFunction))
  invisible()
}

.hmip_HighsMipSolverData_setupDomainPropagation <- function(self) {
  ms <- self$mipsolver; model <- ms$model_; m <- model$a_matrix_
  tr <- .hmip_sparse_transpose(model$num_row_, model$num_col_, m$start_, m$index_, m$value_)
  self$ARstart_ <- tr$ARstart_; self$ARindex_ <- tr$ARindex_; self$ARvalue_ <- tr$ARvalue_
  .hmip_HighsPseudocost_init(self$pseudocost, ms)
  self$maxAbsRowCoef <- vapply(seq_len(model$num_row_), function(i) {
    s <- self$ARstart_[i]; e <- self$ARstart_[i + 1L]
    if (e > s) max(0, abs(self$ARvalue_[(s + 1L):e])) else 0
  }, 0)
  .hmip_HighsDomain_assign(self$domain, .hmip_new_HighsDomain(ms))
  .hmip_HighsDomain_computeRowActivities(self$domain)
  invisible()
}

.hmip_HighsMipSolverData_saveReportMipSolution <- function(self, new_upper_limit = -Inf) {
  ms <- self$mipsolver
  if (ms$submip || new_upper_limit >= self$upper_limit) return(invisible())
  if (isTRUE(ms$options_mip_$mip_improving_solution_save))
    ms$saved_objective_and_solution_[[length(ms$saved_objective_and_solution_) + 1L]] <-
      list(objective = ms$solution_objective_, col_value = ms$solution_)
  invisible()
}

.hmip_HighsMipSolverData_limitsToBounds <- function(self) {
  g <- .hmip_HighsMipSolverData_limitsToGap(self, self$lower_bound, self$upper_bound)
  db <- g$lb; pb <- min(self$mipsolver$options_mip_$objective_bound, g$ub)
  if (self$mipsolver$orig_model_$sense_ == -1L) { db <- -db; pb <- -pb }
  list(dual_bound = db, primal_bound = pb, mip_rel_gap = g$gap)
}

.hmip_HighsMipSolverData_queryExternalSolution <- function(self, ...) invisible()
.hmip_HighsMipSolverData_terminatorConcurrency <- function(self) self$mipsolver$terminator_$num_instance
.hmip_HighsMipSolverData_terminatorActive <- function(self) .hmip_HighsMipSolverData_terminatorConcurrency(self) > 0L
.hmip_HighsMipSolverData_terminatorTerminated <- function(self) {
  self$mipsolver$termination_status_ != "Notset"
}

.hmip_possInfRelDiff <- function(v0, v1, den) {
  if (abs(v0) == Inf) return(if (abs(v1) == Inf) 0 else Inf)
  if (abs(v1) == Inf) return(Inf)
  abs(v1 - v0) / max(1, abs(den))
}

.hmip_HighsMipSolverData_updatePrimalDualIntegral <- function(self, from_lower_bound, to_lower_bound,
                                                              from_upper_bound, to_upper_bound,
                                                              check_bound_change = TRUE,
                                                              check_prev_data = TRUE) {
  pdi <- self$primal_dual_integral
  from <- .hmip_HighsMipSolverData_limitsToGap(self, from_lower_bound, from_upper_bound)
  to <- .hmip_HighsMipSolverData_limitsToGap(self, to_lower_bound, to_upper_bound)
  if (pdi$value > -Inf) {
    if (to$gap < Inf) {
      time <- .hmip_timer_read(self$mipsolver)
      if (from$gap < Inf) pdi$value <- pdi$value + (time - pdi$prev_time) * pdi$prev_gap
      pdi$prev_time <- time
    }
  } else {
    pdi$value <- 0
  }
  pdi$prev_lb <- to$lb; pdi$prev_ub <- to$ub; pdi$prev_gap <- to$gap
  self$primal_dual_integral <- pdi
  invisible()
}

# ---------------------------------------------------------------------------
# Highs::callSolveMip with withoutSemiVariables

# withoutSemiVariables: each semicontinuous or semi-integer column x gets a
# binary y with rows x - l y >= 0 and x - u y <= 0, and lower bound 0
.hmip_withoutSemiVariables <- function(lp) {
  out <- .hmip_HighsLp_copy(lp)
  semi <- which(lp$integrality_ %in% c(2L, 3L))
  if (!length(semi)) return(out)
  m <- lp$a_matrix_
  n <- lp$num_col_; nr <- lp$num_row_
  start <- integer(n + 1L); index <- integer(); value <- numeric()
  row_num <- nr
  for (j in seq_len(n)) {
    k <- if (m$start_[j + 1L] > m$start_[j]) (m$start_[j] + 1L):m$start_[j + 1L] else integer()
    index <- c(index, m$index_[k]); value <- c(value, m$value_[k])
    if (j %in% semi) {
      index <- c(index, row_num, row_num + 1L); value <- c(value, 1, 1)
      row_num <- row_num + 2L
    }
    start[j + 1L] <- length(index)
  }
  cost <- lp$col_cost_; cl <- lp$col_lower_; cu <- lp$col_upper_
  rl <- lp$row_lower_; ru <- lp$row_upper_; integ <- lp$integrality_
  row_num <- nr
  for (j in semi) {
    cost <- c(cost, 0); cl <- c(cl, 0); cu <- c(cu, 1)
    rl <- c(rl, 0, -Inf); ru <- c(ru, Inf, 0)
    index <- c(index, row_num, row_num + 1L); value <- c(value, -lp$col_lower_[j], -lp$col_upper_[j])
    row_num <- row_num + 2L
    start <- c(start, length(index))
    integ <- c(integ, 1L)
    integ[j] <- if (lp$integrality_[j] == 2L) 0L else 1L
    cl[j] <- 0
  }
  .hmip_new_HighsLp(cost, cl, cu, rl, ru, start, index, value, integ, lp$offset_, lp$sense_)
}

# Highs::callSolveMip: returns the model status, the solution in the space
# of lp and the MIP information SciPy reports
.hmip_callSolveMip <- function(lp, options = list()) {
  use_lp <- if (any(lp$integrality_ %in% c(2L, 3L))) .hmip_withoutSemiVariables(lp) else lp
  solver <- .hmip_new_HighsMipSolver(options, use_lp)
  .hmip_HighsMipSolver_run(solver)
  out <- list(model_status = solver$modelstatus_, mip_node_count = solver$node_count_,
              mip_dual_bound = solver$dual_bound_, mip_gap = solver$gap_,
              simplex_iteration_count = solver$total_lp_iterations_,
              objective = solver$solution_objective_, log = solver$log_lines,
              mipsolver = solver)
  if (solver$solution_objective_ != Inf) {
    x <- solver$solution_[seq_len(lp$num_col_)]
    out$col_value <- x
    out$row_value <- .hmip_row_values_quad(lp, x)
  }
  out
}
