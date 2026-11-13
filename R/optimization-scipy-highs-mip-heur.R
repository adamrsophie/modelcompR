# The primal heuristics of the HiGHS 1.12.0 MIP solver, as SciPy 1.18.1
# vendors them: mip/HighsPrimalHeuristics.cpp (RENS, RINS, root reduced
# cost, randomized rounding, shifting, ZI round, central rounding, line
# search rounding, the feasibility pump and flushStatistics),
# mip/HighsFeasibilityJump.cpp and the feasibility jump library HiGHS
# vendors in mip/feasibilityjump.hh, together with the methods of
# util/HighsRandom.h the existing port lacks and the std::mt19937 and
# std::uniform_real_distribution<double> of libstdc++ that the feasibility
# jump draws from. MIT license notice in R/optimization-scipy-highs-util.R.
#
# Conventions of the MIP port (scratchpad API-B1.md): column and row ids
# stored in vectors are 0-based as in C++, read as v[col + 1L]; HighsCDouble
# values are c(hi, lo); solution sources are the integers of .hmip_src.

# HighsRandom::closedFraction: the 53 bit output times the C++ constant
.hs_rng_closed_fraction <- function(r) {
  .hs_rng_advance(r)
  lo <- .hs_lo32(r$state); hi <- .hs_hi32(r$state)
  a <- .hs_u64_shift(.hs_pair_hash(0L, lo, hi), 11L)
  b <- .hs_u64_shift(.hs_pair_hash(1L, lo, hi), 32L)
  x <- .hs_u64_xor(a, b)
  out <- x[1L] + 65536 * (x[2L] + 65536 * (x[3L] + 65536 * x[4L]))
  # the C++ literal 1.1102230246251566e-16 is the double 2^-53
  out * 0x1p-53
}

# HighsRandom::real(a, b) = a + (b - a) * closedFraction()
.hs_rng_real <- function(r, a, b) a + (b - a) * .hs_rng_closed_fraction(r)

# HighsRandom::bit: the top bit of pair_hash<0>
.hs_rng_bit <- function(r) {
  .hs_rng_advance(r)
  .hs_pair_hash(0L, .hs_lo32(r$state), .hs_hi32(r$state))[4L] >= 32768
}

# std::mt19937 on doubles holding 32 bit unsigned values. The generator is a
# closure so that its state vector is updated in place.
.hs_u32_xor <- function(a, b) {
  bitwXor(a %/% 65536, b %/% 65536) * 65536 + bitwXor(a %% 65536, b %% 65536)
}
.hs_u32_and <- function(a, b) {
  bitwAnd(a %/% 65536, b %/% 65536) * 65536 + bitwAnd(a %% 65536, b %% 65536)
}
.hs_mt19937 <- function(seed = 5489) {
  mt <- numeric(624L)
  mt[1L] <- seed %% 4294967296
  for (i in 1:623) {
    x <- .hs_u32_xor(mt[i], mt[i] %/% 1073741824)
    # 1812433253 * x + i modulo 2^32 with 16 bit halves
    xh <- x %/% 65536; xl <- x %% 65536
    v <- ((27655 * xl + 35173 * xh) %% 65536) * 65536 + 35173 * xl + i
    mt[i + 1L] <- v %% 4294967296
  }
  out <- numeric(624L)
  pos <- 625L
  twist <- function(k) {
    # k: 1-based positions whose updates only read values already final
    nxt <- mt[(k %% 624L) + 1L]
    y <- (mt[k] %/% 2147483648) * 2147483648 + nxt %% 2147483648
    v <- .hs_u32_xor(mt[((k + 396L) %% 624L) + 1L], y %/% 2)
    odd <- y %% 2 == 1
    v[odd] <- .hs_u32_xor(v[odd], 2567483615)
    mt[k] <<- v
  }
  refill <- function() {
    twist(1:227); twist(228:454); twist(455:623); twist(624L)
    y <- mt
    y <- .hs_u32_xor(y, y %/% 2048)
    y <- .hs_u32_xor(y, .hs_u32_and((y * 128) %% 4294967296, 2636928640))
    y <- .hs_u32_xor(y, .hs_u32_and((y * 32768) %% 4294967296, 4022730752))
    out <<- .hs_u32_xor(y, y %/% 262144)
    pos <<- 1L
  }
  function() {
    if (pos > 624L) refill()
    v <- out[pos]
    pos <<- pos + 1L
    v
  }
}

# std::uniform_real_distribution<double>(0, 1)(rng) of libstdc++:
# generate_canonical<double, 53> takes two 32 bit draws, sums them as
# g1 + g2 * 2^32 in double, divides by 2^64 and clamps below 1
.hs_uniform01 <- function(rng) {
  g1 <- rng(); g2 <- rng()
  r <- (g1 + g2 * 4294967296) / 18446744073709551616
  if (r >= 1) r <- 0x1.fffffffffffffp-1
  r
}

# ---------------------------------------------------------------------------
# external_feasibilityjump (mip/feasibilityjump.hh). One closure holds the
# Problem, the JumpMove and the FeasibilityJumpSolver state; variables and
# constraints are 1-based here (the library only uses them as positions).
# Sense codes: 0 Equal, 1 Lte, 2 Gte. Logging has no effect on results and
# is left out.
.hmip_fj_new <- function(seed = 0L, equalityTolerance = 1e-5, violationTolerance = 1e-5,
                         weightUpdateDecay = 1.0) {
  # Problem
  vtype <- logical(); vlb <- numeric(); vub <- numeric(); vobj <- numeric()
  vcons <- list(); vcoef <- list()
  csense <- integer(); crhs <- numeric(); cweight <- numeric(); clhs <- numeric(); cviol <- integer()
  cvars <- list(); ccoef <- list()
  incumbent <- numeric()
  violated <- integer(); nviol <- 0L
  usedRelaxContinuous <- FALSE
  nNonzeros <- 0
  incumbentObjective <- NaN
  # JumpMove
  mvalue <- numeric(); mscore <- numeric()
  # FeasibilityJumpSolver
  goodVars <- integer(); ngood <- 0L; goodIdx <- integer()
  rng <- .hs_mt19937(seed)
  bestObjective <- Inf
  objectiveWeight <- 0
  bestViolationScore <- Inf
  effortAtLastCallback <- 0
  effortAtLastImprovement <- 0
  totalEffort <- 0
  weightUpdateIncrement <- 1
  nBumps <- 0
  randomVarProbability <- 0.001
  randomCellProbability <- 0.01
  maxMovesToEvaluate <- 25L

  score <- function(c, lhs) {
    s <- csense[c]
    if (s == 0L) -abs(lhs - crhs[c])
    else if (s == 1L) { d <- lhs - crhs[c]; -(if (0 < d) d else 0) }
    else { d <- crhs[c] - lhs; -(if (0 < d) d else 0) }
  }

  addVar <- function(isInt, lb, ub, objCoeff) {
    goodIdx[length(goodIdx) + 1L] <<- 0L
    k <- length(vlb) + 1L
    vtype[k] <<- isInt; vlb[k] <<- lb; vub[k] <<- ub; vobj[k] <<- objCoeff
    vcons[k] <<- list(integer()); vcoef[k] <<- list(numeric())
    incumbent[k] <<- lb
    k - 1L
  }

  # rowVarIdxs are 0-based as the HiGHS caller passes them
  addConstraint <- function(sense, rhs, rowVarIdxs, rowCoeffs, relax_continuous = 0L) {
    if (relax_continuous) usedRelaxContinuous <<- TRUE
    if (relax_continuous > 0 && sense == 0L && any(!vtype[rowVarIdxs + 1L])) {
      addConstraint(2L, rhs, rowVarIdxs, rowCoeffs, relax_continuous)
      addConstraint(1L, rhs, rowVarIdxs, rowCoeffs, relax_continuous)
      return(.Machine$integer.max)
    }
    idx <- integer(); val <- numeric()
    for (i in seq_along(rowVarIdxs)) {
      v <- rowVarIdxs[i] + 1L
      if (relax_continuous > 0 && !vtype[v]) {
        a <- rowCoeffs[i]
        if (sense == 1L) {
          rhs <- rhs - (if (a >= 0) a * vlb[v] else a * vub[v])
        } else if (sense == 2L) {
          rhs <- rhs - (if (a >= 0) a * vub[v] else a * vlb[v])
        } else return(-.Machine$integer.max - 1)
      } else {
        idx[length(idx) + 1L] <- v; val[length(val) + 1L] <- rowCoeffs[i]
      }
    }
    if (!length(idx)) {
      ok <- if (sense == 1L) 0 <= rhs + equalityTolerance
            else if (sense == 2L) 0 + equalityTolerance >= rhs
            else abs(0 - rhs) < equalityTolerance
      return(if (ok) .Machine$integer.max else -.Machine$integer.max - 1)
    }
    c <- length(crhs) + 1L
    for (t in seq_along(idx)) {
      v <- idx[t]
      vcons[[v]][length(vcons[[v]]) + 1L] <<- c
      vcoef[[v]][length(vcoef[[v]]) + 1L] <<- val[t]
    }
    nNonzeros <<- nNonzeros + length(idx)
    csense[c] <<- sense; crhs[c] <<- rhs; cweight[c] <<- 1; clhs[c] <<- NaN; cviol[c] <<- 0L
    cvars[c] <<- list(idx); ccoef[c] <<- list(val)
    c - 1L
  }

  resetIncumbent <- function(initialValues) {
    if (!is.null(initialValues)) incumbent <<- as.numeric(initialValues)
    o <- 0
    for (i in seq_along(vlb)) o <- o + vobj[i] * incumbent[i]
    incumbentObjective <<- o
    violated <<- integer(); nviol <<- 0L
    for (c in seq_along(crhs)) {
      l <- 0
      vs <- cvars[[c]]; cs <- ccoef[[c]]
      for (t in seq_along(vs)) l <- l + cs[t] * incumbent[vs[t]]
      clhs[c] <<- l
      if (score(c, l) < -violationTolerance) {
        nviol <<- nviol + 1L
        violated[nviol] <<- c
        cviol[c] <<- nviol
      } else cviol[c] <<- 0L
    }
  }

  # JumpMove::updateValue
  updateValue <- function(v) {
    bf <- numeric(); bw <- numeric(); nb <- 0L
    inc <- incumbent[v]
    lower <- vlb[v]; upper <- vub[v]
    currentValue <- lower
    currentSlope <- 0
    cs <- vcons[[v]]; co <- vcoef[[v]]
    for (t in seq_along(cs)) {
      c <- cs[t]; a <- co[t]; w <- cweight[c]; s <- csense[c]; rhs <- crhs[c]
      if (s == 1L) { b1 <- -Inf; b2 <- rhs }
      else if (s == 2L) { b1 <- rhs; b2 <- Inf }
      else { b1 <- c(-Inf, rhs, rhs); b2 <- c(rhs, rhs, Inf) }
      for (k in seq_along(b1)) {
        residual <- clhs[c] - a * inc
        r1 <- (1.0 / a) * (b1[k] - residual)
        r2 <- (1.0 / a) * (b2[k] - residual)
        if (vtype[v]) {
          r1 <- ceiling(r1 - equalityTolerance)
          r2 <- floor(r2 + equalityTolerance)
        }
        if (r1 > r2) next
        if (r1 > currentValue) {
          currentSlope <- currentSlope - w
          if (r1 < upper) { nb <- nb + 1L; bf[nb] <- r1; bw[nb] <- w }
        }
        if (r2 <= currentValue) {
          currentSlope <- currentSlope + w
        } else if (r2 < upper) { nb <- nb + 1L; bf[nb] <- r2; bw[nb] <- w }
      }
    }
    if (is.finite(lower)) { nb <- nb + 1L; bf[nb] <- lower; bw[nb] <- 0 }
    if (is.finite(upper)) { nb <- nb + 1L; bf[nb] <- upper; bw[nb] <- 0 }
    if (nb > 1L) {
      # std::sort on pairs (equal pairs are indistinguishable); the buffers
      # are short, so an insertion sort beats order()
      for (k in 2:nb) {
        f <- bf[k]; w <- bw[k]; p <- k - 1L
        while (p >= 1L && (f < bf[p] || (!(bf[p] < f) && w < bw[p]))) {
          bf[p + 1L] <- bf[p]; bw[p + 1L] <- bw[p]; p <- p - 1L
        }
        bf[p + 1L] <- f; bw[p + 1L] <- w
      }
    }
    currentScore <- 0
    currentValue <- if (nb == 0L) inc else bf[1L]
    bestScore <- currentScore
    bestValue <- currentValue
    for (k in seq_len(nb)) {
      currentScore <- currentScore + (bf[k] - currentValue) * currentSlope
      currentSlope <- currentSlope + bw[k]
      currentValue <- bf[k]
      if (abs(bestValue - inc) < equalityTolerance ||
          (!(abs(currentValue - inc) < equalityTolerance) && currentScore < bestScore)) {
        bestScore <- currentScore
        bestValue <- currentValue
      }
      if (!(abs(bestValue - inc) < equalityTolerance) && currentSlope >= 0) break
    }
    mvalue[v] <<- bestValue
  }

  # bestMove over the single jump move per variable: Move::undef unless the
  # score exceeds -inf
  bestScoreOf <- function(v) { s <- mscore[v]; if (s > -Inf) s else -Inf }

  updateGoodMoves <- function(v) {
    anyGood <- bestScoreOf(v) > 0
    if (anyGood && goodIdx[v] == 0L) {
      ngood <<- ngood + 1L
      goodVars[ngood] <<- v
      goodIdx[v] <<- ngood
    } else if (!anyGood && goodIdx[v] != 0L) {
      lastVar <- goodVars[ngood]
      this <- goodIdx[v]
      goodVars[this] <<- lastVar
      goodVars[ngood] <<- v
      goodIdx[lastVar] <<- this
      goodIdx[v] <<- 0L
      ngood <<- ngood - 1L
    }
  }

  resetMoves <- function(v) {
    cs <- vcons[[v]]; co <- vcoef[[v]]
    totalEffort <<- totalEffort + length(cs)
    updateValue(v)
    d <- mvalue[v] - incumbent[v]
    s <- 0
    s <- s + objectiveWeight * vobj[v] * d
    for (t in seq_along(cs)) {
      c <- cs[t]
      candidateLhs <- clhs[c] + co[t] * d
      s <- s + cweight[c] * (score(c, candidateLhs) - score(c, clhs[c]))
    }
    mscore[v] <<- s
    updateGoodMoves(v)
  }

  init <- function(initialValues) {
    resetIncumbent(initialValues)
    mvalue <<- numeric(length(vlb)); mscore <<- numeric(length(vlb))
    totalEffort <<- totalEffort + nNonzeros
    ngood <<- 0L; goodVars <<- integer()
    for (i in seq_along(vlb)) resetMoves(i)
  }

  updateWeights <- function() {
    nBumps <<- nBumps + 1
    rescaleAllWeights <- FALSE
    dt <- 0
    if (nviol == 0L) {
      objectiveWeight <<- objectiveWeight + weightUpdateIncrement
      if (objectiveWeight > 1.0e20) rescaleAllWeights <- TRUE
      dt <- dt + length(vlb)
      for (v in seq_along(vlb)) {
        mscore[v] <<- mscore[v] + weightUpdateIncrement * vobj[v] * (mvalue[v] - incumbent[v])
      }
    } else {
      for (k in seq_len(nviol)) {
        c <- violated[k]
        cweight[c] <<- cweight[c] + weightUpdateIncrement
        if (cweight[c] > 1.0e20) rescaleAllWeights <- TRUE
        vs <- cvars[[c]]; co <- ccoef[[c]]
        dt <- dt + length(vs)
        for (t in seq_along(vs)) {
          v <- vs[t]
          candidateLhs <- clhs[c] + co[t] * (mvalue[v] - incumbent[v])
          diff <- weightUpdateIncrement * (score(c, candidateLhs) - score(c, clhs[c]))
          mscore[v] <<- mscore[v] + diff
          updateGoodMoves(v)
        }
      }
    }
    weightUpdateIncrement <<- weightUpdateIncrement / weightUpdateDecay
    if (rescaleAllWeights) {
      weightUpdateIncrement <<- weightUpdateIncrement * 1.0e-20
      objectiveWeight <<- objectiveWeight * 1.0e-20
      cweight <<- cweight * 1.0e-20
      dt <- dt + length(crhs)
      for (i in seq_along(vlb)) resetMoves(i)
    }
    totalEffort <<- totalEffort + dt
  }

  selectVariable <- function() {
    if (ngood > 0L) {
      if (.hs_uniform01(rng) < randomVarProbability) return(goodVars[rng() %% ngood + 1])
      sampleSize <- min(maxMovesToEvaluate, ngood)
      totalEffort <<- totalEffort + sampleSize
      bestScore <- -Inf
      bestVar <- NA_integer_
      for (i in seq_len(sampleSize)) {
        v <- goodVars[rng() %% ngood + 1]
        s <- bestScoreOf(v)
        if (s > bestScore) { bestScore <- s; bestVar <- v }
      }
      return(bestVar)
    }
    updateWeights()
    if (nviol > 0L) {
      c <- violated[rng() %% nviol + 1]
      vs <- cvars[[c]]
      if (.hs_uniform01(rng) < randomCellProbability) return(vs[rng() %% length(vs) + 1])
      bestScore <- -Inf
      bestVar <- NA_integer_
      for (v in vs) {
        s <- bestScoreOf(v)
        if (s > bestScore) { bestScore <- s; bestVar <- v }
      }
      return(bestVar)
    }
    as.integer(rng() %% length(vlb) + 1)
  }

  # doVariableMove with Problem::setValue and modifyMove inlined
  doVariableMove <- function(v) {
    newValue <- if (mscore[v] > -Inf) mvalue[v] else NaN
    dt <- 0
    oldValue <- incumbent[v]
    delta <- newValue - oldValue
    incumbent[v] <<- newValue
    incumbentObjective <<- incumbentObjective + vobj[v] * delta
    cs <- vcons[[v]]; co <- vcoef[[v]]
    for (t in seq_along(cs)) {
      c <- cs[t]
      oldLhs <- clhs[c]
      newLhs <- oldLhs + co[t] * delta
      clhs[c] <<- newLhs
      newCost <- score(c, newLhs)
      if (newCost < -violationTolerance && cviol[c] == 0L) {
        nviol <<- nviol + 1L
        violated[nviol] <<- c
        cviol[c] <<- nviol
      }
      if (newCost >= -violationTolerance && cviol[c] != 0L) {
        lastC <- violated[nviol]
        this <- cviol[c]
        violated[this] <<- lastC
        violated[nviol] <<- c
        cviol[lastC] <<- this
        cviol[c] <<- 0L
        nviol <<- nviol - 1L
      }
      vs <- cvars[[c]]; vc <- ccoef[[c]]
      dt <- dt + length(vs)
      w <- cweight[c]
      sOld <- score(c, oldLhs); sNew <- score(c, newLhs)
      se <- csense[c]; rhs <- crhs[c]
      for (k in seq_along(vs)) {
        u <- vs[k]
        if (u == v) next
        d <- vc[k] * (mvalue[u] - incumbent[u])
        # score() inlined for the two modified left hand sides
        l1 <- oldLhs + d; l2 <- newLhs + d
        if (se == 0L) {
          a1 <- -abs(l1 - rhs); a2 <- -abs(l2 - rhs)
        } else if (se == 1L) {
          t1 <- l1 - rhs; t2 <- l2 - rhs
          a1 <- -(if (0 < t1) t1 else 0); a2 <- -(if (0 < t2) t2 else 0)
        } else {
          t1 <- rhs - l1; t2 <- rhs - l2
          a1 <- -(if (0 < t1) t1 else 0); a2 <- -(if (0 < t2) t2 else 0)
        }
        oldScoreTerm <- w * (a1 - sOld)
        newScoreTerm <- w * (a2 - sNew)
        mscore[u] <<- mscore[u] + (newScoreTerm - oldScoreTerm)
        updateGoodMoves(u)
      }
    }
    totalEffort <<- totalEffort + dt
    resetMoves(v)
  }

  user_terminate <- function(callback, solution) {
    if (!is.null(solution) || totalEffort - effortAtLastCallback > 500000) {
      effortAtLastCallback <<- totalEffort
      status <- list(totalEffort = totalEffort,
                     effortSinceLastImprovement = totalEffort - effortAtLastImprovement,
                     solution = solution, numVars = length(vlb),
                     solutionObjectiveValue = incumbentObjective)
      if (isTRUE(callback(status))) return(TRUE)
    }
    FALSE
  }

  # callback(status) returns TRUE for CallbackControlFlow::Terminate; trace,
  # when a function, receives (step, var, value, totalEffort, nviol, ngood,
  # nBumps, objective) after every move
  solve <- function(initialValues, callback, trace = NULL) {
    init(initialValues)
    step <- 0L
    while (step < .Machine$integer.max) {
      if (user_terminate(callback, NULL)) break
      if (nviol < bestViolationScore) {
        effortAtLastImprovement <<- totalEffort
        bestViolationScore <<- nviol
      }
      if (nviol == 0L && incumbentObjective < bestObjective) {
        a <- abs(incumbentObjective)
        rel <- (bestObjective - incumbentObjective) / (if (1.0 < a) a else 1.0)
        if (rel > 1e-4) {
          effortAtLastImprovement <<- totalEffort
          bestObjective <<- incumbentObjective
          if (user_terminate(callback, incumbent)) break
        }
      }
      if (!length(vlb)) break
      v <- selectVariable()
      doVariableMove(v)
      if (!is.null(trace)) trace(step, v - 1L, incumbent[v], totalEffort, nviol, ngood, nBumps, incumbentObjective)
      step <- step + 1L
    }
    0L
  }

  list(addVar = addVar, addConstraint = addConstraint, solve = solve,
       state = function() list(totalEffort = totalEffort, nBumps = nBumps, incumbent = incumbent,
                               incumbentObjective = incumbentObjective, nviol = nviol))
}

# HighsMipSolverData::feasibilityJump (mip/HighsFeasibilityJump.cpp), built
# with 32 bit HighsInt. Returns the model status string ("Notset" or
# "Infeasible").
.hmip_HighsMipSolverData_feasibilityJump <- function(self) {
  mipsolver <- self$mipsolver
  model <- mipsolver$model_
  sense_multiplier <- as.numeric(model$sense_)
  n <- model$num_col_
  feastol <- self$feastol
  found <- FALSE
  col_value <- numeric(n)
  objective_function_value <- NaN
  use_incumbent <- length(self$incumbent) > 0L
  solver <- .hmip_fj_new(mipsolver$options_mip_$random_seed, self$epsilon, feastol)
  for (col in seq_len(n)) {
    lower <- model$col_lower_[col]; upper <- model$col_upper_[col]
    isint <- model$integrality_[col] != 0L
    if (isint) {
      lower <- ceiling(lower - feastol)
      upper <- floor(upper + feastol)
    }
    legal <- !is.nan(lower) && !is.nan(upper) && lower <= upper && lower < Inf && upper > -Inf
    if (!legal) return("Infeasible")
    solver$addVar(isint, lower, upper, sense_multiplier * model$col_cost_[col])
    init <- 0
    if (use_incumbent && is.finite(self$incumbent[col])) {
      x <- self$incumbent[col]
      m <- if (x < upper) x else upper
      init <- if (lower < m) m else lower
    } else if (is.finite(lower)) init <- lower
    else if (is.finite(upper)) init <- upper
    col_value[col] <- init
  }
  # HighsSparseMatrix::createRowwise: rows hold their entries in column order
  st <- model$a_matrix_$start_; ix <- model$a_matrix_$index_; vx <- model$a_matrix_$value_
  m <- model$num_row_
  rows_i <- vector("list", m); rows_v <- vector("list", m)
  for (i in seq_len(m)) { rows_i[[i]] <- integer(); rows_v[[i]] <- numeric() }
  for (col in seq_len(n)) {
    if (st[col + 1L] > st[col]) for (k in (st[col] + 1L):st[col + 1L]) {
      r <- ix[k] + 1L
      rows_i[[r]][length(rows_i[[r]]) + 1L] <- col - 1L
      rows_v[[r]][length(rows_v[[r]]) + 1L] <- vx[k]
    }
  }
  for (i in seq_len(m)) {
    lo <- model$row_lower_[i]; up <- model$row_upper_[i]
    if (is.finite(lo)) solver$addConstraint(2L, lo, rows_i[[i]], rows_v[[i]], 0L)
    if (is.finite(up)) solver$addConstraint(1L, up, rows_i[[i]], rows_v[[i]], 0L)
  }
  nnz <- st[n + 1L]
  kMaxTotalEffort <- nnz * 1024
  kMaxEffortSinceLastImprovement <- nnz * 256
  callback <- function(status) {
    if (!is.null(status$solution)) {
      found <<- TRUE
      col_value <<- status$solution
      objective_function_value <<- model$offset_ + sense_multiplier * status$solutionObjectiveValue
    }
    status$effortSinceLastImprovement > kMaxEffortSinceLastImprovement ||
      status$totalEffort > kMaxTotalEffort
  }
  solver$solve(col_value, callback)
  if (found) .hmip_HighsMipSolverData_trySolution(self, col_value, .hmip_src$FeasibilityJump)
  "Notset"
}

# ---------------------------------------------------------------------------
# HighsPrimalHeuristics (mip/HighsPrimalHeuristics.h/.cpp)

# HighsDomain::Reason::branching()
.hmip_Reason_branching <- function() list(type = -1L, index = 0L)

# member-wise copy of an environment (a C++ copy constructor over plain
# members)
.hmip_heur_copy_env <- function(e) list2env(as.list(e, all.names = TRUE), parent = parent.env(e))

# std::min(a, b) and std::max(a, b): the first argument unless the second
# compares strictly smaller or larger (NaN keeps the first)
.hmip_min2 <- function(a, b) if (b < a) b else a
.hmip_max2 <- function(a, b) if (a < b) b else a

.hmip_new_HighsPrimalHeuristics <- function(mipsolver) {
  h <- new.env(parent = emptyenv())
  h$mipsolver <- mipsolver
  h$total_repair_lp <- 0
  h$total_repair_lp_feasible <- 0
  h$total_repair_lp_iterations <- 0
  h$lp_iterations <- 0
  h$randgen <- .hs_rng(mipsolver$options_mip_$random_seed)
  h$successObservations <- 0
  h$numSuccessObservations <- 0L
  h$infeasObservations <- 0
  h$numInfeasObservations <- 0L
  h$intcols <- integer()
  h
}

.hmip_HighsPrimalHeuristics_setupIntCols <- function(self) {
  md <- self$mipsolver$mipdata_
  ic <- md$integer_cols
  ft <- md$feastol
  lock <- clq <- numeric(length(ic)); hsh <- vector("list", length(ic))
  for (k in seq_along(ic)) {
    c <- ic[k]
    lock[k] <- (ft + md$uplocks[c + 1L]) * (ft + md$downlocks[c + 1L])
    clq[k] <- (ft + .hmip_HighsCliqueTable_getNumImplications2(md$cliquetable, c, 1L)) *
      (ft + .hmip_HighsCliqueTable_getNumImplications2(md$cliquetable, c, 0L))
    hsh[[k]] <- .hmip_hash(c)
  }
  # a strict total order, so the result does not depend on the pdqsort path
  less <- function(a, b) {
    if (lock[a] > lock[b]) return(TRUE)
    if (lock[b] > lock[a]) return(FALSE)
    # std::tuple (cliqueScore, hash, col) compared with >
    if (clq[a] > clq[b]) return(TRUE)
    if (clq[b] > clq[a]) return(FALSE)
    hc <- .hmip_u64_cmp(hsh[[a]], hsh[[b]])
    if (hc != 0L) return(hc > 0L)
    ic[a] > ic[b]
  }
  ord <- .hmip_pdqsort(seq_along(ic), less)
  self$intcols <- ic[unlist(ord)]
  invisible(self)
}

.hmip_HighsPrimalHeuristics_solveSubMip <- function(self, lp, basis, fixingRate, colLower, colUpper,
                                                    maxleaves, maxnodes, stallnodes) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  submipoptions <- mipsolver$options_mip_
  submip <- list2env(as.list(lp, all.names = TRUE), parent = emptyenv())
  submip$col_lower_ <- colLower
  submip$col_upper_ <- colUpper
  submip$integrality_ <- mipsolver$model_$integrality_
  submip$offset_ <- 0
  submipoptions$mip_max_leaves <- maxleaves
  submipoptions$output_flag <- FALSE
  submipoptions$mip_max_nodes <- maxnodes
  submipoptions$mip_max_stall_nodes <- stallnodes
  submipoptions$mip_pscost_minreliable <- 0L
  # time_limit -= timer_.read(): the port has no clock (B1 timer stub), so
  # the sub-MIP keeps the parent's limit
  submipoptions$objective_bound <- md$upper_limit
  if (!isTRUE(mipsolver$submip)) {
    curr_abs_gap <- md$upper_limit - md$lower_bound
    if (curr_abs_gap == Inf) {
      curr_abs_gap <- abs(md$lower_bound)
      if (curr_abs_gap == Inf) curr_abs_gap <- 0.0
    }
    submipoptions$mip_rel_gap <- 0.0
    submipoptions$mip_abs_gap <- md$feastol * .hmip_max2(curr_abs_gap, 1000.0)
  }
  submipoptions$presolve <- if (isTRUE(submipoptions$mip_root_presolve_only)) "off" else "on"
  submipoptions$mip_detect_symmetry <- FALSE
  submipoptions$mip_heuristic_effort <- 0.8
  sub <- .hmip_new_HighsMipSolver(submipoptions, submip, NULL, TRUE, mipsolver$submip_level + 1L)
  .hmip_HighsMipSolver_initialiseTerminator(sub, mipsolver)
  sub$rootbasis <- basis
  sub$pscostinit <- .hmip_new_HighsPseudocostInitialization(md$pseudocost, 1L)
  sub$clqtableinit <- md$cliquetable
  sub$implicinit <- md$implications
  .hmip_HighsMipSolver_run(sub)
  mipsolver$max_submip_level <- max(sub$max_submip_level + 1L, mipsolver$max_submip_level)
  if (!identical(sub$termination_status_, "Notset")) {
    mipsolver$termination_status_ <- sub$termination_status_
    return(FALSE)
  }
  if (!is.null(sub$mipdata_)) {
    nint <- length(md$integral_cols); ncont <- length(md$continuous_cols)
    numUnfixed <- nint + ncont
    adjustmentfactor <- ((1 - fixingRate) * nint + ncont) / .hmip_max2(1.0, numUnfixed)
    self$lp_iterations <- self$lp_iterations + trunc(adjustmentfactor * sub$mipdata_$total_lp_iterations)
    self$total_repair_lp <- self$total_repair_lp + sub$mipdata_$total_repair_lp
    self$total_repair_lp_feasible <- self$total_repair_lp_feasible + sub$mipdata_$total_repair_lp_feasible
    self$total_repair_lp_iterations <- self$total_repair_lp_iterations + sub$mipdata_$total_repair_lp_iterations
    if (isTRUE(mipsolver$submip))
      md$num_nodes <- md$num_nodes + max(1, trunc(adjustmentfactor * sub$node_count_))
  }
  infeasible <- identical(sub$modelstatus_, "Infeasible")
  if (infeasible) {
    self$infeasObservations <- self$infeasObservations + fixingRate
    self$numInfeasObservations <- self$numInfeasObservations + 1L
  }
  if (sub$node_count_ <= 1 && infeasible) return(FALSE)
  oldNumImprovingSols <- md$numImprovingSols
  if (!infeasible && length(sub$solution_))
    .hmip_HighsMipSolverData_trySolution(md, sub$solution_, .hmip_src$SubMip)
  if (md$numImprovingSols != oldNumImprovingSols) {
    self$successObservations <- self$successObservations + fixingRate
    self$numSuccessObservations <- self$numSuccessObservations + 1L
  }
  TRUE
}

.hmip_HighsPrimalHeuristics_determineTargetFixingRate <- function(self) {
  lowFixingRate <- 0.6
  highFixingRate <- 0.6
  if (self$numInfeasObservations != 0L) {
    infeasRate <- self$infeasObservations / self$numInfeasObservations
    highFixingRate <- 0.9 * infeasRate
    lowFixingRate <- .hmip_min2(lowFixingRate, highFixingRate)
  }
  if (self$numSuccessObservations != 0L) {
    successFixingRate <- self$successObservations / self$numSuccessObservations
    lowFixingRate <- .hmip_min2(lowFixingRate, 0.9 * successFixingRate)
    highFixingRate <- .hmip_max2(successFixingRate * 1.1, highFixingRate)
  }
  .hs_rng_real(self$randgen, lowFixingRate, highFixingRate)
}

# HeuristicNeighbourhood: of the hash set fixedCols only the size is read
.hmip_new_HeuristicNeighbourhood <- function(mipsolver, localdom) {
  nb <- new.env(parent = emptyenv())
  nb$localdom <- localdom
  start <- length(.hmip_HighsDomain_getDomainChangeStack(localdom)$column)
  nb$startCheckedChanges <- start
  nb$nCheckedChanges <- start
  ic <- mipsolver$mipdata_$integral_cols
  numFixed <- 0L
  for (i in ic) if (localdom$col_lower_[i + 1L] == localdom$col_upper_[i + 1L]) numFixed <- numFixed + 1L
  nb$numTotal <- length(ic) - numFixed
  nb$fixed <- logical(mipsolver$model_$num_col_)
  nb$nfixed <- 0L
  nb
}
.hmip_HeuristicNeighbourhood_getFixingRate <- function(nb) {
  dom <- nb$localdom
  cols <- .hmip_HighsDomain_getDomainChangeStack(dom)$column
  while (nb$nCheckedChanges < length(cols)) {
    col <- cols[nb$nCheckedChanges + 1L]
    nb$nCheckedChanges <- nb$nCheckedChanges + 1L
    if (.hmip_HighsDomain_variableType(dom, col) == 0L) next
    if (.hmip_HighsDomain_isFixed(dom, col) && !nb$fixed[col + 1L]) {
      nb$fixed[col + 1L] <- TRUE
      nb$nfixed <- nb$nfixed + 1L
    }
  }
  if (nb$numTotal) nb$nfixed / nb$numTotal else 0.0
}
.hmip_HeuristicNeighbourhood_backtracked <- function(nb) {
  nb$nCheckedChanges <- nb$startCheckedChanges
  if (nb$nfixed) { nb$fixed[] <- FALSE; nb$nfixed <- 0L }
}

.hmip_HighsPrimalHeuristics_rootReducedCost <- function(self) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  lurkingBounds <- .hmip_HighsRedcostFixing_getLurkingBounds(md$redcostfixing, mipsolver)
  if (10 * length(lurkingBounds) < length(md$integral_cols)) return(invisible())
  # ties in the key are possible, so the pdqsort path decides their order
  lurkingBounds <- .hmip_pdqsort(lurkingBounds, function(a, b) a$first > b$first)
  localdom <- .hmip_HighsDomain_copy(md$domain)
  neighbourhood <- .hmip_new_HeuristicNeighbourhood(mipsolver, localdom)
  currCutoff <- Inf
  lower_bound <- md$lower_bound + md$feastol
  for (domchg in lurkingBounds) {
    currCutoff <- domchg$first
    if (currCutoff <= lower_bound) break
    if (.hmip_HighsDomain_isActive(localdom, domchg$second)) next
    .hmip_HighsDomain_changeBound(localdom, domchg$second)
    repeat {
      .hmip_HighsDomain_propagate(localdom)
      if (.hmip_HighsDomain_infeasible(localdom)) {
        .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
        prev_lower_bound <- md$lower_bound
        md$lower_bound <- .hmip_max2(md$lower_bound, currCutoff)
        if (!isTRUE(mipsolver$submip) && md$lower_bound != prev_lower_bound)
          .hmip_HighsMipSolverData_updatePrimalDualIntegral(md, prev_lower_bound, md$lower_bound,
                                                            md$upper_bound, md$upper_bound)
        .hmip_HighsDomain_backtrack(localdom)
        if (.hmip_HighsDomain_getBranchDepth(localdom) == 0L) break
        .hmip_HeuristicNeighbourhood_backtracked(neighbourhood)
        next
      }
      break
    }
    if (.hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood) >= 0.5) break
  }
  fixingRate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
  if (fixingRate < 0.3) return(invisible())
  .hmip_HighsPrimalHeuristics_solveSubMip(self, mipsolver$model_, md$firstrootbasis, fixingRate,
                                          localdom$col_lower_, localdom$col_upper_, 500L,
                                          200 + md$num_nodes %/% 20, 12L)
  invisible()
}

# calcFixVal (static in HighsPrimalHeuristics.cpp)
.hmip_calcFixVal <- function(rootchange, fracval, cost) {
  if (rootchange >= 0.4) ceiling(fracval)
  else if (rootchange <= -0.4) floor(fracval)
  else if (cost > 0.0) ceiling(fracval)
  else if (cost < 0.0) floor(fracval)
  else floor(fracval + 0.5)
}

# pdqsort of the first k fractional integers of an LP relaxation by
# (|fixval - value|, hash((col << 32) + size)), written back in place
.hmip_heur_sort_fracints <- function(heurlp, k, fixval) {
  fi <- heurlp$fractionalints
  n <- length(fi$first)
  if (k < 1L) return(invisible())
  d <- numeric(k); hs <- vector("list", k)
  for (t in seq_len(k)) {
    d[t] <- abs(fixval(fi$first[t], fi$second[t]) - fi$second[t])
    hs[[t]] <- .hmip_hash(n, fi$first[t])
  }
  less <- function(a, b) {
    if (d[a] < d[b]) return(TRUE)
    if (d[b] < d[a]) return(FALSE)
    .hmip_u64_cmp(hs[[a]], hs[[b]]) < 0L
  }
  ord <- c(unlist(.hmip_pdqsort(seq_len(k), less)), seq_len(n - k) + k)
  heurlp$fractionalints <- list(first = fi$first[ord], second = fi$second[ord])
  invisible()
}

# the setup RENS and RINS share: a search with its own pseudocost copy
.hmip_heur_dive_setup <- function(self) {
  mipsolver <- self$mipsolver
  pscost <- .hmip_heur_copy_env(mipsolver$mipdata_$pseudocost)
  heur <- .hmip_new_HighsSearch(mipsolver, pscost)
  localdom <- .hmip_HighsSearch_getLocalDomain(heur)
  .hmip_HighsSearch_setHeuristic(heur, TRUE)
  list(heur = heur, localdom = localdom)
}
# ... and its own copy of the LP relaxation
.hmip_heur_dive_lp <- function(self, heur, localdom) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  heurlp <- .hmip_HighsLpRelaxation_copy(md$lp)
  .hmip_HighsLpRelaxation_setObjectiveLimit(heurlp, md$upper_limit)
  .hmip_HighsLpRelaxation_setAdjustSymmetricBranchingCol(heurlp, FALSE)
  .hmip_HighsSearch_setLpRelaxation(heur, heurlp)
  .hmip_Highs_changeColsBoundsRange(.hmip_HighsLpRelaxation_getLpSolver(heurlp), 0L,
                                    mipsolver$model_$num_col_ - 1L, localdom$col_lower_, localdom$col_upper_)
  .hmip_HighsDomain_clearChangedCols(localdom)
  .hmip_HighsSearch_createNewNode(heur)
  heurlp
}
# intcols.erase(remove_if(isFixed in the global domain)); permanent
.hmip_heur_drop_fixed_intcols <- function(self) {
  dom <- self$mipsolver$mipdata_$domain
  keep <- vapply(self$intcols, function(i) !.hmip_HighsDomain_isFixed(dom, i), TRUE)
  self$intcols <- self$intcols[keep]
}

# the end RENS and RINS share after the dive: depth first search or sub-MIP;
# returns action "return", or "retry" with the new maxfixingrate and
# targetdepth (the goto retry of the C++)
.hmip_heur_dive_finish <- function(self, heur, heurlp, localdom, neighbourhood) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  if (!.hmip_HighsSearch_hasNode(heur)) {
    self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
    return(list(action = "return"))
  }
  fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
  if (fixingrate < 0.1 || (isTRUE(mipsolver$submip) && md$numImprovingSols != 0)) {
    .hmip_HighsSearch_setMinReliable(heur, 0L)
    .hmip_HighsSearch_solveDepthFirst(heur, 10L)
    self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
    if (isTRUE(mipsolver$submip)) md$num_nodes <- md$num_nodes + .hmip_HighsSearch_getLocalNodes(heur)
    return(list(action = "return"))
  }
  .hmip_HighsLpRelaxation_removeObsoleteRows(heurlp, FALSE)
  ok <- .hmip_HighsPrimalHeuristics_solveSubMip(
    self, .hmip_HighsLpRelaxation_getLp(heurlp), .hmip_Highs_getBasis(.hmip_HighsLpRelaxation_getLpSolver(heurlp)),
    fixingrate, localdom$col_lower_, localdom$col_upper_, 500L, 200 + md$num_nodes %/% 20, 12L)
  if (.hmip_HighsMipSolverData_terminatorTerminated(md)) return(list(action = "return"))
  if (!ok) {
    new_lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
    if (new_lp_iterations + md$heuristic_lp_iterations >
        100000 + (md$total_lp_iterations - md$heuristic_lp_iterations - md$sb_lp_iterations) %/% 2) {
      self$lp_iterations <- new_lp_iterations
      return(list(action = "return"))
    }
    targetdepth <- .hmip_HighsSearch_getCurrentDepth(heur) %/% 2L
    if (targetdepth <= 1L || .hmip_HighsMipSolverData_checkLimits(md)) {
      self$lp_iterations <- new_lp_iterations
      return(list(action = "return"))
    }
    return(list(action = "retry", maxfixingrate = fixingrate * 0.5, targetdepth = targetdepth))
  }
  self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
  list(action = "return")
}

# branch col up or down to val, then propagate when asked; TRUE when the
# local domain became infeasible (after conflict analysis)
.hmip_heur_branch <- function(md, heur, localdom, up, col, val, point, propagate = TRUE) {
  if (up) .hmip_HighsSearch_branchUpwards(heur, col, val, point)
  else .hmip_HighsSearch_branchDownwards(heur, col, val, point)
  if (propagate) .hmip_HighsDomain_propagate(localdom)
  if (.hmip_HighsDomain_infeasible(localdom)) {
    .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
    return(TRUE)
  }
  FALSE
}

.hmip_HighsPrimalHeuristics_RENS <- function(self, tmp) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  if (.hmip_HighsDomain_infeasible(md$domain)) return(invisible())
  s <- .hmip_heur_dive_setup(self)
  heur <- s$heur; localdom <- s$localdom
  .hmip_heur_drop_fixed_intcols(self)
  heurlp <- .hmip_heur_dive_lp(self, heur, localdom)
  maxfixingrate <- .hmip_HighsPrimalHeuristics_determineTargetFixingRate(self)
  fixingrate <- 0.0
  stop <- FALSE
  targetdepth <- 1L
  nbacktracks <- -1L
  neighbourhood <- .hmip_new_HeuristicNeighbourhood(mipsolver, localdom)
  feastol <- md$feastol
  getFixVal <- function(col, fracval) {
    fixval <- .hmip_calcFixVal(if (!length(md$rootlpsol)) 0.0 else fracval - md$rootlpsol[col + 1L],
                               fracval, mipsolver$model_$col_cost_[col + 1L])
    fixval <- .hmip_min2(localdom$col_upper_[col + 1L], fixval)
    .hmip_max2(localdom$col_lower_[col + 1L], fixval)
  }
  repeat { # retry:
    nbacktracks <- nbacktracks + 1L
    .hmip_HeuristicNeighbourhood_backtracked(neighbourhood)
    if (.hmip_HighsSearch_getCurrentDepth(heur) > targetdepth) {
      if (!.hmip_HighsSearch_backtrackUntilDepth(heur, targetdepth)) {
        self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
        return(invisible())
      }
    }
    repeat {
      .hmip_HighsSearch_evaluateNode(heur)
      if (.hmip_HighsSearch_currentNodePruned(heur)) {
        nbacktracks <- nbacktracks + 1L
        if (.hmip_HighsDomain_infeasible(md$domain)) {
          self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
          return(invisible())
        }
        if (!.hmip_HighsSearch_backtrack(heur)) break
        .hmip_HeuristicNeighbourhood_backtracked(neighbourhood)
        next
      }
      fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
      if (fixingrate >= maxfixingrate) break
      if (stop) break
      if (nbacktracks >= 10L) break
      numBranched <- 0L
      stopFixingRate <- .hmip_min2(1.0 - (1.0 - .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)) * 0.9,
                                   maxfixingrate)
      relaxationsol <- .hmip_HighsLpRelaxation_getSolution(heurlp)$col_value
      for (i in self$intcols) {
        if (localdom$col_lower_[i + 1L] == localdom$col_upper_[i + 1L]) next
        downval <- floor(relaxationsol[i + 1L] + feastol)
        upval <- ceiling(relaxationsol[i + 1L] - feastol)
        downval <- .hmip_min2(downval, localdom$col_upper_[i + 1L])
        upval <- .hmip_max2(upval, localdom$col_lower_[i + 1L])
        if (localdom$col_lower_[i + 1L] < downval) {
          numBranched <- numBranched + 1L
          if (.hmip_heur_branch(md, heur, localdom, TRUE, i, downval, downval - 0.5)) break
        }
        if (localdom$col_upper_[i + 1L] > upval) {
          numBranched <- numBranched + 1L
          if (.hmip_heur_branch(md, heur, localdom, FALSE, i, upval, upval + 0.5)) break
        }
        if (.hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood) >= stopFixingRate) break
      }
      if (numBranched == 0L) {
        .hmip_heur_sort_fracints(heurlp, length(heurlp$fractionalints$first), getFixVal)
        fi <- heurlp$fractionalints
        change <- 0.0
        for (t in seq_along(fi$first)) {
          col <- fi$first[t]; val <- fi$second[t]
          fixval <- getFixVal(col, val)
          if (localdom$col_lower_[col + 1L] < fixval) {
            numBranched <- numBranched + 1L
            if (.hmip_heur_branch(md, heur, localdom, TRUE, col, fixval, val)) break
            fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
          }
          if (localdom$col_upper_[col + 1L] > fixval) {
            numBranched <- numBranched + 1L
            if (.hmip_heur_branch(md, heur, localdom, FALSE, col, fixval, val)) break
            fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
          }
          if (fixingrate >= maxfixingrate) break
          change <- change + abs(fixval - val)
          if (change >= 0.5) break
        }
      }
      if (numBranched == 0L) break
      .hmip_HighsLpRelaxation_flushDomain(heurlp, localdom)
    }
    r <- .hmip_heur_dive_finish(self, heur, heurlp, localdom, neighbourhood)
    if (r$action == "return") return(invisible())
    maxfixingrate <- r$maxfixingrate
    targetdepth <- r$targetdepth
  }
}

.hmip_HighsPrimalHeuristics_RINS <- function(self, relaxationsol) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  if (.hmip_HighsDomain_infeasible(md$domain)) return(invisible())
  if (length(relaxationsol) != mipsolver$model_$num_col_) return(invisible())
  .hmip_heur_drop_fixed_intcols(self)
  s <- .hmip_heur_dive_setup(self)
  heur <- s$heur; localdom <- s$localdom
  heurlp <- .hmip_heur_dive_lp(self, heur, localdom)
  maxfixingrate <- .hmip_HighsPrimalHeuristics_determineTargetFixingRate(self)
  minfixingrate <- 0.25
  fixingrate <- 0.0
  stop <- FALSE
  nbacktracks <- -1L
  targetdepth <- 1L
  neighbourhood <- .hmip_new_HeuristicNeighbourhood(mipsolver, localdom)
  feastol <- md$feastol
  fixtolpsol <- TRUE
  getFixVal <- function(col, fracval) {
    fixval <- if (fixtolpsol) floor(relaxationsol[col + 1L] + 0.5)
              else .hmip_calcFixVal(fracval - md$rootlpsol[col + 1L], fracval,
                                    mipsolver$model_$col_cost_[col + 1L])
    fixval <- .hmip_min2(localdom$col_upper_[col + 1L], fixval)
    .hmip_max2(localdom$col_lower_[col + 1L], fixval)
  }
  repeat { # retry:
    nbacktracks <- nbacktracks + 1L
    .hmip_HeuristicNeighbourhood_backtracked(neighbourhood)
    if (.hmip_HighsSearch_getCurrentDepth(heur) > targetdepth) {
      if (!.hmip_HighsSearch_backtrackUntilDepth(heur, targetdepth)) {
        self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
        return(invisible())
      }
    }
    repeat {
      .hmip_HighsSearch_evaluateNode(heur)
      if (.hmip_HighsSearch_currentNodePruned(heur)) {
        nbacktracks <- nbacktracks + 1L
        if (.hmip_HighsDomain_infeasible(md$domain)) {
          self$lp_iterations <- self$lp_iterations + .hmip_HighsSearch_getLocalLpIterations(heur)
          return(invisible())
        }
        if (!.hmip_HighsSearch_backtrack(heur)) break
        .hmip_HeuristicNeighbourhood_backtracked(neighbourhood)
        next
      }
      fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
      if (stop) break
      if (fixingrate >= maxfixingrate) break
      if (nbacktracks >= 10L) break
      # std::partition: the fractional integers whose LP value agrees with
      # the incumbent go first
      fi <- heurlp$fractionalints
      part <- .hmip_std_partition(seq_along(fi$first), function(t) {
        abs(relaxationsol[fi$first[t] + 1L] - md$incumbent[fi$first[t] + 1L]) <= feastol
      })
      ord <- unlist(part$v)
      heurlp$fractionalints <- list(first = fi$first[ord], second = fi$second[ord])
      fixcandend <- part$k
      fixtolpsol <- TRUE
      numBranched <- 0L
      if (fixcandend == 0L) {
        fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
        stopFixingRate <- .hmip_min2(maxfixingrate, 1.0 - (1.0 - fixingrate) * 0.9)
        currlpsol <- .hmip_HighsLpRelaxation_getSolution(heurlp)$col_value
        for (i in self$intcols) {
          if (localdom$col_lower_[i + 1L] == localdom$col_upper_[i + 1L]) next
          if (abs(currlpsol[i + 1L] - md$incumbent[i + 1L]) <= feastol) {
            fixval <- .hmip_nearestInteger(currlpsol[i + 1L])
            if (localdom$col_lower_[i + 1L] < fixval) {
              numBranched <- numBranched + 1L
              if (.hmip_heur_branch(md, heur, localdom, TRUE, i, fixval, fixval - 0.5)) break
              fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
            }
            if (localdom$col_upper_[i + 1L] > fixval) {
              numBranched <- numBranched + 1L
              if (.hmip_heur_branch(md, heur, localdom, FALSE, i, fixval, fixval + 0.5)) break
              fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
            }
            if (fixingrate >= stopFixingRate) break
          }
        }
        if (numBranched != 0L) {
          .hmip_HighsLpRelaxation_flushDomain(heurlp, localdom)
          next
        }
        if (fixingrate >= minfixingrate) break
        fixcandend <- length(heurlp$fractionalints$first)
        fixtolpsol <- FALSE
      }
      .hmip_heur_sort_fracints(heurlp, fixcandend, getFixVal)
      fi <- heurlp$fractionalints
      change <- 0.0
      for (t in seq_len(fixcandend)) {
        col <- fi$first[t]; val <- fi$second[t]
        fixval <- getFixVal(col, val)
        # RINS branches here without propagating, as the C++ does
        if (localdom$col_lower_[col + 1L] < fixval) {
          numBranched <- numBranched + 1L
          if (.hmip_heur_branch(md, heur, localdom, TRUE, col, fixval, val, FALSE)) break
          fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
        }
        if (localdom$col_upper_[col + 1L] > fixval) {
          numBranched <- numBranched + 1L
          if (.hmip_heur_branch(md, heur, localdom, FALSE, col, fixval, val, FALSE)) break
          fixingrate <- .hmip_HeuristicNeighbourhood_getFixingRate(neighbourhood)
        }
        if (fixingrate >= maxfixingrate) break
        change <- change + abs(fixval - val)
        if (change >= 0.5) break
      }
      if (numBranched == 0L) break
      .hmip_HighsLpRelaxation_flushDomain(heurlp, localdom)
    }
    r <- .hmip_heur_dive_finish(self, heur, heurlp, localdom, neighbourhood)
    if (r$action == "return") return(invisible())
    maxfixingrate <- r$maxfixingrate
    targetdepth <- r$targetdepth
  }
}

# fix the integer columns of a copy of the global domain to value(col), in
# intcols order, propagating after each (value NULL skips the column);
# NULL when the domain becomes infeasible
.hmip_heur_fix_and_propagate <- function(md, intcols, value) {
  localdom <- .hmip_HighsDomain_copy(md$domain)
  for (col in intcols) {
    intval <- value(col)
    if (is.null(intval)) next
    intval <- .hmip_min2(localdom$col_upper_[col + 1L], intval)
    intval <- .hmip_max2(localdom$col_lower_[col + 1L], intval)
    .hmip_HighsDomain_fixCol(localdom, col, intval, .hmip_Reason_branching())
    if (.hmip_HighsDomain_infeasible(localdom)) {
      .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
      return(NULL)
    }
    .hmip_HighsDomain_propagate(localdom)
    if (.hmip_HighsDomain_infeasible(localdom)) {
      .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
      return(NULL)
    }
  }
  localdom
}

# the LP over the fixed domain that tryRoundedPoint and randomizedRounding
# solve, with the dual proof conflict when it is infeasible; returns
# list(lprelax, st)
.hmip_heur_fixed_lp <- function(self, localdom, numintcols) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  ncol <- mipsolver$model_$num_col_
  lprelax <- .hmip_new_HighsLpRelaxation(mipsolver)
  .hmip_HighsLpRelaxation_loadModel(lprelax)
  .hmip_HighsLpRelaxation_setIterationLimit(lprelax, max(10000, 2 * md$firstrootlpiters))
  lps <- .hmip_HighsLpRelaxation_getLpSolver(lprelax)
  .hmip_Highs_changeColsBoundsRange(lps, 0L, ncol - 1L, localdom$col_lower_, localdom$col_upper_)
  rpo <- isTRUE(mipsolver$options_mip_$mip_root_presolve_only)
  if (rpo) .hmip_Highs_setOptionValue(lps, "presolve", "off")
  if (!rpo && (5 * numintcols) %/% ncol >= 1) .hmip_Highs_setOptionValue(lps, "presolve", "on")
  else .hmip_Highs_setBasis(lps, md$firstrootbasis)
  st <- .hmip_HighsLpRelaxation_resolveLp(lprelax)
  if (st == "Infeasible") {
    proof <- .hmip_HighsLpRelaxation_computeDualInfProof(lprelax, md$domain)
    if (isTRUE(proof$ok)) {
      cg <- .hmip_new_HighsCutGeneration(lprelax, md$cutpool)
      .hmip_HighsCutGeneration_generateConflict(cg, localdom, proof$inds, proof$vals, proof$rhs)
    }
  }
  list(lprelax = lprelax, st = st)
}

.hmip_HighsPrimalHeuristics_tryRoundedPoint <- function(self, point, solution_source) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  integerFeasible <- TRUE
  numintcols <- length(self$intcols)
  localdom <- .hmip_heur_fix_and_propagate(md, self$intcols, function(col) {
    intval <- point[col + 1L]
    rounded <- .hmip_round(intval)
    feasible <- abs(intval - rounded) <= md$feastol
    integerFeasible <<- integerFeasible && feasible
    if (feasible) rounded else NULL
  })
  if (is.null(localdom)) return(FALSE)
  if (numintcols != mipsolver$model_$num_col_) {
    r <- .hmip_heur_fixed_lp(self, localdom, numintcols)
    if (r$st == "Infeasible") return(FALSE)
    if (.hmip_HighsLpRelaxation_unscaledPrimalFeasible(r$st)) {
      lpsol <- .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(r$lprelax))$col_value
      if (!integerFeasible) {
        .hmip_HighsPrimalHeuristics_ziRound(self, lpsol)
        return(.hmip_HighsMipSolverData_trySolution(md, lpsol, solution_source))
      }
      .hmip_HighsMipSolverData_addIncumbent(md, lpsol, .hmip_HighsLpRelaxation_getObjective(r$lprelax),
                                            solution_source)
      return(TRUE)
    }
  }
  .hmip_HighsMipSolverData_trySolution(md, localdom$col_lower_, solution_source)
}

.hmip_HighsPrimalHeuristics_linesearchRounding <- function(self, point1, point2, solution_source) {
  md <- self$mipsolver$mipdata_
  feastol <- md$feastol
  roundedpoint <- numeric(self$mipsolver$model_$num_col_)
  alpha <- 0.0
  while (alpha < 1.0) {
    nextalpha <- 1.0
    reachedpoint2 <- TRUE
    for (col in self$intcols) {
      j <- col + 1L
      if (md$uplocks[j] == 0) {
        roundedpoint[j] <- ceiling(.hmip_max2(point1[j], point2[j]) - feastol)
        next
      }
      if (md$downlocks[j] == 0) {
        roundedpoint[j] <- floor(.hmip_min2(point1[j], point2[j]) + feastol)
        next
      }
      convexcomb <- (1.0 - alpha) * point1[j] + alpha * point2[j]
      intpoint2 <- floor(point2[j] + 0.5)
      roundedpoint[j] <- floor(convexcomb + 0.5)
      if (roundedpoint[j] == intpoint2) next
      reachedpoint2 <- FALSE
      tmpalpha <- (roundedpoint[j] + 0.5 + feastol - point1[j]) / abs(point2[j] - point1[j])
      if (tmpalpha < nextalpha && tmpalpha > alpha + 1e-2) nextalpha <- tmpalpha
    }
    if (.hmip_HighsPrimalHeuristics_tryRoundedPoint(self, roundedpoint, solution_source)) return(TRUE)
    if (reachedpoint2) return(FALSE)
    alpha <- nextalpha
  }
  FALSE
}

.hmip_HighsPrimalHeuristics_randomizedRounding <- function(self, relaxationsol) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  ncol <- mipsolver$model_$num_col_
  if (length(relaxationsol) != ncol) return(invisible())
  localdom <- .hmip_heur_fix_and_propagate(md, self$intcols, function(i) {
    if (md$uplocks[i + 1L] == 0) ceiling(relaxationsol[i + 1L] - md$feastol)
    else if (md$downlocks[i + 1L] == 0) floor(relaxationsol[i + 1L] + md$feastol)
    else floor(relaxationsol[i + 1L] + .hs_rng_real(self$randgen, 0.1, 0.9))
  })
  if (is.null(localdom)) return(invisible())
  if (length(md$integer_cols) != ncol) {
    r <- .hmip_heur_fixed_lp(self, localdom, length(self$intcols))
    if (r$st != "Infeasible" && .hmip_HighsLpRelaxation_unscaledPrimalFeasible(r$st))
      .hmip_HighsMipSolverData_addIncumbent(
        md, .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(r$lprelax))$col_value,
        .hmip_HighsLpRelaxation_getObjective(r$lprelax), .hmip_src$RandomizedRounding)
  } else {
    .hmip_HighsMipSolverData_trySolution(md, localdom$col_lower_, .hmip_src$RandomizedRounding)
  }
  invisible()
}

.hmip_HighsPrimalHeuristics_shifting <- function(self, relaxationsol) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  if (length(relaxationsol) != mipsolver$model_$num_col_) return(invisible())
  feastol <- md$feastol
  cur <- relaxationsol
  t <- 0L
  lp <- mipsolver$model_
  # the C++ copies the LP relaxation (HighsLpRelaxation lprelax(mipdata_->lp))
  # only to read its fractional integers
  fi <- .hmip_HighsLpRelaxation_getFractionalIntegers(md$lp)
  fcols <- fi$first; fvals <- fi$second
  rows <- .hmip_HighsMipSolverData_getInfeasibleRows(md, cur)
  previous_size <- length(rows)
  hasInfeasibleConstraints <- length(rows) != 0L
  iterationsWithoutReductions <- 0L
  maxIterationsWithoutReductions <- 5L
  shifts <- new.env(parent = emptyenv())
  isMaximization <- mipsolver$orig_model_$sense_ == -1
  HMAX <- .Machine$integer.max
  while ((length(fcols) > 0L || hasInfeasibleConstraints) &&
         iterationsWithoutReductions <= maxIterationsWithoutReductions &&
         t <= length(md$integer_cols)) {
    t <- t + 1L
    fractionalIntegersReduced <- FALSE
    iterationsWithoutReductions <- iterationsWithoutReductions + 1L
    if (hasInfeasibleConstraints) {
      fractionalIntegerFound <- FALSE
      rIndex <- 0L
      while (!fractionalIntegerFound && rIndex != length(rows)) {
        r <- rows[[rIndex + 1L]][[1L]]
        st <- md$ARstart_[r + 1L]; en <- md$ARstart_[r + 2L]
        if (en > st) for (jInd in (st + 1L):en) {
          fractionalIntegerFound <- any(fcols == md$ARindex_[jInd])
          if (fractionalIntegerFound) break
        }
        rIndex <- rIndex + 1L
      }
      row_index <- rIndex - 1L
      if (!fractionalIntegerFound) row_index <- .hs_rng_integer(self$randgen, length(rows))
      rr <- rows[[row_index + 1L]]
      row <- rr[[1L]]; row_sense <- rr[[2L]]; infeasibility <- rr[[3L]]
      score_min <- Inf; j_min <- HMAX; x_j_min <- Inf; aij_min <- 0.0; moveValueUp <- FALSE
      repair <- function(col, direction, row_sense, numLocks, coef, cost, isInteger, isAtBound) {
        if ((row_sense < 0 || direction * coef > 0) && (row_sense > 0 || direction * coef < 0)) return()
        if (isAtBound) return()
        score <- Inf
        if (any(fcols == col)) {
          score <- -1.0 + 1.0 / (numLocks + 1)
        } else {
          sh <- shifts[[as.character(col)]]
          if (is.null(sh)) {
            score <- direction * (if (isMaximization) -cost else cost)
          } else {
            score <- 0.0
            # pow(1.1, k) of libm; R's ^ calls pow for these arguments
            for (s in sh) if (direction * s > 0) score <- score + 1.1^(direction * s - t)
          }
          if (isInteger) score <- score + 1
        }
        if (score < score_min) {
          score_min <<- score; j_min <<- col; aij_min <<- coef
          x_j_min <<- cur[col + 1L]; moveValueUp <<- direction > 0
        }
      }
      st <- md$ARstart_[row + 1L]; en <- md$ARstart_[row + 2L]
      if (en > st) for (jInd in (st + 1L):en) {
        j <- md$ARindex_[jInd]
        if (lp$col_lower_[j + 1L] == lp$col_upper_[j + 1L]) next
        isInt <- lp$integrality_[j + 1L] == 1L
        repair(j, 1L, row_sense, md$downlocks[j + 1L], md$ARvalue_[jInd], lp$col_cost_[j + 1L], isInt,
               abs(lp$col_upper_[j + 1L] - cur[j + 1L]) <= feastol)
        repair(j, -1L, row_sense, md$uplocks[j + 1L], md$ARvalue_[jInd], lp$col_cost_[j + 1L], isInt,
               abs(cur[j + 1L] - lp$col_lower_[j + 1L]) <= feastol)
      }
      if (j_min != HMAX) {
        p <- match(j_min, fcols)
        if (!is.na(p)) {
          fcols <- fcols[-p]; fvals <- fvals[-p]
          fractionalIntegersReduced <- TRUE
        }
        key <- as.character(j_min)
        isInt <- lp$integrality_[j_min + 1L] == 1L
        if (moveValueUp) {
          if (fractionalIntegersReduced) {
            cur[j_min + 1L] <- ceiling(x_j_min - feastol)
          } else {
            cur[j_min + 1L] <- if (isInt) x_j_min + 1.0
              else .hmip_min2(x_j_min + infeasibility / abs(aij_min), lp$col_upper_[j_min + 1L] + feastol)
            shifts[[key]] <- c(shifts[[key]], t)
          }
        } else {
          if (fractionalIntegersReduced) {
            cur[j_min + 1L] <- floor(x_j_min + feastol)
          } else {
            cur[j_min + 1L] <- if (isInt) x_j_min - 1.0
              else .hmip_max2(x_j_min - infeasibility / abs(aij_min), lp$col_lower_[j_min + 1L] - feastol)
            shifts[[key]] <- c(shifts[[key]], -t)
          }
        }
      }
    } else {
      xi_max <- -1; delta_c_min <- Inf; pind_j_min <- HMAX; j_min <- HMAX; x_j_min <- Inf; sigma <- 0L
      for (i in seq_along(fcols)) {
        col <- fcols[i]; v <- fvals[i]
        isBetter <- function(xi, roundedval, direction) {
          c_min <- lp$col_cost_[col + 1L] * (roundedval - v)
          if (xi > xi_max || (xi == xi_max && c_min < delta_c_min)) {
            xi_max <<- xi; delta_c_min <<- c_min; pind_j_min <<- i - 1L
            j_min <<- col; x_j_min <<- roundedval; sigma <<- direction
          }
        }
        isBetter(md$uplocks[col + 1L], floor(v + feastol), -1L)
        isBetter(md$downlocks[col + 1L], ceiling(v - feastol), 1L)
      }
      if (sigma != 0L) cur[j_min + 1L] <- x_j_min
      if (pind_j_min != HMAX) {
        fcols <- fcols[-(pind_j_min + 1L)]; fvals <- fvals[-(pind_j_min + 1L)]
        fractionalIntegersReduced <- TRUE
      }
    }
    rows <- .hmip_HighsMipSolverData_getInfeasibleRows(md, cur)
    hasInfeasibleConstraints <- length(rows) != 0L
    if (length(rows) < previous_size || fractionalIntegersReduced) iterationsWithoutReductions <- 0L
    previous_size <- length(rows)
  }
  if (hasInfeasibleConstraints) {
    .hmip_HighsPrimalHeuristics_tryRoundedPoint(self, cur, .hmip_src$Shifting)
  } else if (length(fcols) > 0L) {
    .hmip_HighsPrimalHeuristics_ziRound(self, cur)
  } else {
    .hmip_HighsMipSolverData_trySolution(md, cur, .hmip_src$Shifting)
  }
  invisible()
}

.hmip_HighsPrimalHeuristics_ziRound <- function(self, relaxationsol) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  if (length(relaxationsol) != mipsolver$model_$num_col_) return(invisible())
  feastol <- md$feastol
  cur <- relaxationsol
  zi <- function(x) .hmip_min2(ceiling(x - feastol) - x, x - floor(x + feastol))
  zi_total <- c(0, 0)
  for (i in self$intcols) zi_total <- .cd_add(zi_total, zi(cur[i + 1L]))
  if (.cd_val(zi_total) <= feastol) return(invisible())
  lp <- mipsolver$model_
  nrow <- lp$num_row_
  st <- lp$a_matrix_$start_; ix <- lp$a_matrix_$index_; vx <- lp$a_matrix_$value_
  XrowLower <- lp$row_lower_; XrowUpper <- lp$row_upper_
  loop_count <- 0L
  max_loop_count <- 5L
  improvement <- c(Inf, 0)
  while (.cd_val(zi_total) > feastol && .cd_val(improvement) > feastol && loop_count <= max_loop_count) {
    previous_zi_total <- zi_total
    loop_count <- loop_count + 1L
    for (j in self$intcols) {
      x <- cur[j + 1L]
      if (.hmip_fractionality(x) <= feastol) next
      # calculateRowValuesQuad
      rq_hi <- numeric(nrow); rq_lo <- numeric(nrow)
      for (col in seq_len(lp$num_col_)) {
        if (st[col + 1L] > st[col]) for (el in (st[col] + 1L):st[col + 1L]) {
          r <- ix[el] + 1L
          q <- .cd_add(c(rq_hi[r], rq_lo[r]), cur[col] * vx[el])
          rq_hi[r] <- q[1L]; rq_lo[r] <- q[2L]
        }
      }
      rowActivities <- rq_hi + rq_lo
      min_up <- Inf; min_lo <- Inf
      if (st[j + 2L] > st[j + 1L]) for (el in (st[j + 1L] + 1L):st[j + 2L]) {
        i <- ix[el] + 1L
        aij <- vx[el]
        slack_upper <- XrowUpper[i] - rowActivities[i]
        slack_lower <- rowActivities[i] - XrowLower[i]
        min_up <- .hmip_min2(min_up, (if (aij > 0) slack_upper else -slack_lower) / aij)
        min_lo <- .hmip_min2(min_lo, (if (aij > 0) slack_lower else -slack_upper) / aij)
      }
      upper_bound <- .hmip_min2(lp$col_upper_[j + 1L] - x, min_up)
      lower_bound <- .hmip_min2(x - lp$col_lower_[j + 1L], min_lo)
      change <- NULL
      zu <- zi(x + upper_bound); zl <- zi(x - lower_bound); z0 <- zi(x)
      if (abs(zu - zl) <= feastol && zu < z0) {
        cost <- lp$col_cost_[j + 1L]
        ubObjChangeSmaller <- cost * (x + upper_bound) <= cost * (x - lower_bound)
        isMinimization <- lp$sense_ == 1
        change <- if ((isMinimization && ubObjChangeSmaller) || (!isMinimization && !ubObjChangeSmaller))
          upper_bound else -lower_bound
      } else if (zu < zl && zu < z0) {
        change <- upper_bound
      } else if (zu > zl && zl < z0) {
        change <- -lower_bound
      }
      if (!is.null(change)) {
        old <- cur[j + 1L]
        cur[j + 1L] <- cur[j + 1L] + change
        zi_total <- .cd_add(.cd_sub(zi_total, zi(old)), zi(cur[j + 1L]))
      }
    }
    improvement <- .cd_sub_cd(previous_zi_total, zi_total)
  }
  .hmip_HighsMipSolverData_trySolution(md, cur, .hmip_src$ZiRound)
  invisible()
}

# HighsPrimalHeuristics::feasibilityPump (HighsPrimalHeuristics.cpp:1476);
# the unordered_set of reference points only answers membership, so an
# environment keyed by the joined values stands in for it
.hmip_HighsPrimalHeuristics_feasibilityPump <- function(self) {
  mipsolver <- self$mipsolver
  md <- mipsolver$mipdata_
  ncol <- mipsolver$model_$num_col_
  lprelax <- .hmip_HighsLpRelaxation_copy(md$lp)
  referencepoints <- new.env(parent = emptyenv())
  emplace <- function(rp) {
    key <- paste(rp, collapse = ",")
    if (exists(key, envir = referencepoints, inherits = FALSE)) return(FALSE)
    assign(key, TRUE, envir = referencepoints)
    TRUE
  }
  st <- .hmip_HighsLpRelaxation_resolveLp(lprelax)
  self$lp_iterations <- self$lp_iterations + .hmip_HighsLpRelaxation_getNumLpIterations(lprelax)
  mask <- rep(1L, ncol)
  cost <- numeric(ncol)
  lps <- .hmip_HighsLpRelaxation_getLpSolver(lprelax)
  .hmip_Highs_setOptionValue(lps, "simplex_strategy", 4L)  # kSimplexStrategyPrimal
  .hmip_HighsLpRelaxation_setObjectiveLimit(lprelax)
  .hmip_Highs_setOptionValue(lps, "primal_simplex_bound_perturbation_multiplier", 0.0)
  # setIterationLimit takes a HighsInt: the double converts by truncation
  .hmip_HighsLpRelaxation_setIterationLimit(lprelax, trunc(5 * md$avgrootlpiters))
  intcols <- md$integer_cols
  nint <- length(intcols)
  while (length(.hmip_HighsLpRelaxation_getFractionalIntegers(lprelax)$first) != 0L) {
    lpsol <- .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(lprelax))$col_value
    roundedsol <- lpsol
    referencepoint <- numeric(nint)
    localdom <- .hmip_HighsDomain_copy(md$domain)
    for (k in seq_len(nint)) {
      i <- intcols[k]
      intval <- floor(roundedsol[i + 1L] + .hs_rng_real(self$randgen, 0.4, 0.6))
      intval <- .hmip_max2(intval, localdom$col_lower_[i + 1L])
      intval <- .hmip_min2(intval, localdom$col_upper_[i + 1L])
      roundedsol[i + 1L] <- intval
      referencepoint[k] <- trunc(intval)
      if (!.hmip_HighsDomain_infeasible(localdom)) {
        .hmip_HighsDomain_fixCol(localdom, i, intval, .hmip_Reason_branching())
        if (.hmip_HighsDomain_infeasible(localdom)) {
          .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
          next
        }
        .hmip_HighsDomain_propagate(localdom)
        if (.hmip_HighsDomain_infeasible(localdom)) {
          .hmip_HighsDomain_conflictAnalysis(localdom, md$conflictPool)
          next
        }
      }
    }
    havecycle <- !emplace(referencepoint)
    k <- 0L
    while (havecycle && k < 2L) {
      for (t in 1:10) {
        flippos <- .hs_rng_integer(self$randgen, nint)
        col <- intcols[flippos + 1L]
        if (roundedsol[col + 1L] > lpsol[col + 1L])
          roundedsol[col + 1L] <- floor(lpsol[col + 1L])
        else if (roundedsol[col + 1L] < lpsol[col + 1L])
          roundedsol[col + 1L] <- ceiling(lpsol[col + 1L])
        else if (roundedsol[col + 1L] < md$domain$col_upper_[col + 1L])
          roundedsol[col + 1L] <- md$domain$col_upper_[col + 1L]
        else
          roundedsol[col + 1L] <- md$domain$col_lower_[col + 1L]
        referencepoint[flippos + 1L] <- trunc(roundedsol[col + 1L])
      }
      havecycle <- !emplace(referencepoint)
      k <- k + 1L
    }
    if (havecycle) return(invisible())
    if (.hmip_HighsPrimalHeuristics_linesearchRounding(self, lpsol, roundedsol, .hmip_src$FeasibilityPump))
      return(invisible())
    if (.hmip_HighsLpRelaxation_getNumLpIterations(lprelax) >= 1000 + md$avgrootlpiters * 5) break
    for (i in intcols) {
      if (md$uplocks[i + 1L] == 0 || md$downlocks[i + 1L] == 0)
        cost[i + 1L] <- 0.0
      else if (lpsol[i + 1L] > roundedsol[i + 1L] - md$feastol)
        cost[i + 1L] <- -1.0 + .hs_rng_real(self$randgen, -1e-4, 1e-4)
      else
        cost[i + 1L] <- 1.0 + .hs_rng_real(self$randgen, -1e-4, 1e-4)
    }
    .hmip_Highs_changeColsCost(.hmip_HighsLpRelaxation_getLpSolver(lprelax), mask, cost)
    niters <- -.hmip_HighsLpRelaxation_getNumLpIterations(lprelax)
    st <- .hmip_HighsLpRelaxation_resolveLp(lprelax)
    niters <- niters + .hmip_HighsLpRelaxation_getNumLpIterations(lprelax)
    if (niters == 0) break
    self$lp_iterations <- self$lp_iterations + niters
  }
  if (length(.hmip_HighsLpRelaxation_getFractionalIntegers(lprelax)$first) == 0L &&
      .hmip_HighsLpRelaxation_unscaledPrimalFeasible(st))
    .hmip_HighsMipSolverData_addIncumbent(
      md, .hmip_Highs_getSolution(.hmip_HighsLpRelaxation_getLpSolver(lprelax))$col_value,
      .hmip_HighsLpRelaxation_getObjective(lprelax), .hmip_src$FeasibilityPump)
  invisible()
}

.hmip_HighsPrimalHeuristics_centralRounding <- function(self) {
  md <- self$mipsolver$mipdata_
  if (length(md$analyticCenter) != self$mipsolver$model_$num_col_) return(invisible())
  start <- if (length(md$firstlpsol)) md$firstlpsol
           else if (length(md$rootlpsol)) md$rootlpsol
           else md$analyticCenter
  .hmip_HighsPrimalHeuristics_linesearchRounding(self, start, md$analyticCenter, .hmip_src$CentralRounding)
  invisible()
}

.hmip_HighsPrimalHeuristics_flushStatistics <- function(self) {
  md <- self$mipsolver$mipdata_
  md$total_repair_lp <- md$total_repair_lp + self$total_repair_lp
  md$total_repair_lp_feasible <- md$total_repair_lp_feasible + self$total_repair_lp_feasible
  md$total_repair_lp_iterations <- md$total_repair_lp_iterations + self$total_repair_lp_iterations
  self$total_repair_lp <- 0
  self$total_repair_lp_feasible <- 0
  self$total_repair_lp_iterations <- 0
  md$heuristic_lp_iterations <- md$heuristic_lp_iterations + self$lp_iterations
  md$total_lp_iterations <- md$total_lp_iterations + self$lp_iterations
  self$lp_iterations <- 0
  invisible()
}
