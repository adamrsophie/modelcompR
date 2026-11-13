# COBYLA of NLopt 2.11.0, ported to R from src/algs/cobyla/cobyla.c: the C
# version by Jean-Sebastien Roy of Powell's COBYLA2, with the changes of
# Steven G. Johnson (explicit bound constraints, pseudorandom simplex steps
# from a linear congruential generator, the growth of RHO after good
# predictions, the stopping criteria of NLopt). This is a different code
# from the COBYLA of PRIMA (R/optimization-prima-cobyla.R).
#
# The NLopt library that NLopt.jl loads (NLopt_jll) rounds every product
# and sum separately (no fused multiply add), so plain R arithmetic in the
# order of the C expressions reproduces it. The float constants of the C
# code (1e-6f, .1f, .2f) keep their single precision values.
#
# Copyright (c) 1992, Michael J. D. Powell; Copyright (c) 2004, Jean-Sebastien
# Roy; Copyright (c) 2007-2014 Massachusetts Institute of Technology. Released
# under the MIT license (see R/optimization-nlopt-util.R).

# MIN2 and MAX2 of cobyla.c
.ncob_min <- function(a, b) if (isTRUE(a <= b)) a else b
.ncob_max <- function(a, b) if (isTRUE(a >= b)) a else b

.ncob_f1e6 <- 9.999999974752427e-07 # 1e-6f
.ncob_f01 <- 0.10000000149011612 # .1f
.ncob_f02 <- 0.20000000298023224 # .2f

# a fused dot product (temp += x[i] * y[i])
.ncob_dot <- function(x, y) {
  s <- 0
  for (i in seq_along(x)) s <- (x[i] * y[i] + s)
  s
}

# trstlp: the step dx of the linear programming subproblem (a: n by m + 1,
# b: length m + 1). Returns list(rc, dx, ifull).
.ncob_trstlp <- function(n, m, a, b, rho) {
  ifull <- 1L
  mcon <- m
  nact <- 0L
  resmax <- 0
  icon <- 0L
  nactx <- 0L
  resold <- 0
  kk <- 0L
  z <- diag(1, n)
  dx <- numeric(n)
  iact <- integer(m + 1L)
  vmultc <- numeric(m + 1L)
  vmultd <- numeric(m + 1L)
  zdota <- numeric(n)
  sdirn <- numeric(n)
  dxnew <- numeric(n)
  optold <- 0
  icount <- 0L
  if (m >= 1L) {
    for (k in seq_len(m)) if (b[k] > resmax) { resmax <- b[k]; icon <- k }
    for (k in seq_len(m)) { iact[k] <- k; vmultc[k] <- resmax - b[k] }
  }
  # the rotation of columns k and kp of z that moves the active constraint
  # with gradient a[, kw] from position kp to k (L170, L270 and L210)
  rotate <- function(k, kp, kw) {
    sp <- .ncob_dot(z[, k], a[, kw])
    temp <- sqrt((sp * sp + zdota[kp] * zdota[kp]))
    alpha <- zdota[kp] / temp
    beta <- sp / temp
    zdota[kp] <<- alpha * zdota[k]
    zdota[k] <<- temp
    t <- (alpha * z[, kp] + beta * z[, k])
    z[, kp] <<- (alpha * z[, k] - beta * z[, kp])
    z[, k] <<- t
  }
  lab <- if (resmax == 0) "L480" else "L60"
  repeat {
    if (lab == "L60") {
      optold <- 0
      icount <- 0L
      lab <- "L70"
    }
    if (lab == "L70") {
      if (mcon == m) {
        optnew <- resmax
      } else {
        optnew <- 0
        for (i in seq_len(n)) optnew <- (-dx[i] * a[i, mcon] + optnew)
      }
      if (icount == 0L || optnew < optold) {
        optold <- optnew
        nactx <- nact
        icount <- 3L
      } else if (nact > nactx) {
        nactx <- nact
        icount <- 3L
      } else {
        icount <- icount - 1L
        if (icount == 0L) { lab <- "L490"; next }
      }
      if (icon <= nact) { lab <- "L260"; next }
      kk <- iact[icon]
      dxnew <- a[, kk]
      tot <- 0
      k <- n
      while (k > nact) {
        sp <- 0
        spabs <- 0
        for (i in seq_len(n)) {
          temp <- z[i, k] * dxnew[i]
          sp <- sp + temp
          spabs <- spabs + abs(temp)
        }
        acca <- (abs(sp) * 0.1 + spabs)
        accb <- (abs(sp) * 0.2 + spabs)
        if (spabs >= acca || acca >= accb) sp <- 0
        if (tot == 0) {
          tot <- sp
        } else {
          kp <- k + 1L
          temp <- sqrt((sp * sp + tot * tot))
          alpha <- sp / temp
          beta <- tot / temp
          tot <- temp
          t <- (alpha * z[, k] + beta * z[, kp])
          z[, kp] <- (alpha * z[, kp] - beta * z[, k])
          z[, k] <- t
        }
        k <- k - 1L
      }
      if (tot != 0) {
        nact <- nact + 1L
        zdota[nact] <- tot
        vmultc[icon] <- vmultc[nact]
        vmultc[nact] <- 0
        lab <- "L210"
      } else {
        # a deletion makes room for the new active constraint
        ratio <- -1
        for (k in rev(seq_len(nact))) {
          zdotv <- 0
          zdvabs <- 0
          for (i in seq_len(n)) {
            temp <- z[i, k] * dxnew[i]
            zdotv <- zdotv + temp
            zdvabs <- zdvabs + abs(temp)
          }
          acca <- (abs(zdotv) * 0.1 + zdvabs)
          accb <- (abs(zdotv) * 0.2 + zdvabs)
          if (zdvabs < acca && acca < accb) {
            temp <- zdotv / zdota[k]
            if (temp > 0 && iact[k] <= m) {
              tempa <- vmultc[k] / temp
              if (ratio < 0 || tempa < ratio) ratio <- tempa
            }
            if (k >= 2L) {
              kw <- iact[k]
              dxnew <- (-temp * a[, kw] + dxnew)
            }
            vmultd[k] <- temp
          } else {
            vmultd[k] <- 0
          }
        }
        if (ratio < 0) { lab <- "L490"; next }
        for (k in seq_len(nact)) vmultc[k] <- .ncob_max(0, (-ratio * vmultd[k] + vmultc[k]))
        if (icon < nact) {
          isave <- iact[icon]
          vsave <- vmultc[icon]
          k <- icon
          repeat {
            kp <- k + 1L
            kw <- iact[kp]
            rotate(k, kp, kw)
            iact[k] <- kw
            vmultc[k] <- vmultc[kp]
            k <- kp
            if (k >= nact) break
          }
          iact[k] <- isave
          vmultc[k] <- vsave
        }
        temp <- .ncob_dot(z[, nact], a[, kk])
        if (temp == 0) { lab <- "L490"; next }
        zdota[nact] <- temp
        vmultc[icon] <- 0
        vmultc[nact] <- ratio
        lab <- "L210"
      }
    }
    if (lab == "L210") {
      iact[icon] <- iact[nact]
      iact[nact] <- kk
      if (mcon > m && kk != mcon) {
        k <- nact - 1L
        sp <- .ncob_dot(z[, k], a[, kk])
        temp <- sqrt((sp * sp + zdota[nact] * zdota[nact]))
        alpha <- zdota[nact] / temp
        beta <- sp / temp
        zdota[nact] <- alpha * zdota[k]
        zdota[k] <- temp
        t <- (alpha * z[, nact] + beta * z[, k])
        z[, nact] <- (alpha * z[, k] - beta * z[, nact])
        z[, k] <- t
        iact[nact] <- iact[k]
        iact[k] <- kk
        temp <- vmultc[k]
        vmultc[k] <- vmultc[nact]
        vmultc[nact] <- temp
      }
      if (mcon > m) {
        lab <- "L320"
      } else {
        kk <- iact[nact]
        temp <- .ncob_dot(sdirn, a[, kk])
        temp <- temp + -1
        temp <- temp / zdota[nact]
        sdirn <- (-temp * z[, nact] + sdirn)
        lab <- "L340"
      }
    }
    if (lab == "L260") {
      # delete the constraint with index iact[icon] from the active set
      if (icon < nact) {
        isave <- iact[icon]
        vsave <- vmultc[icon]
        k <- icon
        repeat {
          kp <- k + 1L
          kk <- iact[kp]
          rotate(k, kp, kk)
          iact[k] <- kk
          vmultc[k] <- vmultc[kp]
          k <- kp
          if (k >= nact) break
        }
        iact[k] <- isave
        vmultc[k] <- vsave
      }
      nact <- nact - 1L
      if (mcon > m) {
        lab <- "L320"
      } else {
        temp <- .ncob_dot(sdirn, z[, nact + 1L])
        sdirn <- (-temp * z[, nact + 1L] + sdirn)
        lab <- "L340"
      }
    }
    if (lab == "L320") {
      temp <- 1 / zdota[nact]
      sdirn <- temp * z[, nact]
      lab <- "L340"
    }
    if (lab == "L340") {
      dd <- rho * rho
      sd <- 0
      ss <- 0
      for (i in seq_len(n)) {
        if (abs(dx[i]) >= rho * .ncob_f1e6) dd <- (-dx[i] * dx[i] + dd)
        sd <- (dx[i] * sdirn[i] + sd)
        ss <- (sdirn[i] * sdirn[i] + ss)
      }
      if (dd <= 0) { lab <- "L490"; next }
      temp <- sqrt(ss * dd)
      if (abs(sd) >= temp * .ncob_f1e6) temp <- sqrt((ss * dd + sd * sd))
      stpful <- dd / (temp + sd)
      step <- stpful
      if (mcon == m) {
        acca <- (resmax * 0.1 + step)
        accb <- (resmax * 0.2 + step)
        if (step >= acca || acca >= accb) { lab <- "L480"; next }
        step <- .ncob_min(step, resmax)
      }
      if (.nlopt_isinf(step)) return(list(rc = "ROUNDOFF_LIMITED", dx = dx, ifull = ifull))
      dxnew <- (step * sdirn + dx)
      if (mcon == m) {
        resold <- resmax
        resmax <- 0
        for (k in seq_len(nact)) {
          kk <- iact[k]
          temp <- b[kk]
          for (i in seq_len(n)) temp <- (-a[i, kk] * dxnew[i] + temp)
          resmax <- .ncob_max(resmax, temp)
        }
      }
      k <- nact
      while (k >= 1L) {
        zdotw <- 0
        zdwabs <- 0
        for (i in seq_len(n)) {
          temp <- z[i, k] * dxnew[i]
          zdotw <- zdotw + temp
          zdwabs <- zdwabs + abs(temp)
        }
        acca <- (abs(zdotw) * 0.1 + zdwabs)
        accb <- (abs(zdotw) * 0.2 + zdwabs)
        if (zdwabs >= acca || acca >= accb) zdotw <- 0
        vmultd[k] <- zdotw / zdota[k]
        if (k < 2L) break
        kk <- iact[k]
        dxnew <- (-vmultd[k] * a[, kk] + dxnew)
        k <- k - 1L
      }
      if (mcon > m) vmultd[nact] <- .ncob_max(0, vmultd[nact])
      dxnew <- (step * sdirn + dx)
      if (mcon > nact) {
        for (k in (nact + 1L):mcon) {
          kk <- iact[k]
          sum <- resmax - b[kk]
          sumabs <- resmax + abs(b[kk])
          for (i in seq_len(n)) {
            temp <- a[i, kk] * dxnew[i]
            sum <- sum + temp
            sumabs <- sumabs + abs(temp)
          }
          acca <- (abs(sum) * .ncob_f01 + sumabs)
          accb <- (abs(sum) * .ncob_f02 + sumabs)
          if (sumabs >= acca || acca >= accb) sum <- 0
          vmultd[k] <- sum
        }
      }
      ratio <- 1
      icon <- 0L
      for (k in seq_len(mcon)) {
        if (vmultd[k] < 0) {
          temp <- vmultc[k] / (vmultc[k] - vmultd[k])
          if (temp < ratio) { ratio <- temp; icon <- k }
        }
      }
      temp <- 1 - ratio
      dx <- (temp * dx + ratio * dxnew)
      for (k in seq_len(mcon)) vmultc[k] <- .ncob_max(0, (temp * vmultc[k] + ratio * vmultd[k]))
      if (mcon == m) resmax <- (ratio * (resmax - resold) + resold)
      if (icon > 0L) { lab <- "L70"; next }
      if (step == stpful) return(list(rc = "SUCCESS", dx = dx, ifull = ifull))
      lab <- "L480"
    }
    if (lab == "L480") {
      mcon <- m + 1L
      icon <- mcon
      iact[mcon] <- mcon
      vmultc[mcon] <- 0
      lab <- "L60"
      next
    }
    if (lab == "L490") {
      if (mcon == m) { lab <- "L480"; next }
      ifull <- 0L
      return(list(rc = "SUCCESS", dx = dx, ifull = ifull))
    }
  }
}

# cobylb: the main loop of COBYLA on the scaled variables. `calcfc(x)`
# returns list(f, con) (con: the m constraints, nonnegative when
# satisfied) or NULL after a forced stop; con_tol are the tolerances of the
# feasibility test of stopval.
.ncob_cobylb <- function(n, m, x, rhobeg, rhoend, stop, lb, ub, calcfc, con_tol) {
  mpp <- m + 2L
  np <- n + 1L
  mp <- m + 1L
  alpha <- 0.25; beta <- 2.1; gamma_ <- 0.5; delta <- 1.1
  rho <- rhobeg
  parmu <- 0
  minf <- Inf
  rc <- "SUCCESS"
  seed <- (n + m) %% 4294967296 # the deterministic seed of the LCG
  lcg_urand <- function(a, b) {
    seed <<- (.nl32_mul(seed, 1103515245) + 12345) %% 4294967296
    a + seed * (b - a) / 4294967295
  }
  sim <- matrix(0, n, np)
  simi <- matrix(0, n, n)
  datmat <- matrix(0, mpp, np)
  A <- matrix(0, n, mp)
  con <- numeric(mpp)
  vsig <- numeric(n); veta <- numeric(n); sigbar <- numeric(n); dx <- numeric(n)
  f <- 0; resmax <- 0; ifull <- 0L; iflag <- 0L; prerem <- 0; prerec <- 0; parsig <- 0
  for (i in seq_len(n)) {
    sim[i, np] <- x[i]
    rhocur <- rho
    if (x[i] + rhocur > ub[i]) {
      if (x[i] - rhocur >= lb[i]) rhocur <- -rhocur
      else if (ub[i] - x[i] > x[i] - lb[i]) rhocur <- 0.5 * (ub[i] - x[i])
      else rhocur <- 0.5 * (x[i] - lb[i])
    }
    sim[i, i] <- rhocur
    simi[i, i] <- 1 / rhocur
  }
  jdrop <- np
  ibrnch <- 0L
  # the revision of the simplex with the step dx replacing vertex jdrop
  revise <- function() {
    temp <- 0
    for (i in seq_len(n)) {
      sim[i, jdrop] <<- dx[i]
      temp <- (simi[jdrop, i] * dx[i] + temp)
    }
    simi[jdrop, ] <<- simi[jdrop, ] / temp
    for (j in seq_len(n)) if (j != jdrop) {
      temp <- .ncob_dot(simi[j, ], dx)
      simi[j, ] <<- (-temp * simi[jdrop, ] + simi[j, ])
    }
  }
  lab <- "L40"
  repeat {
    if (lab == "L40") {
      if (.nlopt_stop_forced(stop)) rc <- "FORCED_STOP"
      else if (stop$nevals > 0L) {
        if (.nlopt_stop_evals(stop)) rc <- "MAXEVAL_REACHED"
        else if (.nlopt_stop_time(stop)) rc <- "MAXTIME_REACHED"
      }
      if (rc != "SUCCESS") { lab <- "L600"; next }
      stop$nevals <- stop$nevals + 1L
      r <- calcfc(x)
      if (is.null(r)) { rc <- "FORCED_STOP"; lab <- "L600"; next }
      f <- r$f
      if (m > 0L) con[seq_len(m)] <- r$con
      resmax <- 0
      feasible <- TRUE
      for (k in seq_len(m)) {
        resmax <- .ncob_max(resmax, -con[k])
        if (-con[k] > con_tol[k]) feasible <- FALSE
      }
      if (f < stop$minf_max && feasible) { rc <- "STOPVAL_REACHED"; lab <- "L620"; next }
      con[mp] <- f
      con[mpp] <- resmax
      if (ibrnch == 1L) { lab <- "L440"; next }
      datmat[, jdrop] <- con
      if (stop$nevals <= np) {
        if (jdrop <= n) {
          if (datmat[mp, np] <= f) {
            x[jdrop] <- sim[jdrop, np]
          } else {
            rhocur <- x[jdrop] - sim[jdrop, np]
            sim[jdrop, np] <- x[jdrop]
            for (k in seq_len(mpp)) {
              datmat[k, jdrop] <- datmat[k, np]
              datmat[k, np] <- con[k]
            }
            for (k in seq_len(jdrop)) {
              sim[jdrop, k] <- -rhocur
              temp <- 0
              for (i in k:jdrop) temp <- temp - simi[i, k]
              simi[jdrop, k] <- temp
            }
          }
        }
        if (stop$nevals <= n) {
          jdrop <- stop$nevals
          x[jdrop] <- x[jdrop] + sim[jdrop, jdrop]
          next
        }
      }
      ibrnch <- 1L
      lab <- "L140"
    }
    if (lab == "L140") {
      # the optimal vertex of the current simplex
      phimin <- (parmu * datmat[mpp, np] + datmat[mp, np])
      nbest <- np
      for (j in seq_len(n)) {
        temp <- (parmu * datmat[mpp, j] + datmat[mp, j])
        if (temp < phimin) {
          nbest <- j
          phimin <- temp
        } else if (temp == phimin && parmu == 0) {
          if (datmat[mpp, j] < datmat[mpp, nbest]) nbest <- j
        }
      }
      if (nbest <= n) {
        t <- datmat[, np]; datmat[, np] <- datmat[, nbest]; datmat[, nbest] <- t
        for (i in seq_len(n)) {
          temp <- sim[i, nbest]
          sim[i, nbest] <- 0
          sim[i, np] <- sim[i, np] + temp
          tempa <- 0
          for (k in seq_len(n)) {
            sim[i, k] <- sim[i, k] - temp
            tempa <- tempa - simi[k, i]
          }
          simi[nbest, i] <- tempa
        }
      }
      # an error return if simi is a poor inverse of sim
      error <- 0
      for (i in seq_len(n)) for (j in seq_len(n)) {
        temp <- if (i == j) -1 else 0
        for (k in seq_len(n)) if (sim[k, j] != 0) temp <- (simi[i, k] * sim[k, j] + temp)
        error <- .ncob_max(error, abs(temp))
      }
      if (error > 0.1) { rc <- "ROUNDOFF_LIMITED"; lab <- "L600"; next }
      # the linear approximations of the objective and constraints
      for (k in seq_len(mp)) {
        con[k] <- -datmat[k, np]
        w <- datmat[k, seq_len(n)] + con[k]
        for (i in seq_len(n)) {
          temp <- .ncob_dot(w, simi[, i])
          if (k == mp) temp <- -temp
          A[i, k] <- temp
        }
      }
      iflag <- 1L
      parsig <- alpha * rho
      pareta <- beta * rho
      for (j in seq_len(n)) {
        wsig <- 0
        weta <- 0
        for (i in seq_len(n)) {
          wsig <- (simi[j, i] * simi[j, i] + wsig)
          weta <- (sim[i, j] * sim[i, j] + weta)
        }
        vsig[j] <- 1 / sqrt(wsig)
        veta[j] <- sqrt(weta)
        if (vsig[j] < parsig || veta[j] > pareta) iflag <- 0L
      }
      if (ibrnch == 1L || iflag == 1L) { lab <- "L370"; next }
      # a new vertex to improve the acceptability of the simplex
      jdrop <- 0L
      temp <- pareta
      for (j in seq_len(n)) if (veta[j] > temp) { jdrop <- j; temp <- veta[j] }
      if (jdrop == 0L) for (j in seq_len(n)) if (vsig[j] < temp) { jdrop <- j; temp <- vsig[j] }
      temp <- gamma_ * rho * vsig[jdrop]
      dx <- temp * simi[jdrop, ]
      cvmaxp <- 0
      cvmaxm <- 0
      sum <- 0
      for (k in seq_len(mp)) {
        sum <- .ncob_dot(A[, k], dx)
        if (k < mp) {
          temp <- datmat[k, np]
          cvmaxp <- .ncob_max(cvmaxp, -sum - temp)
          cvmaxm <- .ncob_max(cvmaxm, sum - temp)
        }
      }
      dxsign <- if (parmu * (cvmaxp - cvmaxm) > sum + sum) -1 else 1
      temp <- 0
      for (i in seq_len(n)) {
        dx[i] <- dxsign * dx[i] * lcg_urand(0.01, 1)
        xi <- sim[i, np]
        repeat {
          if (xi + dx[i] > ub[i]) dx[i] <- -dx[i]
          if (xi + dx[i] < lb[i]) {
            if (xi - dx[i] <= ub[i]) dx[i] <- -dx[i]
            else { dx[i] <- dx[i] * 0.5; next }
          }
          break
        }
        sim[i, jdrop] <- dx[i]
        temp <- (simi[jdrop, i] * dx[i] + temp)
      }
      simi[jdrop, ] <- simi[jdrop, ] / temp
      for (j in seq_len(n)) {
        if (j != jdrop) {
          temp <- .ncob_dot(simi[j, ], dx)
          simi[j, ] <- (-temp * simi[jdrop, ] + simi[j, ])
        }
        x[j] <- sim[j, np] + dx[j]
      }
      lab <- "L40"
      next
    }
    if (lab == "L370") {
      tr <- .ncob_trstlp(n, m, A, con[seq_len(mp)], rho)
      if (tr$rc != "SUCCESS") { rc <- tr$rc; lab <- "L600"; next }
      dx <- tr$dx
      ifull <- tr$ifull
      for (i in seq_len(n)) {
        xi <- sim[i, np]
        if (xi + dx[i] > ub[i]) dx[i] <- ub[i] - xi
        if (xi + dx[i] < lb[i]) dx[i] <- xi - lb[i]
      }
      if (ifull == 0L) {
        temp <- .ncob_dot(dx, dx)
        if (temp < rho * 0.25 * rho) { ibrnch <- 1L; lab <- "L550"; next }
      }
      # the predicted changes of f and the constraint violation
      resnew <- 0
      con[mp] <- 0
      sum <- 0
      for (k in seq_len(mp)) {
        sum <- con[k]
        for (i in seq_len(n)) sum <- (-A[i, k] * dx[i] + sum)
        if (k < mp) resnew <- .ncob_max(resnew, sum)
      }
      barmu <- 0
      prerec <- datmat[mpp, np] - resnew
      if (prerec > 0) barmu <- sum / prerec
      if (parmu < barmu * 1.5) {
        parmu <- barmu * 2
        phi <- (parmu * datmat[mpp, np] + datmat[mp, np])
        back <- FALSE
        for (j in seq_len(n)) {
          temp <- (parmu * datmat[mpp, j] + datmat[mp, j])
          if (temp < phi) { back <- TRUE; break }
          if (temp == phi && parmu == 0 && datmat[mpp, j] < datmat[mpp, np]) { back <- TRUE; break }
        }
        if (back) { lab <- "L140"; next }
      }
      prerem <- (parmu * prerec + -sum)
      x <- sim[, np] + dx
      ibrnch <- 1L
      lab <- "L40"
      next
    }
    if (lab == "L440") {
      vmold <- (parmu * datmat[mpp, np] + datmat[mp, np])
      vmnew <- (parmu * resmax + f)
      trured <- vmold - vmnew
      if (parmu == 0 && f == datmat[mp, np]) {
        prerem <- prerec
        trured <- datmat[mpp, np] - resmax
      }
      ratio <- if (trured <= 0) 1 else 0
      jdrop <- 0L
      for (j in seq_len(n)) {
        temp <- abs(.ncob_dot(simi[j, ], dx))
        if (temp > ratio) { jdrop <- j; ratio <- temp }
        sigbar[j] <- temp * vsig[j]
      }
      edgmax <- delta * rho
      l <- 0L
      for (j in seq_len(n)) {
        if (sigbar[j] >= parsig || sigbar[j] >= vsig[j]) {
          temp <- veta[j]
          if (trured > 0) {
            temp <- 0
            for (i in seq_len(n)) {
              d1 <- dx[i] - sim[i, j]
              temp <- (d1 * d1 + temp)
            }
            temp <- sqrt(temp)
          }
          if (temp > edgmax) { l <- j; edgmax <- temp }
        }
      }
      if (l > 0L) jdrop <- l
      if (jdrop == 0L) { lab <- "L550"; next }
      revise()
      datmat[, jdrop] <- con
      if (trured > 0 && trured >= prerem * 0.1) {
        if (trured >= prerem * 0.9 && trured <= prerem * 1.1 && iflag) rho <- rho * 2
        lab <- "L140"
        next
      }
      lab <- "L550"
    }
    if (lab == "L550") {
      if (iflag == 0L) { ibrnch <- 0L; lab <- "L140"; next }
      fbest <- if (ifull == 1L) f else datmat[mp, np]
      if (fbest < minf && .nlopt_stop_ftol(stop, fbest, minf)) { rc <- "FTOL_REACHED"; lab <- "L600"; next }
      minf <- fbest
      if (rho > rhoend) {
        rho <- rho * 0.5
        if (rho <= rhoend * 1.5) rho <- rhoend
        if (parmu > 0) {
          denom <- 0
          cmin <- 0; cmax <- 0
          for (k in seq_len(mp)) {
            cmin <- datmat[k, np]
            cmax <- cmin
            for (i in seq_len(n)) {
              cmin <- .ncob_min(cmin, datmat[k, i])
              cmax <- .ncob_max(cmax, datmat[k, i])
            }
            if (k <= m && cmin < cmax * 0.5) {
              temp <- .ncob_max(cmax, 0) - cmin
              denom <- if (denom <= 0) temp else .ncob_min(denom, temp)
            }
          }
          if (denom == 0) parmu <- 0
          else if (cmax - cmin < parmu * denom) parmu <- (cmax - cmin) / denom
        }
        lab <- "L140"
        next
      }
      rc <- if (rhoend > 0) "XTOL_REACHED" else "ROUNDOFF_LIMITED"
      lab <- if (ifull == 1L) "L620" else "L600"
    }
    if (lab == "L600") {
      x <- sim[, np]
      f <- datmat[mp, np]
      lab <- "L620"
    }
    if (lab == "L620") break
  }
  list(ret = rc, x = x, minf = f)
}

# cobyla_minimize: the NLopt interface (o: the optimizer description with
# its objective and constraints, dx: the initial step)
.nlopt_cobyla_minimize <- function(o, x, stop, dx) {
  n <- length(x)
  lb0 <- o$lb; ub0 <- o$ub
  scale <- .nlopt_compute_rescaling(n, dx)
  if (any(scale == 0 | !is.finite(scale))) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  bnd <- .nlopt_reorder_bounds(lb0 / scale, ub0 / scale)
  lb <- bnd$lb; ub <- bnd$ub
  rhobeg <- abs(dx[1L] / scale[1L])
  rhoend <- stop$xtol_rel * rhobeg
  if (!is.null(stop$xtol_abs)) for (j in seq_len(n)) {
    t <- rep_len(stop$xtol_abs, n)[j] / abs(scale[j])
    if (rhoend < t) rhoend <- t
  }
  fc <- o$fc; h <- o$h
  finite_lb <- !vapply(lb, .nlopt_isinf, TRUE)
  finite_ub <- !vapply(ub, .nlopt_isinf, TRUE)
  m <- .nlopt_count_constraints(fc) + 2L * .nlopt_count_constraints(h) +
    sum(!vapply(lb0, .nlopt_isinf, TRUE)) + sum(!vapply(ub0, .nlopt_isinf, TRUE))
  con_tol <- c(unlist(lapply(fc, function(cc) rep_len(cc$tol, cc$m))),
               unlist(lapply(h, function(cc) rep_len(rep_len(cc$tol, cc$m), 2L * cc$m))))
  con_tol <- c(con_tol, numeric(m - length(con_tol)))
  calcfc <- function(xs) {
    xt <- ifelse(xs < lb, lb, ifelse(xs > ub, ub, xs))
    xt <- xt * scale
    fv <- o$f(xt)
    if (.nlopt_stop_forced(stop)) return(NULL)
    con <- numeric(0)
    for (cc in fc) con <- c(con, -as.numeric(cc$f(xt)))
    for (cc in h) { v <- as.numeric(cc$f(xt)); con <- c(con, v, -v) }
    for (j in seq_len(n)) {
      if (finite_lb[j]) con <- c(con, xs[j] - lb[j])
      if (finite_ub[j]) con <- c(con, ub[j] - xs[j])
    }
    list(f = fv, con = con)
  }
  stop$nevals <- 0L
  res <- .ncob_cobylb(n, m, x / scale, rhobeg, rhoend, stop, lb, ub, calcfc, con_tol)
  xr <- res$x * scale
  xr <- ifelse(xr < lb0, lb0, xr)
  xr <- ifelse(xr > ub0, ub0, xr)
  list(ret = res$ret, x = xr, minf = res$minf)
}
