# The Luksan algorithms of NLopt 2.11.0, ported to R from src/algs/luksan
# (the C translation by Steven G. Johnson of Ladislav Luksan's Fortran
# codes): the limited memory BFGS method based on the Strang recurrences
# (LD_LBFGS, plis.c), the shifted limited memory variable metric methods of
# rank one and two (LD_VAR1, LD_VAR2, plip.c), the truncated Newton methods
# with and without restarts and preconditioning (the LD_TNEWTON family,
# pnet.c) and their shared subroutines (mssubs.c and pssubs.c).
#
# The NLopt library that NLopt.jl loads (NLopt_jll) rounds every product
# and sum separately (no fused multiply add), so plain R arithmetic in the
# order of the C expressions reproduces it.
#
# The Luksan codes are Copyright (c) Ladislav Luksan, released by him under
# the GNU LGPL and included in NLopt under the MIT license by permission
# (see src/algs/luksan/COPYRIGHT of NLopt); NLopt is Copyright (c)
# 2007-2014 Massachusetts Institute of Technology, released under the MIT
# license (see R/optimization-nlopt-util.R).

# MAX2 and MIN2 of the C code
.nlu_max <- function(a, b) if (isTRUE(a > b)) a else b
.nlu_min <- function(a, b) if (isTRUE(a < b)) a else b

# mxudot (job > 0, the only case NLopt uses: the indices with ix < 0 are
# skipped) and mxvdot (ix NULL): the dot product, fused term by term
.nlu_dot <- function(x, y, ix = NULL) {
  s <- 0
  for (i in seq_along(x)) if (is.null(ix) || ix[i] >= 0L) s <- (x[i] * y[i] + s)
  s
}

# mxudir: z := y + a * x on the free indices (the others keep z)
.nlu_udir <- function(a, x, y, z, ix) {
  k <- ix >= 0L
  z[k] <- (a * x[k] + y[k])
  z
}

# mxuneg: y := -x on the free indices, 0 elsewhere
.nlu_uneg <- function(x, ix) ifelse(ix >= 0L, -x, 0)

# mxdrcb: the backward part of the Strang formula (A, B: n by m columns)
.nlu_drcb <- function(m, A, B, u, v, x, ix) {
  for (i in seq_len(m)) {
    v[i] <- u[i] * .nlu_dot(x, A[, i], ix)
    x <- .nlu_udir(-v[i], B[, i], x, x, ix)
  }
  list(v = v, x = x)
}

# mxdrcf: the forward part of the Strang formula
.nlu_drcf <- function(m, A, B, u, v, x, ix) {
  for (i in rev(seq_len(m))) {
    temp <- u[i] * .nlu_dot(x, B[, i], ix)
    x <- .nlu_udir(v[i] - temp, A[, i], x, x, ix)
  }
  x
}

# mxdrsu: shift the first m - 1 columns of A and B and elements of u one
# place to the right
.nlu_drsu <- function(m, A, B, u) {
  if (m >= 2L) for (i in (m - 1L):1L) {
    A[, i + 1L] <- A[, i]
    B[, i + 1L] <- B[, i]
    u[i + 1L] <- u[i]
  }
  list(A = A, B = B, u = u)
}

# a matrix with at least k columns (the columns NLopt allocates up front,
# all zero, are added as they are reached)
.nlu_cols <- function(A, k) {
  if (ncol(A) < k) A <- cbind(A, matrix(0, nrow(A), max(k, 2L * ncol(A)) - ncol(A)))
  A
}

# mxdrmm: y_j = column j of A (rowwise row j) times x, fused
.nlu_drmm <- function(m, A, x) vapply(seq_len(m), function(j) .nlu_dot(A[, j], x), 0)

# pcbs04: snap the variables within eps9 of an active bound onto it
.nlu_pcbs04 <- function(x, ix, xl, xu, eps9, kbf) {
  if (kbf > 0L) for (i in seq_along(x)) {
    ixi <- abs(ix[i])
    if ((ixi == 1L || ixi == 3L || ixi == 4L) && x[i] <= (eps9 * .nlu_max(abs(xl[i]), 1) + xl[i])) x[i] <- xl[i]
    if ((ixi == 2L || ixi == 3L || ixi == 4L) && x[i] >= (-eps9 * .nlu_max(abs(xu[i]), 1) + xu[i])) x[i] <- xu[i]
  }
  x
}

# pyadc0: add the box constraints the variables reach (e: the driver)
.nlu_pyadc0 <- function(e) {
  nf <- length(e$x)
  e$nact <- nf
  e$inew <- 0L
  for (i in seq_len(nf)) {
    ii <- e$ix[i]
    ixi <- abs(ii)
    if (ixi >= 5L) {
      e$ix[i] <- -ixi
    } else if ((ixi == 1L || ixi == 3L || ixi == 4L) && e$x[i] <= e$xl[i]) {
      e$x[i] <- e$xl[i]
      e$ix[i] <- if (ixi == 4L) -3L else -ixi
      e$nact <- e$nact - 1L
      if (ii > 0L) e$inew <- e$inew + 1L
    } else if ((ixi == 2L || ixi == 3L || ixi == 4L) && e$x[i] >= e$xu[i]) {
      e$x[i] <- e$xu[i]
      e$ix[i] <- if (ixi == 3L) -4L else -ixi
      e$nact <- e$nact - 1L
      if (ii > 0L) e$inew <- e$inew + 1L
    }
  }
  invisible()
}

# pytrcg: the projected gradient (umax, gmax, iold) of gradient g
.nlu_pytrcg <- function(e, g) {
  if (e$kbf > 0L) {
    gmax <- 0; umax <- 0; iold <- 0L
    for (i in seq_along(g)) {
      temp <- g[i]
      ixi <- e$ix[i]
      if (ixi >= 0L) {
        gmax <- .nlu_max(gmax, abs(temp))
      } else if (ixi <= -5L) {
      } else if ((ixi == -1L || ixi == -3L) && umax + temp >= 0) {
      } else if ((ixi == -2L || ixi == -4L) && umax - temp >= 0) {
      } else {
        iold <- i
        umax <- abs(temp)
      }
    }
    e$gmax <- gmax; e$umax <- umax; e$iold <- iold
  } else {
    e$umax <- 0
    m <- 0
    for (gi in g) m <- .nlu_max(m, abs(gi))
    e$gmax <- m
  }
  e$nact <- length(g)
  invisible()
}

# pyfut1: the termination tests
.nlu_pyfut1 <- function(e) {
  stop <- e$stop
  if (e$iterm < 0L) return(invisible())
  tests <- e$ites > 0L && e$iters != 0L
  if (tests) {
    if (e$nit <= 0L) e$fo <- e$f + .nlu_min(sqrt(abs(e$f)), abs(e$f) / 10)
    if (.nlopt_stop_forced(stop)) { e$iterm <- -999L; return(invisible()) }
    if (e$f <= stop$minf_max) { e$iterm <- 3L; return(invisible()) }
    if (e$kd > 0L && e$gmax <= e$tolg && e$umax <= e$tolg) { e$iterm <- 4L; return(invisible()) }
    if (e$nit <= 0L) { e$ntesx <- 0L; e$ntesf <- 0L }
    if (e$xstop) {
      e$iterm <- 1L
      e$ntesx <- e$ntesx + 1L
      if (e$ntesx >= e$mtesx) return(invisible())
    } else e$ntesx <- 0L
    if (.nlopt_stop_ftol(stop, e$f, e$fo)) {
      e$iterm <- 2L
      e$ntesf <- e$ntesf + 1L
      if (e$ntesf >= e$mtesf) return(invisible())
    } else e$ntesf <- 0L
  }
  if (e$nit >= e$mit) { e$iterm <- 11L; return(invisible()) }
  if (.nlopt_stop_evals(stop)) { e$iterm <- 12L; return(invisible()) }
  if (e$nfg >= e$mfg) { e$iterm <- 13L; return(invisible()) }
  e$iterm <- 0L
  nf <- length(e$x)
  if (nf > 0L && e$nit - e$kit >= e$ires1 * nf + e$ires2) e$irest <- max(e$irest, 1L)
  e$nit <- e$nit + 1L
  invisible()
}

# pyrmc0: delete the box constraints whose multipliers have the wrong sign
.nlu_pyrmc0 <- function(e, g) {
  if (e$nact == 0L || e$rmax > 0) {
    if (e$umax > e$eps8 * e$gmax) {
      e$iold <- 0L
      for (i in seq_along(g)) {
        ixi <- e$ix[i]
        if (ixi >= 0L) {
        } else if (ixi <= -5L) {
        } else if ((ixi == -1L || ixi == -3L) && -g[i] <= 0) {
        } else if ((ixi == -2L || ixi == -4L) && g[i] <= 0) {
        } else {
          e$iold <- e$iold + 1L
          e$ix[i] <- min(abs(e$ix[i]), 3L)
          if (e$rmax == 0) break
        }
      }
      if (e$iold > 1L) e$irest <- max(e$irest, 1L)
    }
  }
  invisible()
}

# pytrcd: the differences xo := x - xo and go := g - go of the accepted
# step (or the restore of the old point), and dmax
.nlu_pytrcd <- function(e) {
  if (e$iters > 0L) {
    e$xo[] <- e$x - e$xo
    e$go[] <- e$gf - e$go
    e$po <- e$r * e$po
    e$p <- e$r * e$p
  } else {
    e$f <- e$fo
    e$p <- e$po
    t <- e$xo; e$xo[] <- e$x - e$xo; e$x <- t
    t <- e$go; e$go[] <- e$gf - e$go; e$gf <- t
    e$ld <- e$kd
  }
  dmax <- 0
  for (i in seq_along(e$x)) {
    if (e$kbf > 0L && e$ix[i] < 0L) {
      e$xo[i] <- 0; e$go[i] <- 0
      next
    }
    dmax <- .nlu_max(dmax, abs(e$xo[i]) / .nlu_max(abs(e$x[i]), 1))
  }
  e$dmax <- dmax
  invisible()
}

# pytrcs: save the point (xo, go, fo, po) and limit rmax by the bounds
.nlu_pytrcs <- function(e) {
  e$fp <- e$fo
  e$ro <- 0
  e$fo <- e$f
  e$po <- e$p
  e$xo[] <- e$x
  e$go[] <- e$gf
  if (e$kbf > 0L) for (i in seq_along(e$x)) {
    ixi <- e$ix[i]
    if (ixi < 0L) {
      e$s[i] <- 0
    } else {
      if ((ixi == 1L || ixi >= 3L) && e$s[i] < -1 / e$eta9)
        e$rmax <- .nlu_min(e$rmax, (e$xl[i] - e$x[i]) / e$s[i])
      if ((ixi == 2L || ixi >= 3L) && e$s[i] > 1 / e$eta9)
        e$rmax <- .nlu_min(e$rmax, (e$xu[i] - e$x[i]) / e$s[i])
    }
  }
  invisible()
}

# pnint1: extrapolation or interpolation of the step with directional
# derivatives (the state ls holds rl, ru, fl, fu, pl, pu, mode, mtyp; sets
# ls$r and ls$merr)
.nlu_pnint1 <- function(ls) {
  ls$merr <- 0L
  if (ls$mode <= 0L) return(invisible())
  if (ls$pl >= 0) { ls$merr <- 2L; return(invisible()) }
  if (ls$ru <= ls$rl) { ls$merr <- 3L; return(invisible()) }
  rl <- ls$rl; ru <- ls$ru
  a <- 0; b <- 0; den <- 0
  for (ntyp in rev(seq_len(ls$mtyp))) {
    if (ntyp == 1L) {
      ls$r <- if (ls$mode == 1L) ru * 4 else (rl + ru) * 0.5
      return(invisible())
    } else if (ntyp == ls$mtyp) {
      a <- (ls$fu - ls$fl) / (ls$pl * (ru - rl))
      b <- ls$pu / ls$pl
    }
    ok <- TRUE
    if (ntyp == 2L) {
      den <- (1 - a) * 2
    } else if (ntyp == 3L) {
      den <- 1 - b
    } else if (ntyp == 4L) {
      cc <- (-a * 2 + b) + 1
      d <- (-a * 3 + b) + 2
      dis <- (d * d - cc * 3)
      if (dis < 0) ok <- FALSE else den <- d + sqrt(dis)
    } else if (ntyp == 5L) {
      dis <- (a * a + -b)
      if (dis < 0) ok <- FALSE else {
        den <- a + sqrt(dis)
        if (den <= 0) ok <- FALSE else {
          d1 <- 1 / den
          den <- (-b * (d1 * (d1 * d1)) + 1)
        }
      }
    }
    if (!ok) next
    if (ls$mode == 1L && den > 0 && den < 1) {
      r <- rl + (ru - rl) / den
      r <- .nlu_max(r, ru * 1.1)
      ls$r <- .nlu_min(r, ru * 1e3)
      return(invisible())
    } else if (ls$mode == 2L && den > 1) {
      r <- rl + (ru - rl) / den
      r <- if (rl == 0) .nlu_max(r, ((ru - rl) * 0.01 + rl)) else .nlu_max(r, ((ru - rl) * 0.1 + rl))
      ls$r <- .nlu_min(r, ((ru - rl) * 0.9 + rl))
      return(invisible())
    }
  }
  invisible()
}

# ps1l01: the standard line search with directional derivatives, as a
# reverse communication routine (e: the driver, whose fields r, rp, f, fo,
# fp, p, po, pp, rmin, rmax, kd, ld, nit, kit, nred, maxst, iters, isys it
# reads and sets; ls: the saved state of the search)
.nlu_ps1l01 <- function(e, ls) {
  tols <- e$tols; tolp <- e$tolp; minf <- e$minf_est
  finish <- function() { e$isys <- 0L; invisible() }
  if (e$isys != 1L) {
    ls$mes1 <- 2L; ls$mes2 <- 2L; ls$mes3 <- 2L
    e$iters <- 0L
    if (e$po >= 0) { e$r <- 0; e$iters <- -2L; return(finish()) }
    if (e$rmax <= 0) { e$iters <- 0L; return(finish()) }
    rtemp <- if (e$inits > 0L) minf - e$f else if (e$iest == 0L) e$f - e$fp else
      .nlu_max(e$f - e$fp, minf - e$f)
    init1 <- abs(e$inits)
    e$rp <- 0
    e$fp <- e$fo
    e$pp <- e$po
    if (init1 == 0L) {
    } else if (init1 == 1L || (e$inits >= 1L && e$iest == 0L)) {
      e$r <- 1
    } else if (init1 == 2L) {
      e$r <- .nlu_min(1, rtemp * 4 / e$po)
    } else if (init1 == 3L) {
      e$r <- .nlu_min(1, rtemp * 2 / e$po)
    } else if (init1 == 4L) {
      e$r <- rtemp * 2 / e$po
    }
    e$r <- .nlu_max(e$r, e$rmin)
    e$r <- .nlu_min(e$r, e$rmax)
    ls$mode <- 0L
    ls$ru <- 0
    ls$fu <- e$fo
    ls$pu <- e$po
  } else {
    # L3
    if (ls$mode == 0L) {
      e$par1 <- e$p / e$po
      e$par2 <- e$f - e$fo
    }
    if (e$iters != 0L) return(finish())
    if (e$f <= minf) { e$iters <- 7L; return(finish()) }
    l1 <- e$r <= e$rmin && e$nit != e$kit
    l2 <- e$r >= e$rmax
    l3 <- e$f - e$fo <= tols * e$r * e$po
    l5 <- e$p >= tolp * e$po || (ls$mes2 == 2L && ls$mode == 2L)
    l7 <- ls$mes2 <= 2L || ls$mode != 0L
    m3 <- l3
    if (ls$mes3 >= 1L) {
      m1 <- abs(e$p) <= abs(e$po) * 0.01 && e$fo - e$f >= abs(e$fo) * 0x1.5fd7fe1796495p-37
      l3 <- l3 || m1
    }
    if (ls$mes3 >= 2L) {
      m2 <- abs(e$p) <= abs(e$po) * 0.5 && abs(e$fo - e$f) <= abs(e$fo) * 2.0000000000000001e-13
      l3 <- l3 || m2
    }
    e$maxst <- if (l2) 1L else 0L
    if (l1 && !l3) {
      e$iters <- 0L; return(finish())
    } else if (l2 && l3 && !l5) {
      e$iters <- 7L; return(finish())
    } else if (m3 && ls$mes1 == 3L) {
      e$iters <- 5L; return(finish())
    } else if (l3 && l5 && l7) {
      e$iters <- 4L; return(finish())
    } else if (e$kters < 0L || (e$kters == 6L && l7)) {
      e$iters <- 6L; return(finish())
    } else if (abs(e$nred) >= e$mred) {
      e$iters <- -1L; return(finish())
    } else {
      e$rp <- e$r
      e$fp <- e$f
      e$pp <- e$p
      ls$mode <- max(ls$mode, 1L)
      ls$mtyp <- abs(e$mes)
      if (e$f >= e$maxf) ls$mtyp <- 1L
    }
    if (ls$mode == 1L) {
      # interval change after extrapolation
      ls$rl <- ls$ru; ls$fl <- ls$fu; ls$pl <- ls$pu
      ls$ru <- e$r; ls$fu <- e$f; ls$pu <- e$p
      if (!l3) {
        e$nred <- 0L
        ls$mode <- 2L
      } else if (ls$mes1 == 1L) {
        ls$mtyp <- 1L
      }
    } else {
      # interval change after interpolation
      if (!l3) {
        ls$ru <- e$r; ls$fu <- e$f; ls$pu <- e$p
      } else {
        ls$rl <- e$r; ls$fl <- e$f; ls$pl <- e$p
      }
    }
  }
  # L2: the new stepsize
  ls$r <- e$r
  .nlu_pnint1(ls)
  e$r <- ls$r
  if (ls$merr > 0L) {
    e$iters <- -ls$merr
    return(finish())
  } else if (ls$mode == 1L) {
    e$nred <- e$nred - 1L
    e$r <- .nlu_min(e$r, e$rmax)
  } else if (ls$mode == 2L) {
    e$nred <- e$nred + 1L
  }
  e$kd <- 1L
  e$ld <- -1L
  e$isys <- 1L
  invisible()
}

.nlu_ls_state <- function() {
  ls <- new.env(parent = emptyenv())
  ls$fl <- 0; ls$fu <- 0; ls$pl <- 0; ls$rl <- 0; ls$pu <- 0; ls$ru <- 0
  ls$mes1 <- 0L; ls$mes2 <- 0L; ls$mes3 <- 0L; ls$mode <- 0L; ls$mtyp <- 0L
  ls$r <- 0; ls$merr <- 0L
  ls
}

# the common initiation of plis_, plip_ and pnet_ in the environment e of
# the driver: the bound types, the counters and the tolerances
.nlu_init <- function(e, n, lb, ub, x, stop, tolg, tolp) {
  e$stop <- stop
  e$ix <- vapply(seq_len(n), function(i) {
    lbu <- lb[i] <= -0.99 * Inf
    ubu <- ub[i] >= 0.99 * Inf
    if (lbu) (if (ubu) 0L else 2L) else if (ubu) 1L else if (lb[i] == ub[i]) 5L else 3L
  }, 0L)
  e$xl <- lb; e$xu <- ub
  e$x <- x
  e$kbf <- 2L # nb = 1
  e$nres <- 0L; e$nit <- 0L; e$nfg <- 0L
  e$isys <- 0L; e$ites <- 1L; e$mtesx <- 2L; e$mtesf <- 2L; e$inits <- 2L; e$iterm <- 0L
  e$iterd <- 0L; e$iters <- 2L; e$kters <- 3L; e$irest <- 0L; e$ires1 <- 999L; e$ires2 <- 0L
  e$mred <- 10L; e$mes <- 4L; e$eta9 <- 0x1.8c8dac6a0342ap+398; e$eps8 <- 1; e$eps9 <- 1e-8; e$alf1 <- 1e-10
  e$alf2 <- 1e10; e$rmax <- e$eta9; e$dmax <- e$eta9; e$maxf <- 1e20
  e$minf_est <- -Inf
  e$iest <- 0L
  e$xmax <- 1e16
  # the tolerances are passed by pointer into the stopping criteria
  if (stop$xtol_rel <= 0) stop$xtol_rel <- 1e-16
  if (stop$ftol_rel <= 0) stop$ftol_rel <- 1e-14
  e$tolg <- if (tolg <= 0) 1e-8 else tolg
  e$told <- 1e-4; e$tols <- 1e-4; e$tolp <- tolp
  e$mit <- .Machine$integer.max
  e$mfg <- if (stop$maxeval > 0L) stop$maxeval else .Machine$integer.max
  e$kd <- 1L; e$ld <- -1L
  e$kit <- -(e$ires1 * n + e$ires2)
  e$fo <- e$minf_est
  e$xstop <- FALSE
  e$ntesx <- 0L; e$ntesf <- 0L
  e$r <- 0; e$rp <- 0; e$fp <- 0; e$pp <- 0; e$p <- 0; e$po <- 0; e$ro <- 0; e$rmin <- 0
  e$nred <- 0L; e$maxst <- 0L; e$par1 <- 0; e$par2 <- 0
  e$umax <- 0; e$gmax <- 0; e$iold <- 0L; e$nact <- n; e$inew <- 0L
  if (e$kbf > 0L) {
    for (i in seq_len(n)) {
      if ((e$ix[i] == 3L || e$ix[i] == 4L) && e$xu[i] <= e$xl[i]) {
        e$xu[i] <- e$xl[i]; e$ix[i] <- 5L
      } else if (e$ix[i] == 5L || e$ix[i] == 6L) {
        e$xl[i] <- e$x[i]; e$xu[i] <- e$x[i]; e$ix[i] <- 5L
      }
    }
    e$x <- .nlu_pcbs04(e$x, e$ix, e$xl, e$xu, e$eps9, e$kbf)
    .nlu_pyadc0(e)
  }
  invisible()
}

# one objective and gradient evaluation at e$x into e$f and e$gf
.nlu_eval <- function(e, fg) {
  r <- fg(e$x)
  e$stop$nevals <- e$stop$nevals + 1L
  e$nfg <- e$nfg + 1L
  e$f <- r$val
  e$gf <- as.numeric(r$grad)
  invisible()
}

# the default number of stored steps (mf <= 0) of luksan_plis, plip, pnet
.nlu_mf <- function(mf, n, stop) {
  if (mf <= 0L) {
    mf <- max(1310720L %/% n, 10L)
    if (stop$maxeval > 0L && stop$maxeval <= mf) mf <- max(stop$maxeval, 1L)
  }
  as.integer(mf)
}

.nlu_result <- function(e) {
  ret <- switch(as.character(e$iterm), "1" = "XTOL_REACHED", "2" = "FTOL_REACHED", "3" = "STOPVAL_REACHED",
                "4" = , "6" = "SUCCESS", "12" = , "13" = "MAXEVAL_REACHED", "100" = "MAXTIME_REACHED",
                "-999" = "FORCED_STOP", "FAILURE")
  list(ret = ret, x = e$x, minf = e$f)
}

# the line search loop shared by the three drivers (L11170 of plis_):
# returns after the search ended (isys == 0)
.nlu_linesearch <- function(e, ls, fg) {
  repeat {
    .nlu_ps1l01(e, ls)
    if (e$isys == 0L) return(invisible())
    e$x <- .nlu_udir(e$r, e$s, e$xo, e$x, e$ix)
    e$x <- .nlu_pcbs04(e$x, e$ix, e$xl, e$xu, e$eps9, e$kbf)
    .nlu_eval(e, fg)
    e$ld <- e$kd
    e$p <- .nlu_dot(e$gf, e$s, e$ix)
  }
}

# the test on the descent direction and the preparation of the line search
.nlu_descent <- function(e) {
  if (e$iterd < 0L) {
    e$iterm <- e$iterd
  } else {
    if (e$snorm <= 0) {
      e$irest <- max(e$irest, 1L)
    } else if ((e$told * e$gnorm * e$snorm + e$p) <= 0) {
      e$irest <- 0L
    } else {
      e$irest <- max(e$irest, 1L)
    }
    if (e$irest == 0L) {
      e$nred <- 0L
      e$rmin <- e$alf1 * e$gnorm / e$snorm
      e$rmax <- .nlu_min(e$alf2 * e$gnorm / e$snorm, e$xmax / e$snorm)
    }
  }
  invisible()
}

# luksan_plis: LD_LBFGS
.nlopt_luksan_plis <- function(fg, lb, ub, x, stop, mf, tolg) {
  n <- length(x)
  mf <- .nlu_mf(mf, n, stop)
  e <- environment()
  ls <- .nlu_ls_state()
  .nlu_init(e, n, lb, ub, x, stop, tolg, tolp = 0.8)
  XO <- matrix(0, n, 1L); GO <- matrix(0, n, 1L)
  xo <- numeric(n); go <- numeric(n)
  uo <- numeric(1L); vo <- numeric(1L)
  s <- numeric(n)
  snorm <- 0; gnorm <- 0
  if (iterm != 0L) return(.nlu_result(e))
  .nlu_eval(e, fg)
  if (.nlopt_stop_time(stop)) { iterm <- 100L; return(.nlu_result(e)) }
  repeat {
    # L11120
    .nlu_pytrcg(e, gf)
    .nlu_pyfut1(e)
    if (iterm != 0L) break
    if (.nlopt_stop_time(stop)) { iterm <- 100L; break }
    if (kbf > 0L && rmax > 0) .nlu_pyrmc0(e, gf)
    done <- FALSE
    repeat {
      # L11130: direction determination
      gnorm <- sqrt(.nlu_dot(gf, gf, ix))
      XO[, 1L] <- xo; GO[, 1L] <- go
      if (irest == 0L) {
        k <- min(nit - kit, mf)
        if (k <= 0L) {
          irest <- max(irest, 1L)
        } else {
          b <- .nlu_dot(xo, go, ix)
          if (b <= 0) {
            irest <- max(irest, 1L)
          } else {
            if (length(uo) < k + 1L) { uo <- c(uo, numeric(k + 1L - length(uo))); vo <- c(vo, numeric(k + 1L - length(vo))) }
            XO <- .nlu_cols(XO, k); GO <- .nlu_cols(GO, k)
            uo[1L] <- 1 / b
            s <- .nlu_uneg(gf, ix)
            r1 <- .nlu_drcb(k, XO, GO, uo, vo, s, ix)
            vo <- r1$v; s <- r1$x
            a <- .nlu_dot(go, go, ix)
            if (a > 0) s <- (b / a) * s
            s <- .nlu_drcf(k, XO, GO, uo, vo, s, ix)
            snorm <- sqrt(.nlu_dot(s, s, ix))
            k <- min(k + 1L, mf)
            XO <- .nlu_cols(XO, k); GO <- .nlu_cols(GO, k)
            r2 <- .nlu_drsu(k, XO, GO, uo)
            XO <- r2$A; GO <- r2$B; uo <- r2$u
            xo <- XO[, 1L]; go <- GO[, 1L]
          }
        }
      }
      # L12620
      iterd <- 0L
      if (irest != 0L) {
        s <- .nlu_uneg(gf, ix)
        snorm <- gnorm
        if (kit < nit) {
          nres <- nres + 1L
          kit <- nit
        } else {
          iterm <- -10L
          if (iters < 0L) iterm <- iters - 5L
        }
      }
      if (kd > 0L) p <- .nlu_dot(gf, s, ix)
      .nlu_descent(e)
      if (iterm != 0L) { done <- TRUE; break }
      if (.nlopt_stop_time(stop)) { iterm <- 100L; done <- TRUE; break }
      if (irest != 0L) next
      .nlu_pytrcs(e)
      if (rmax == 0) break
      .nlu_linesearch(e, ls, fg)
      if (iters <= 0L) {
        r <- 0
        f <- fo
        p <- po
        x <- xo
        gf <- go
        irest <- max(irest, 1L)
        ld <- kd
        next
      }
      .nlu_pytrcd(e)
      xstop <- .nlopt_stop_dx(stop, x, xo)
      break
    }
    if (done) break
    # L11175
    if (kbf > 0L) {
      ix <- abs(ix)
      .nlu_pyadc0(e)
    }
  }
  .nlu_result(e)
}

# pulsp3: the limited storage variable metric update (shifted BFGS) while
# the matrix XM has fewer than mf columns (e: the driver, using nn, XM, gr,
# xo, go, par, iterh)
.nlu_pulsp3 <- function(e, met3) {
  n <- length(e$xo)
  if (e$nn >= e$mf) return(invisible())
  b <- .nlu_dot(e$xo, e$go)
  if (b <= 0) { e$iterh <- 2L; return(invisible()) }
  m <- e$nn
  gr <- .nlu_drmm(m, e$XM, e$go)
  ah <- .nlu_dot(e$go, e$go)
  aa <- .nlu_dot(gr, gr)
  sig <- e$par
  a <- (ah * sig + aa)
  par <- 1
  pom <- b / ah
  if (a > 0) {
    den <- .nlu_dot(e$xo, e$xo)
    sig <- if (met3 <= 4L) {
      sqrt(.nlu_max(0, 1 - aa / a)) / (sqrt(.nlu_max(0, 1 - b * b / (den * ah))) + 1) * pom
    } else {
      sqrt(.nlu_max(0, sig * ah / a)) / (sqrt(.nlu_max(0, 1 - b * b / (den * ah))) + 1) * pom
    }
    sig <- .nlu_max(sig, pom * 0.2)
    sig <- .nlu_min(sig, pom * 0.8)
  } else {
    sig <- pom * 0.25
  }
  bb <- (-ah * sig + b)
  e$xo <- (-sig * e$go + e$xo)
  XM <- e$XM
  if (m > 0L) for (j in seq_len(m)) XM[, j] <- ((-1 / bb) * gr[j] * e$xo + XM[, j])
  XM <- .nlu_cols(XM, m + 1L)
  XM[, m + 1L] <- sqrt(par / bb) * e$xo
  e$XM <- XM
  e$nn <- m + 1L
  e$par <- sig
  if (length(gr)) e$gr[seq_along(gr)] <- gr
  e$iterh <- 0L
  invisible()
}

# pulvp3: the shifted limited memory variable metric update of rank one
# (met5 = 1) or two once XM is full
.nlu_pulvp3 <- function(e, met2, met3, met5) {
  zet <- e$par
  b <- .nlu_dot(e$xo, e$go)
  if (b <= 0) { e$iterh <- 2L; return(invisible()) }
  m <- e$nn
  XM <- e$XM
  gr <- .nlu_drmm(m, XM, e$go)
  s <- e$r * e$s
  xr <- .nlu_drmm(m, XM, s)
  sig <- e$par
  so <- (-sig * s + e$xo)
  ah <- .nlu_dot(e$go, e$go)
  aa <- .nlu_dot(gr, gr)
  bb <- .nlu_dot(gr, xr)
  cc <- .nlu_dot(xr, xr)
  a <- (ah * sig + aa)
  pom <- b / ah
  if (a > 0) {
    den <- .nlu_dot(e$xo, e$xo)
    sig <- if (met3 <= 4L) {
      sqrt(.nlu_max(0, 1 - aa / a)) / (sqrt(.nlu_max(0, 1 - b * b / (den * ah))) + 1) * pom
    } else {
      sqrt(.nlu_max(0, sig * ah / a)) / (sqrt(.nlu_max(0, 1 - b * b / (den * ah))) + 1) * pom
    }
    sig <- .nlu_max(sig, pom * 0.2)
    sig <- .nlu_min(sig, pom * 0.8)
  } else {
    sig <- pom * 0.25
  }
  b <- (-ah * sig + b)
  xo <- (-sig * e$go + e$xo)
  par <- if (met2 <= 1L) 1 else if (met2 == 2L) sig * ah / b else if (met2 == 3L) sqrt(1 - aa / a) else
    if (met2 == 4L) sqrt(sqrt(1 - aa / a) * (sig * ah / b)) else zet / (zet + sig)
  d1 <- sqrt(par * b / cc)
  pom <- if (bb < 0 || (bb == 0 && 1 / bb < 0)) -abs(d1) else abs(d1)
  if (met5 == 1L) {
    xr <- (pom * xr + gr)
    s <- (par * xo + pom * so)
    d1 <- -1 / (par * b + pom * bb)
    for (j in seq_len(m)) XM[, j] <- (d1 * xr[j] * s + XM[, j])
  } else {
    s <- ((par / pom - bb / b) * xo + so)
    ta <- -1 / b; tb <- -1 / cc
    for (j in seq_len(m)) XM[, j] <- (tb * xr[j] * s + (ta * gr[j] * xo + XM[, j]))
  }
  e$XM <- XM
  e$xo <- xo; e$s <- s; e$so <- so
  e$par <- sig
  e$iterh <- 0L
  invisible()
}

# luksan_plip: LD_VAR1 (method 1) and LD_VAR2 (method 2)
.nlopt_luksan_plip <- function(fg, lb, ub, x, stop, mf, method, tolg) {
  n <- length(x)
  mf <- .nlu_mf(mf, n, stop)
  e <- environment()
  ls <- .nlu_ls_state()
  .nlu_init(e, n, lb, ub, x, stop, tolg, tolp = 0.9)
  iterh <- 0L
  met3 <- 4L; mec <- 4L
  met <- if (method <= 0L) 2L else as.integer(method)
  xo <- numeric(n); go <- numeric(n); so <- numeric(n); s <- numeric(n)
  XM <- matrix(0, n, 1L); gr <- numeric(1L)
  nn <- 0L; par <- 1
  snorm <- 0; gnorm <- 0
  if (iterm != 0L) return(.nlu_result(e))
  .nlu_eval(e, fg)
  if (.nlopt_stop_time(stop)) { iterm <- 100L; return(.nlu_result(e)) }
  repeat {
    # L11120
    .nlu_pytrcg(e, gf)
    .nlu_pyfut1(e)
    if (iterm != 0L) break
    if (.nlopt_stop_time(stop)) { iterm <- 100L; break }
    if (kbf > 0L && rmax > 0) .nlu_pyrmc0(e, gf)
    done <- FALSE
    repeat {
      # L11130
      if (irest > 0L) {
        nn <- 0L
        par <- 1
        ld <- min(ld, 1L)
        if (kit < nit) {
          nres <- nres + 1L
          kit <- nit
        } else {
          iterm <- -10L
          if (iters < 0L) iterm <- iters - 5L
        }
      }
      if (iterm != 0L) { done <- TRUE; break }
      if (.nlopt_stop_time(stop)) { iterm <- 100L; done <- TRUE; break }
      gnorm <- sqrt(.nlu_dot(gf, gf, ix))
      # Newton like step
      s <- .nlu_uneg(gf, ix)
      grn <- .nlu_drmm(nn, XM, s)
      z <- par * s
      for (j in seq_len(nn)) z <- (grn[j] * XM[, j] + z)
      s <- z
      if (nn > 0L) gr[seq_len(nn)] <- grn
      s[ix < 0L] <- 0
      iterd <- 1L
      snorm <- sqrt(.nlu_dot(s, s, ix))
      if (kd > 0L) p <- .nlu_dot(gf, s, ix)
      .nlu_descent(e)
      if (iterm != 0L) { done <- TRUE; break }
      if (.nlopt_stop_time(stop)) { iterm <- 100L; done <- TRUE; break }
      if (irest != 0L) next
      .nlu_pytrcs(e)
      if (rmax == 0) break
      .nlu_linesearch(e, ls, fg)
      if (iters <= 0L) {
        r <- 0
        f <- fo
        p <- po
        x <- xo
        gf <- go
        irest <- max(irest, 1L)
        ld <- kd
        next
      }
      s <- .nlu_uneg(go, ix)
      .nlu_pytrcd(e)
      xstop <- .nlopt_stop_dx(stop, x, xo)
      so <- ifelse(ix >= 0L, gf, 0)
      if (nn < mf) .nlu_pulsp3(e, met3) else .nlu_pulvp3(e, mec, met3, met)
      break
    }
    if (done) break
    # L11175
    if (iterh != 0L) irest <- max(irest, 1L)
    if (kbf > 0L) .nlu_pyadc0(e)
  }
  .nlu_result(e)
}

# luksan_pnet: the LD_TNEWTON family (mos1 = 2: restarts by the steepest
# descent direction, mos2 = 2: the limited memory BFGS preconditioner)
.nlopt_luksan_pnet <- function(fg, lb, ub, x, stop, mf, mos1, mos2, tolg) {
  n <- length(x)
  mf <- .nlu_mf(mf, n, stop)
  e <- environment()
  ls <- .nlu_ls_state()
  .nlu_init(e, n, lb, ub, x, stop, tolg, tolp = 0.9)
  mfg <- .Machine$integer.max
  eps <- 0.8; eta0 <- 1e-15
  if (mos1 <= 0L) mos1 <- 1L
  if (mos2 <= 0L) mos2 <- 1L
  xo <- numeric(n); go <- numeric(n); s <- numeric(n); gn <- numeric(n); xs <- numeric(n); gs <- numeric(n)
  XM <- matrix(0, n, 1L); GM <- matrix(0, n, 1L); u1 <- numeric(1L); u2 <- numeric(1L)
  mx <- 0L; a <- 0; b <- 0
  snorm <- 0; gnorm <- 0
  .nlu_eval(e, fg)
  if (.nlopt_stop_time(stop)) { iterm <- 100L; return(.nlu_result(e)) }
  ld <- kd
  repeat {
    # L11020
    .nlu_pytrcg(e, gf)
    gn <- gf
    .nlu_pyfut1(e)
    if (iterm != 0L) break
    if (.nlopt_stop_time(stop)) { iterm <- 100L; break }
    if (kbf > 0L) {
      .nlu_pyrmc0(e, gn)
      if (umax > eps8 * gmax) irest <- max(irest, 1L)
    }
    xo <- x
    done <- FALSE
    repeat {
      # L11040: direction determination
      skip_cg <- FALSE
      if (irest != 0L) {
        if (kit < nit) {
          mx <- 0L
          nres <- nres + 1L
          kit <- nit
        } else {
          iterm <- -10L
          if (iters < 0L) iterm <- iters - 5L
          done <- TRUE
          break
        }
        if (mos1 > 1L) {
          s <- -gn
          gnorm <- sqrt(.nlu_dot(gn, gn, ix))
          snorm <- gnorm
          skip_cg <- TRUE
        }
      }
      if (!skip_cg) {
        rho1 <- .nlu_dot(gn, gn, ix)
        gnorm <- sqrt(rho1)
        par <- .nlu_min(eps, sqrt(gnorm))
        if (par > 0.01) par <- .nlu_min(par, 1 / nit)
        par <- par * par
        # CG initiation
        rho <- rho1
        snorm <- 0
        s <- numeric(n)
        gs <- -gn
        xs <- gs
        if (mos2 > 1L) {
          b <- if (mx == 0L) 0 else .nlu_dot(XM[, 1L], GM[, 1L], ix)
          if (b > 0) {
            if (length(u1) < mx) { u1 <- c(u1, numeric(mx - length(u1))); u2 <- c(u2, numeric(mx - length(u2))) }
            u1[1L] <- 1 / b
            r1 <- .nlu_drcb(mx, XM, GM, u1, u2, xs, ix)
            u2 <- r1$v; xs <- r1$x
            a <- .nlu_dot(GM[, 1L], GM[, 1L], ix)
            if (a > 0) xs <- (b / a) * xs
            xs <- .nlu_drcf(mx, XM, GM, u1, u2, xs, ix)
          }
        }
        rho <- .nlu_dot(gs, xs, ix)
        mmx <- n + 3L
        nred <- 0L
        repeat {
          # L12520
          nred <- nred + 1L
          if (nred > mmx) break
          fo <- f
          pp <- sqrt(eta0 / .nlu_dot(xs, xs, ix))
          ld <- 0L
          x <- .nlu_udir(pp, xs, xo, x, ix)
          .nlu_eval(e, fg)
          ld <- kd
          go <- gf - gn
          f <- fo
          go <- (1 / pp) * go
          alf <- .nlu_dot(xs, go, ix)
          if (alf <= 1 / eta9) {
            # CG fails (the matrix is not positive definite)
            if (nred == 1L) {
              s <- -gn
              snorm <- gnorm
            }
            iterd <- 0L
            break
          } else {
            iterd <- 2L
          }
          # CG step
          alf <- rho / alf
          s <- .nlu_udir(alf, xs, s, s, ix)
          gs <- .nlu_udir(-alf, go, gs, gs, ix)
          rho2 <- .nlu_dot(gs, gs, ix)
          snorm <- sqrt(.nlu_dot(s, s, ix))
          if (rho2 <= par * rho1) break
          if (nred >= mmx) break
          if (mos2 > 1L && b > 0) {
            go <- gs
            r1 <- .nlu_drcb(mx, XM, GM, u1, u2, go, ix)
            u2 <- r1$v; go <- r1$x
            if (a > 0) go <- (b / a) * go
            go <- .nlu_drcf(mx, XM, GM, u1, u2, go, ix)
            rho2 <- .nlu_dot(gs, go, ix)
            alf <- rho2 / rho
            xs <- .nlu_udir(alf, xs, go, xs, ix)
          } else {
            alf <- rho2 / rho
            xs <- .nlu_udir(alf, xs, gs, xs, ix)
          }
          rho <- rho2
        }
      }
      # L12560: end of the direction determination
      x <- xo
      gf <- gn
      if (kd > 0L) p <- .nlu_dot(gn, s, ix)
      .nlu_descent(e)
      ld <- kd
      if (iterm != 0L) { done <- TRUE; break }
      if (.nlopt_stop_time(stop)) { iterm <- 100L; done <- TRUE; break }
      if (irest != 0L) next
      .nlu_pytrcs(e)
      if (rmax == 0) break
      .nlu_linesearch(e, ls, fg)
      if (iters <= 0L) {
        r <- 0
        f <- fo
        p <- po
        x <- xo
        gf <- go
        irest <- max(irest, 1L)
        ld <- kd
        next
      }
      .nlu_pytrcd(e)
      xstop <- .nlopt_stop_dx(stop, x, xo)
      if (mos2 > 1L) {
        mx <- min(mx + 1L, mf)
        XM <- .nlu_cols(XM, mx); GM <- .nlu_cols(GM, mx)
        if (length(u1) < mx) { u1 <- c(u1, numeric(mx - length(u1))); u2 <- c(u2, numeric(mx - length(u2))) }
        r2 <- .nlu_drsu(mx, XM, GM, u1)
        XM <- r2$A; GM <- r2$B; u1 <- r2$u
        XM[, 1L] <- xo; GM[, 1L] <- go
      }
      break
    }
    if (done) break
    # L11075
    if (kbf > 0L) {
      .nlu_pyadc0(e)
      if (inew > 0L) irest <- max(irest, 1L)
    }
  }
  .nlu_result(e)
}

# the Luksan drivers evaluate in the state environment of .nlu_*, which
# holds these variables
utils::globalVariables(c("eps8", "eta9", "fo", "gmax", "inew", "iters", "ix", "kbf", "kd", "nit",
                         "po", "rmax", "umax"))
