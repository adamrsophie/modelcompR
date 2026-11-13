# Gablonsky's Fortran DIRECT of NLopt 2.11.0 (GN_ORIG_DIRECT and
# GN_ORIG_DIRECT_L, src/algs/direct: DIRect.c, DIRsubrout.c, DIRserial.c and
# direct_wrap.c, the f2c translation of DIRECT 2.0.4 by J. M. Gablonsky),
# ported to R with the f_direct wrapper of optimize.c, which marks a point
# infeasible when f is not finite or an inequality constraint is positive.
# The arrays keep the layout of the Fortran code: the centers c(n, maxfunc),
# the side length levels length(n, maxfunc), f(2, maxfunc) (value and
# feasibility flag), the linked lists point(maxfunc) and anchor(-1:maxdeep).
#
# Copyright (c) 2007-2014 Massachusetts Institute of Technology, released
# under the MIT license (see R/optimization-nlopt-util.R). DIRECT 2.0.4 is
# Copyright (c) 1999-2001 Joerg M. Gablonsky, under the MIT license of
# src/algs/direct/COPYING.

# direct_optimize: returns list(ierror, x, minf)
.nlopt_direct_optimize <- function(fcn, n, lower, upper, x, max_feval, max_iter, stop, magic_eps, magic_eps_abs,
                                   volume_reltol, sigma_reltol, fglobal, fglobal_reltol, algorithm,
                                   upstream_replaceinf = FALSE) {
  algmethod <- as.integer(algorithm == "GABLONSKY")
  volume_reltol <- volume_reltol * 100
  sigma_reltol <- sigma_reltol * 100
  fglobal_reltol <- fglobal_reltol * 100
  if (volume_reltol <= 0) volume_reltol <- -1
  if (sigma_reltol <= 0) sigma_reltol <- -1
  if (fglobal == -Inf) fglobal_reltol <- 0.0
  if (n < 1L) return(list(ierror = -101L, x = x, minf = Inf))
  volper <- volume_reltol; sigmaper <- sigma_reltol; fglper <- fglobal_reltol
  eps <- magic_eps; epsabs <- magic_eps_abs
  maxf <- as.integer(max_feval); maxt <- as.integer(max_iter)
  MAXFUNC <- if (maxf <= 0L) 101000L else maxf + 1000L + maxf %/% 2L
  MAXDEEP <- if (maxt <= 0L) MAXFUNC %/% 5L else maxt + 1000L
  MAXDIV <- 5000L
  if (maxf <= 0L) maxf <- MAXFUNC - 1000L
  if (maxt <= 0L) maxt <- MAXDEEP
  l <- as.numeric(lower); u <- as.numeric(upper)
  oldl <- l; oldu <- u
  cheat <- 0L
  kmax <- 1e10
  mdeep <- MAXDEEP
  jones <- algmethod

  Cm <- matrix(0, n, MAXFUNC)
  Lm <- matrix(0L, n, MAXFUNC)
  Fm <- matrix(0, 2L, MAXFUNC)
  point <- integer(MAXFUNC)
  AN <- integer(MAXDEEP + 2L) # anchor(j) is AN[j + 2], j = -1..maxdeep
  S <- matrix(0L, MAXDIV, 2L)
  thirds <- numeric(MAXDEEP + 1L) # thirds(i) is thirds[i + 1]
  levels <- numeric(MAXDEEP + 1L)
  w <- numeric(n)
  list2 <- matrix(0L, n, 2L)
  arrayi <- integer(n)
  xs <- x # the scratch point of direct_direct_
  minf <- Inf; minpos <- 1L; fmax <- 0
  ifeasiblef <- 0L; iinfesiblef <- 0L
  ifree <- 1L; ifreeold <- 0L; maxi <- 0L
  epsfix <- 0; iepschange <- 0L
  ierror <- 0L

  # ---- the subroutines (DIRsubrout.c, DIRserial.c)
  getlevel <- function(pos, jon) {
    lp <- Lm[, pos]
    if (jon == 0L) {
      help <- lp[1L]; k <- help; p <- 1L
      for (i in seq_len(n)[-1L]) {
        if (lp[i] < k) k <- lp[i]
        if (lp[i] == help) p <- p + 1L
      }
      if (k == help) k * n + n - p else k * n + p
    } else min(lp)
  }
  getmaxdeep <- function(pos) min(Lm[, pos])
  get_i <- function(pos) {
    lp <- Lm[, pos]
    a <- which(lp == min(lp))
    arrayi[seq_along(a)] <<- a
    maxi <<- length(a)
  }
  infcn <- function(xx) {
    # dirinfcn: (x + c2) * c1 to the original domain, then back
    xx <- (xx + u) * l
    r <- fcn(xx)
    xs <<- xx / l - u
    r
  }
  samplepoints <- function(sample, delta) {
    pos <- ifree
    start <- ifree
    for (k in seq_len(maxi + maxi)) {
      Lm[, ifree] <<- Lm[, sample]
      Cm[, ifree] <<- Cm[, sample]
      pos <- ifree
      ifree <<- point[ifree]
      if (ifree == 0L) return(list(oops = 1L, start = start))
    }
    point[pos] <<- 0L
    pos <- start
    for (j in seq_len(maxi)) {
      ai <- arrayi[j]
      Cm[ai, pos] <<- Cm[ai, sample] + delta
      pos <- point[pos]
      Cm[ai, pos] <<- Cm[ai, sample] - delta
      pos <- point[pos]
    }
    list(oops = 0L, start = start)
  }
  samplef <- function(new) {
    pos <- new
    for (j in seq_len(maxi + maxi)) {
      xs <<- Cm[, pos]
      if (.nlopt_stop_forced(stop)) {
        Fm[1L, pos] <<- fmax
        kret <- -1L
      } else {
        r <- infcn(xs)
        Fm[1L, pos] <<- r$f
        kret <- r$flag
        if (.nlopt_stop_forced(stop)) kret <- -1L
      }
      if (kret > iinfesiblef) iinfesiblef <<- kret
      if (kret == 0L) {
        Fm[2L, pos] <<- 0
        ifeasiblef <<- 0L
        if (Fm[1L, pos] > fmax) fmax <<- Fm[1L, pos]
      }
      if (kret >= 1L) {
        Fm[2L, pos] <<- 2
        Fm[1L, pos] <<- fmax
      }
      if (kret == -1L) Fm[2L, pos] <<- -1
      pos <- point[pos]
    }
    pos <- new
    for (j in seq_len(maxi + maxi)) {
      if (Fm[1L, pos] < minf && Fm[2L, pos] == 0) {
        minf <<- Fm[1L, pos]
        minpos <<- pos
      }
      pos <- point[pos]
    }
  }
  insertlist_2 <- function(start, j, k) {
    pos <- start
    if (start == 0L) {
      list2[j, 1L] <<- 0L
      start <- j
    } else if (w[start] > w[j]) {
      list2[j, 1L] <<- start
      start <- j
    } else {
      for (i in seq_len(maxi)) {
        if (list2[pos, 1L] == 0L) {
          list2[j, 1L] <<- 0L
          list2[pos, 1L] <<- j
          break
        } else if (w[j] < w[list2[pos, 1L]]) {
          list2[j, 1L] <<- list2[pos, 1L]
          list2[pos, 1L] <<- j
          break
        }
        pos <- list2[pos, 1L]
      }
    }
    list2[j, 2L] <<- k
    start
  }
  divide <- function(new, currentlength, sample) {
    start <- 0L
    pos <- new
    for (i in seq_len(maxi)) {
      j <- arrayi[i]
      w[j] <<- Fm[1L, pos]
      k <- pos
      pos <- point[pos]
      d1 <- Fm[1L, pos]
      if (d1 < w[j]) w[j] <<- d1
      pos <- point[pos]
      start <- insertlist_2(start, j, k)
    }
    for (j in seq_len(maxi)) {
      # dirsearchmin
      k <- start
      pos <- list2[start, 2L]
      start <- list2[start, 1L]
      pos2 <- start
      Lm[k, sample] <<- currentlength + 1L
      for (i in seq_len(maxi - j + 1L)) {
        Lm[k, pos] <<- currentlength + 1L
        pos <- point[pos]
        Lm[k, pos] <<- currentlength + 1L
        if (pos2 > 0L) {
          pos <- list2[pos2, 2L]
          pos2 <- list2[pos2, 1L]
        }
      }
    }
  }
  # dirinsert: returns the advanced start
  dinsert <- function(start, ins) {
    for (i in seq_len(MAXFUNC)) {
      if (point[start] == 0L) {
        point[start] <<- ins
        point[ins] <<- 0L
        return(start)
      } else if (Fm[1L, ins] < Fm[1L, point[start]]) {
        help <- point[start]
        point[start] <<- ins
        point[ins] <<- help
        return(start)
      }
      start <- point[start]
    }
    start
  }
  insertlist <- function(new, samp) {
    for (j in seq_len(maxi)) {
      pos1 <- new
      pos2 <- point[pos1]
      new <- point[pos2]
      deep <- getlevel(pos1, jones)
      if (AN[deep + 2L] == 0L) {
        if (Fm[1L, pos2] < Fm[1L, pos1]) {
          AN[deep + 2L] <<- pos2
          point[pos2] <<- pos1
          point[pos1] <<- 0L
        } else {
          AN[deep + 2L] <<- pos1
          point[pos2] <<- 0L
        }
      } else {
        pos <- AN[deep + 2L]
        if (Fm[1L, pos2] < Fm[1L, pos1]) {
          if (Fm[1L, pos2] < Fm[1L, pos]) {
            AN[deep + 2L] <<- pos2
            if (Fm[1L, pos1] < Fm[1L, pos]) {
              point[pos2] <<- pos1
              point[pos1] <<- pos
            } else {
              point[pos2] <<- pos
              pos <- dinsert(pos, pos1)
            }
          } else {
            pos <- dinsert(pos, pos2)
            pos <- dinsert(pos, pos1)
          }
        } else {
          if (Fm[1L, pos1] < Fm[1L, pos]) {
            AN[deep + 2L] <<- pos1
            if (Fm[1L, pos] < Fm[1L, pos2]) {
              point[pos1] <<- pos
              pos <- dinsert(pos, pos2)
            } else {
              point[pos1] <<- pos2
              point[pos2] <<- pos
            }
          } else {
            pos <- dinsert(pos, pos1)
            pos <- dinsert(pos, pos2)
          }
        }
      }
    }
    deep <- getlevel(samp, jones)
    pos <- AN[deep + 2L]
    if (Fm[1L, samp] < Fm[1L, pos]) {
      AN[deep + 2L] <<- samp
      point[samp] <<- pos
    } else {
      dinsert(pos, samp)
    }
  }
  dir_choose <- function() {
    k <- 1L
    if (ifeasiblef >= 1L) {
      for (j in 0:MAXDEEP) {
        if (AN[j + 2L] > 0L) {
          S[k, 1L] <<- AN[j + 2L]
          S[k, 2L] <<- getlevel(S[k, 1L], jones)
          break
        }
      }
      return(1L)
    }
    for (j in which(AN[2L:(MAXDEEP + 2L)] > 0L) - 1L) {
      S[k, 1L] <<- AN[j + 2L]
      S[k, 2L] <<- getlevel(S[k, 1L], jones)
      k <- k + 1L
    }
    novalue <- 0L; novaluedeep <- 0L
    if (AN[1L] > 0L) {
      novalue <- AN[1L]
      novaluedeep <- getlevel(novalue, jones)
    }
    maxpos <- k - 1L
    if (k - 1L <= MAXDEEP && k <= MAXDIV) S[k, 1L] <<- 0L
    for (j in rev(seq_len(maxpos))) {
      helplower <- Inf
      helpgreater <- 0
      jj <- S[j, 1L]
      drop <- FALSE
      for (i in seq_len(j - 1L)) {
        ii <- S[i, 1L]
        if (ii > 0L && i != j && Fm[2L, ii] <= 1) {
          help2 <- levels[S[i, 2L] + 1L] - levels[S[j, 2L] + 1L]
          help2 <- (Fm[1L, ii] - Fm[1L, jj]) / help2
          if (isTRUE(help2 <= 0)) { drop <- TRUE; break }
          if (isTRUE(help2 < helplower)) helplower <- help2
        }
      }
      if (!drop && j + 1L <= maxpos) {
        for (i in (j + 1L):maxpos) {
          ii <- S[i, 1L]
          if (ii > 0L && i != j && Fm[2L, ii] <= 1) {
            help2 <- levels[S[i, 2L] + 1L] - levels[S[j, 2L] + 1L]
            help2 <- (Fm[1L, ii] - Fm[1L, jj]) / help2
            if (isTRUE(help2 <= 0)) { drop <- TRUE; break }
            if (isTRUE(help2 > helpgreater)) helpgreater <- help2
          }
        }
      }
      if (!drop) {
        if (helpgreater <= helplower) {
          if (cheat == 1L && helplower > kmax) helplower <- kmax
          a <- minf - eps * abs(minf); b <- minf - epsabs
          if (isTRUE(Fm[1L, jj] - helplower * levels[S[j, 2L] + 1L] > (if (isTRUE(a < b)) a else b))) drop <- TRUE
        } else drop <- TRUE
      }
      if (drop) S[j, 1L] <<- 0L
    }
    if (novalue > 0L) {
      maxpos <- maxpos + 1L
      S[maxpos, 1L] <<- novalue
      S[maxpos, 2L] <<- novaluedeep
    }
    maxpos
  }
  doubleinsert <- function(maxpos) {
    oldmaxpos <- maxpos
    for (i in seq_len(oldmaxpos)) {
      if (S[i, 1L] > 0L) {
        actdeep <- S[i, 2L]
        help <- AN[actdeep + 2L]
        pos <- point[help]
        iflag <- 0L
        while (pos > 0L && iflag == 0L) {
          if (isTRUE(Fm[1L, pos] - Fm[1L, help] <= 1e-13)) {
            if (maxpos < MAXDIV) {
              maxpos <- maxpos + 1L
              S[maxpos, 1L] <<- pos
              S[maxpos, 2L] <<- actdeep
              pos <- point[pos]
            } else {
              return(list(maxpos = maxpos, ierror = -6L))
            }
          } else iflag <- 1L
        }
      }
    }
    list(maxpos = maxpos, ierror = 0L)
  }
  resortlist <- function(replace) {
    lv <- getlevel(replace, jones)
    start <- AN[lv + 2L]
    if (replace != start) {
      pos <- start
      for (i in seq_len(MAXFUNC)) {
        if (point[pos] == replace) {
          point[pos] <<- point[replace]
          break
        } else pos <- point[pos]
        if (pos == 0L) break
      }
      if (Fm[1L, start] > Fm[1L, replace]) {
        AN[lv + 2L] <<- replace
        point[replace] <<- start
      } else {
        pos <- start
        for (i in seq_len(MAXFUNC)) {
          if (point[pos] == 0L) {
            point[replace] <<- point[pos]
            point[pos] <<- replace
            break
          } else {
            if (Fm[1L, point[pos]] > Fm[1L, replace]) {
              point[replace] <<- point[pos]
              point[pos] <<- replace
              break
            }
            pos <- point[pos]
          }
        }
      }
    }
  }
  tiny6 <- 0x1.0c6f7ap-20 # the float constant 1e-6f
  # dirreplaceinf: an undefined center (an infeasible point, or f NaN or
  # infinite) takes the smallest value of the defined centers in the box
  # around it, of half sides thirds(level of rectangle i in dimension j).
  # DIRsubrout.c reads that level as length[i__ + j * length_dim1], the
  # element (dimension i, rectangle j), while the same loop reads the center
  # as c__[j + i__ * c_dim1] and every other routine indexes length as
  # [dimension + rectangle * n]: the R port reads (dimension j, rectangle i),
  # and `upstream_replaceinf` reproduces the transposed reading of NLopt.
  replaceinf <- function() {
    m <- ifree - 1L
    feas <- which(Fm[2L, seq_len(m)] == 0)
    Cf <- Cm[, feas, drop = FALSE]
    for (i in which(Fm[2L, seq_len(m)] > 0)) {
      lv <- if (upstream_replaceinf) Lm[i + (seq_len(n) - 1L) * n] else Lm[, i]
      sidelength <- thirds[lv + 1L]
      a <- Cm[, i] - sidelength
      b <- Cm[, i] + sidelength
      Fm[1L, i] <<- Inf
      Fm[2L, i] <<- 2
      inbox <- feas[colSums(Cf < a | Cf > b) == 0]
      for (k in inbox) if (!(Fm[1L, i] < Fm[1L, k])) Fm[1L, i] <<- Fm[1L, k]
      if (length(inbox)) Fm[2L, i] <<- 1
      if (Fm[2L, i] == 1) {
        Fm[1L, i] <<- Fm[1L, i] + abs(Fm[1L, i]) * tiny6
        resortlist(i)
      } else if (!(fmax == Fm[1L, i])) {
        d1 <- fmax + 1
        Fm[1L, i] <<- if (d1 > Fm[1L, i]) d1 else Fm[1L, i]
      }
    }
  }
  # the exits of direct_direct_: cleanup leaves x as the scratch point, L100
  # sets it to the best center
  finish <- function(ierr) list(ierror = ierr, x = xs, minf = minf)
  l100 <- function(ierr) list(ierror = ierr, x = Cm[, minpos] * l + l * u, minf = minf)

  # ---- direct_direct_
  # dirheader
  if (eps < 0) {
    iepschange <- 1L
    epsfix <- -eps
    eps <- -eps
  } else {
    iepschange <- 0L
    epsfix <- 0x1.249ad2594c37dp+332 # 1e100
  }
  if (any(u <= l)) ierror <- -1L
  if (maxf + 20L > MAXFUNC) ierror <- -2L
  if (ierror < 0L) return(finish(ierror))
  divfactor <- if (fglobal == 0) 1 else abs(fglobal)
  oldmaxf <- maxf
  increase <- 0L
  # dirinitlist
  point <- c(seq_len(MAXFUNC - 1L) + 1L, 0L)
  ifree <- 1L
  # dirpreprc: l becomes the widths (c1), u the offsets l / width (c2)
  if (any(u <= l)) return(finish(-3L))
  help <- u - l
  u <- l / help
  l <- help
  # dirinit
  minf <- Inf
  if (jones == 0L) {
    for (j in 0:(n - 1L)) w[j + 1L] <- sqrt(n - j + j / 9) * 0.5
    help2 <- 1
    for (i in seq_len(MAXDEEP %/% n)) {
      for (j in 0:(n - 1L)) levels[(i - 1L) * n + j + 1L] <- w[j + 1L] / help2
      help2 <- help2 * 3
    }
  } else {
    help2 <- 3
    for (i in seq_len(MAXDEEP)) {
      levels[i + 1L] <- 1 / help2
      help2 <- help2 * 3
    }
    levels[1L] <- 1
  }
  help2 <- 3
  for (i in seq_len(MAXDEEP)) {
    thirds[i + 1L] <- 1 / help2
    help2 <- help2 * 3
  }
  thirds[1L] <- 1
  Cm[, 1L] <- 0.5
  xs <- rep(0.5, n)
  Lm[, 1L] <- 0L
  r <- infcn(xs)
  if (.nlopt_stop_forced(stop)) return(l100(-102L))
  Fm[1L, 1L] <- r$f
  Fm[2L, 1L] <- r$flag
  iinfesiblef <- r$flag
  fmax <- Fm[1L, 1L]
  if (Fm[2L, 1L] > 0) {
    Fm[1L, 1L] <- Inf
    fmax <- Inf
    ifeasiblef <- 1L
  } else ifeasiblef <- 0L
  minf <- Fm[1L, 1L]
  minpos <- 1L
  point[1L] <- 0L
  ifree <- 2L
  delta <- thirds[2L]
  if (.nlopt_stop_time(stop)) return(l100(6L))
  get_i(1L)
  sp <- samplepoints(1L, delta)
  if (sp$oops > 0L) return(finish(-4L))
  new <- sp$start
  samplef(new)
  if (.nlopt_stop_forced(stop)) return(l100(-102L))
  if (.nlopt_stop_time(stop)) return(l100(6L))
  divide(new, 0L, 1L)
  insertlist(new, 1L)

  numfunc <- maxi + 1L + maxi
  actmaxdeep <- 1L
  ierror <- 2L
  for (t in seq.int(2L, length.out = max(0L, maxt - 1L))) {
    actdeep <- actmaxdeep
    maxpos <- dir_choose()
    if (algmethod == 0L) {
      di <- doubleinsert(maxpos)
      maxpos <- di$maxpos
      if (di$ierror == -6L) return(finish(-6L))
    }
    oldpos <- minpos
    for (j in seq_len(maxpos)) {
      actdeep <- S[j, 2L]
      if (S[j, 1L] > 0L) {
        actdeep_div <- getmaxdeep(S[j, 1L])
        delta <- thirds[actdeep_div + 2L]
        actdeep <- S[j, 2L]
        if (actdeep + 1L >= mdeep) return(l100(-6L))
        if (actdeep > actmaxdeep) actmaxdeep <- actdeep
        help <- S[j, 1L]
        if (AN[actdeep + 2L] != help) {
          pos1 <- AN[actdeep + 2L]
          while (point[pos1] != help) pos1 <- point[pos1]
          point[pos1] <- point[help]
        } else {
          AN[actdeep + 2L] <- point[help]
        }
        get_i(help)
        sp <- samplepoints(help, delta)
        if (sp$oops > 0L) return(finish(-4L))
        start <- sp$start
        samplef(start)
        if (.nlopt_stop_forced(stop)) return(l100(-102L))
        if (.nlopt_stop_time(stop)) return(l100(6L))
        divide(start, actdeep_div, help)
        insertlist(start, help)
        numfunc <- numfunc + maxi + maxi
      }
    }
    actdeep_div <- getlevel(minpos, 0L)
    delta <- thirds[actdeep_div + 1L] * 100
    if (delta <= volper) { ierror <- 4L; break }
    actdeep_div <- getlevel(minpos, jones)
    delta <- levels[actdeep_div + 1L]
    if (delta <= sigmaper) { ierror <- 5L; break }
    # NaN (an unknown fglobal of -Inf) compares false, as in C
    if (isTRUE((minf - fglobal) * 100 / divfactor <= fglper)) { ierror <- 3L; break }
    if (iinfesiblef > 0L) replaceinf()
    ifreeold <- ifree
    if (iepschange == 1L) {
      d1 <- abs(minf) * 1e-4
      eps <- if (d1 > epsfix) d1 else epsfix
    }
    if (increase == 1L) {
      maxf <- numfunc + oldmaxf
      if (ifeasiblef == 0L) increase <- 0L
    }
    if (numfunc > maxf) {
      if (ifeasiblef == 0L) { ierror <- 1L; break }
      increase <- 1L
      maxf <- numfunc + oldmaxf
    }
  }
  l100(ierror)
}

# the GN_ORIG_DIRECT cases of nlopt_optimize_ (with f_direct)
.nlopt_orig_direct_dispatch <- function(o, x, stop) {
  n <- length(x)
  if (!.nlopt_finite_domain(o$lb, o$ub)) return(list(ret = "INVALID_ARGS", x = x, minf = Inf))
  fc <- o$fc
  fcn <- function(xx) {
    v <- o$f(xx)
    stop$nevals <- stop$nevals + 1L
    undefined <- is.nan(v) || .nlopt_isinf(v)
    if (!.nlopt_stop_forced(stop)) {
      for (cc in fc) {
        if (undefined) break
        cv <- cc$f(xx)
        if (.nlopt_stop_forced(stop)) break
        if (any(cv > 0)) undefined <- TRUE
      }
    }
    list(f = v, flag = as.integer(undefined))
  }
  r <- .nlopt_direct_optimize(fcn, n, o$lb, o$ub, x, stop$maxeval, -1L, stop,
                              .nlopt_get_param(o, "magic_eps", 0), .nlopt_get_param(o, "magic_eps_abs", 0),
                              stop$xtol_rel^n, .nlopt_get_param(o, "sigma_reltol", -1), stop$minf_max,
                              .nlopt_get_param(o, "fglobal_reltol", 0),
                              if (o$algorithm == "GN_ORIG_DIRECT") "ORIGINAL" else "GABLONSKY",
                              upstream_replaceinf = isTRUE(o$upstream_replaceinf))
  ret <- switch(as.character(r$ierror),
                "-1" = , "-2" = , "-101" = "INVALID_ARGS",
                "-3" = , "-4" = , "-5" = "FAILURE",
                "1" = , "2" = "MAXEVAL_REACHED",
                "6" = "MAXTIME_REACHED",
                "3" = "STOPVAL_REACHED",
                "4" = , "5" = "XTOL_REACHED",
                "-100" = "OUT_OF_MEMORY",
                "-102" = "FORCED_STOP",
                "SUCCESS")
  if (ret %in% c("INVALID_ARGS", "FAILURE")) return(list(ret = ret, x = x, minf = Inf))
  list(ret = ret, x = r$x, minf = r$minf)
}
