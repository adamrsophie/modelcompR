# DIRECT of SciPy 1.18.1 (BSD 3-clause license): the driver direct() of
# scipy/optimize/_direct_py.py and the C code of scipy/optimize/_direct
# (DIRect.c, DIRsubrout.c, DIRserial.c, direct_wrap.c), the translation by
# Steven G. Johnson of the DIRECT 2.0.4 Fortran code of Joerg M. Gablonsky
# (North Carolina State University, MIT license). SciPy compiles it with
# clang, which fuses a * b + c within one expression (-ffp-contract=on);
# .spdr_fma() stands at each such site.
#
# Deliberate difference: DIRChoose of DIRsubrout.c drops a rectangle off the
# lower convex hull (helpgreater > helplower) only when a log file is open,
# because the `goto L60` stands under the `if (logfile)` of a commented out
# fprintf. direct() never opens the log (disp = FALSE), so SciPy divides
# rectangles that DIRECT does not select. Here they are dropped, as SciPy
# does with its log on (tools/scipy-probe.py direct, the drf_ cases).
#
# SciPy's objective wrapper (direct_dirinfcn_) never flags a point as
# infeasible and nothing sets force_stop, so the branches of DIRECT for
# infeasible points (DIRreplaceInf, DIRresortlist, the level -1 of points
# without a value, the increased budget) cannot be reached; they are left
# out, and with them the point flags (always 0).
#
# The points are kept as in the C code: cc and ln are the n x maxfunc
# arrays of centers and side length exponents, fv the function values,
# point the linked lists (0 ends a list); anchor[d + 2] heads the list of
# level d (d = -1, ..., maxdeep), thirds[d + 1] and levels[d + 1] hold the
# side lengths of level d and S the potentially optimal rectangles.

.spdr_error_messages <- c("Number of function evaluations done is larger than maxfun={}",
                          "Number of iterations is larger than maxiter={}", "u[i] < l[i] for some i",
                          "maxfun is too large", "Initialization failed",
                          "There was an error in the creation of the sample points",
                          "An error occurred while the function was sampled",
                          "Maximum number of levels has been reached.", "Forced stop", "Invalid arguments",
                          "Out of memory")

.spdr_success_messages <- c(paste("The best function value found is within a relative error={}",
                                  "of the (known) global optimum f_min"),
                            paste("The volume of the hyperrectangle containing the lowest function value",
                                  "found is below vol_tol={}"),
                            paste("The side length measure of the hyperrectangle containing the lowest",
                                  "function value found is below len_tol={}"))

# the message of _direct_py.py for a return code; vals are the formatted
# maxfun, maxiter, f_min_rtol, vol_tol and len_tol
.spdr_message <- function(rc, vals) {
  if (rc > 2L) return(sub("{}", vals[rc], .spdr_success_messages[rc - 2L], fixed = TRUE))
  if (rc > 0L) return(sub("{}", vals[rc], .spdr_error_messages[rc], fixed = TRUE))
  if (rc > -100L) return(.spdr_error_messages[abs(rc) + 2L])
  .spdr_error_messages[length(.spdr_error_messages) + rc + 100L]
}

# a fused multiply add of scalars that also passes infinities through
.spdr_fma <- function(a, b, c) {
  p <- a * b
  if (is.finite(p) && is.finite(c)) .fma(a, b, c) else p + c
}

# direct_direct_ of DIRect.c with its subroutines: fcn is the objective of
# the unscaled point, l and u the bounds, jones 1 for DIRECT_L (locally
# biased, Gablonsky) and 0 for the original DIRECT of Jones
.spdr_direct <- function(fcn, l, u, eps, maxf, maxt, jones, fglobal, fglper, volper, sigmaper, callback) {
  n <- length(l)
  maxdiv <- 5000L
  memmax <- 1073741824
  maxfunc <- if (maxf <= 0L) 101000L else maxf + 1000L + maxf %/% 2L
  fmd <- n * 12 + 20
  if (maxfunc * fmd > memmax) maxfunc <- as.integer(memmax %/% fmd)
  maxf <- if (maxf <= 0L) maxfunc - 1000L else (2L * (maxfunc - 1000L)) %/% 3L
  maxdeep <- if (maxt <= 0L) maxfunc %/% 5L else maxt + 1000L
  if (maxdeep * 20 + 24 > memmax) maxdeep <- as.integer((memmax - 24) %/% 20)
  maxt <- if (maxt <= 0L) maxdeep else maxdeep - 1000L
  # DIRheader
  ierror <- 0L
  if (eps < 0) {
    iepschange <- 1L
    epsfix <- -eps
    eps <- -eps
  } else {
    iepschange <- 0L
    epsfix <- 0x1.249ad2594c37dp+332 # 1e100
  }
  if (any(!(u > l))) ierror <- -1L
  if (maxf + 20L > maxfunc) ierror <- -2L
  # direct_optimize returns NULL without a Python error: _directmodule.c
  # raises a RuntimeError
  if (ierror < 0L) stop(sprintf("DIRECT optimization failed (code %d)", ierror), call. = FALSE)
  divfactor <- if (fglobal == 0) 1 else abs(fglobal)
  # DIRInitList and DIRpreprc
  anchor <- integer(maxdeep + 2L)
  fv <- numeric(maxfunc)
  point <- c(seq_len(maxfunc - 1L) + 1L, 0L)
  ifree <- 1L
  xs1 <- u - l
  xs2 <- l / xs1
  cc <- matrix(0, n, maxfunc)
  ln <- matrix(0L, n, maxfunc)
  thirds <- numeric(maxdeep + 1L)
  levels <- numeric(maxdeep + 1L)
  w <- numeric(n)
  list2 <- matrix(0L, n, 2L)
  arrayi <- integer(n)
  maxi <- 0L
  S <- matrix(0L, maxdiv, 2L)
  maxpos <- 0L
  minf <- Inf
  minpos <- 1L
  fmax <- 0
  xseq <- numeric(n)

  # DIRinfcn: the objective at the unscaled point (x + xs2) * xs1
  infcn <- function(x) {
    xseq <<- (x + xs2) * xs1
    v <- as.numeric(fcn(xseq))
    if (length(v) != 1L) stop("can only convert an array of size 1 to a Python scalar", call. = FALSE)
    v
  }
  getlevel <- function(pos, jones) {
    v <- ln[, pos]
    if (jones == 0L) {
      help <- v[1L]
      k <- min(v)
      p <- sum(v == help)
      if (k == help) k * n + n - p else k * n + p
    } else {
      min(v)
    }
  }
  get_i <- function(pos) {
    v <- ln[, pos]
    a <- which(v == min(v))
    maxi <<- length(a)
    arrayi[seq_len(maxi)] <<- a
  }
  # DIRChoose: the potentially optimal rectangles (the lower right convex
  # hull of the level heads), epsabs = 0 as direct_optimize passes it; a
  # rectangle with helpgreater > helplower is dropped (see the header)
  choose <- function() {
    k <- 1L
    for (li in which(anchor[seq_len(maxdeep + 1L) + 1L] > 0L)) {
      S[k, 1L] <<- anchor[li + 1L]
      S[k, 2L] <<- getlevel(S[k, 1L], jones)
      k <- k + 1L
    }
    maxpos <<- k - 1L
    if (k - 1L <= maxdeep) S[k, 1L] <<- 0L
    epsabs <- 0
    for (j in rev(seq_len(maxpos))) {
      helplower <- Inf
      helpgreater <- 0
      jj <- S[j, 1L]
      lj <- levels[S[j, 2L] + 1L]
      fj <- fv[jj]
      drop <- FALSE
      for (i in seq_len(j - 1L)) {
        ii <- S[i, 1L]
        if (ii > 0L) {
          h <- (fv[ii] - fj) / (levels[S[i, 2L] + 1L] - lj)
          if (isTRUE(h <= 0)) {
            drop <- TRUE
            break
          }
          if (isTRUE(h < helplower)) helplower <- h
        }
      }
      if (!drop) {
        for (i in seq_len(maxpos - j) + j) {
          ii <- S[i, 1L]
          if (ii > 0L) {
            h <- (fv[ii] - fj) / (levels[S[i, 2L] + 1L] - lj)
            if (isTRUE(h <= 0)) {
              drop <- TRUE
              break
            }
            if (isTRUE(h > helpgreater)) helpgreater <- h
          }
        }
      }
      if (!drop) {
        if (isTRUE(helpgreater <= helplower)) {
          lo <- .spdr_fma(-eps, abs(minf), minf)
          hi <- minf - epsabs
          if (isTRUE(.spdr_fma(-helplower, lj, fj) > (if (isTRUE(lo < hi)) lo else hi))) drop <- TRUE
        } else {
          drop <- TRUE
        }
      }
      if (drop) S[j, 1L] <<- 0L
    }
  }
  # DIRDoubleInsert (original DIRECT): the rectangles of the same level
  # whose value ties with the level head are divided too
  double_insert <- function() {
    oldmaxpos <- maxpos
    for (i in seq_len(oldmaxpos)) {
      if (S[i, 1L] > 0L) {
        actdeep <- S[i, 2L]
        help <- anchor[actdeep + 2L]
        pos <- point[help]
        iflag <- 0L
        while (pos > 0L && iflag == 0L) {
          if (isTRUE(fv[pos] - fv[help] <= 1e-13)) {
            if (maxpos < maxdiv) {
              maxpos <<- maxpos + 1L
              S[maxpos, 1L] <<- pos
              S[maxpos, 2L] <<- actdeep
              pos <- point[pos]
            } else {
              ierror <<- -6L
              return(invisible())
            }
          } else {
            iflag <- 1L
          }
        }
      }
    }
  }
  # DIRSamplepoints: 2 maxi new points from the free list (NA when it runs
  # out), returning the first
  sample_points <- function(sample, delta) {
    pos <- ifree
    start <- ifree
    for (k in seq_len(2L * maxi)) {
      ln[, ifree] <<- ln[, sample]
      cc[, ifree] <<- cc[, sample]
      pos <- ifree
      ifree <<- point[ifree]
      if (ifree == 0L) return(NA_integer_)
    }
    point[pos] <<- 0L
    pos <- start
    for (j in seq_len(maxi)) {
      a <- arrayi[j]
      cc[a, pos] <<- cc[a, sample] + delta
      pos <- point[pos]
      cc[a, pos] <<- cc[a, sample] - delta
      pos <- point[pos]
    }
    start
  }
  # DIRSamplef (DIRserial.c)
  sample_f <- function(start) {
    pos <- start
    for (j in seq_len(2L * maxi)) {
      v <- infcn(cc[, pos])
      fv[pos] <<- v
      if (isTRUE(v > fmax)) fmax <<- v
      pos <- point[pos]
    }
    pos <- start
    for (j in seq_len(2L * maxi)) {
      if (isTRUE(fv[pos] < minf)) {
        minf <<- fv[pos]
        minpos <<- pos
      }
      pos <- point[pos]
    }
  }
  # DIRInsertList_2: dimension j into the list sorted by w
  insert_list2 <- function(start, j, k) {
    pos <- start
    if (start == 0L) {
      list2[j, 1L] <<- 0L
      start <- j
    } else if (isTRUE(w[start] > w[j])) {
      list2[j, 1L] <<- start
      start <- j
    } else {
      for (i in seq_len(maxi)) {
        if (list2[pos, 1L] == 0L) {
          list2[j, 1L] <<- 0L
          list2[pos, 1L] <<- j
          break
        }
        if (isTRUE(w[j] < w[list2[pos, 1L]])) {
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
  # DIRDivide: divide along the dimensions in the order of their best value
  divide <- function(new, currentlength, sample) {
    start <- 0L
    pos <- new
    for (i in seq_len(maxi)) {
      j <- arrayi[i]
      w[j] <<- fv[pos]
      k <- pos
      pos <- point[pos]
      if (isTRUE(fv[pos] < w[j])) w[j] <<- fv[pos]
      pos <- point[pos]
      start <- insert_list2(start, j, k)
    }
    for (j in seq_len(maxi)) {
      k <- start
      pos <- list2[start, 2L]
      start <- list2[start, 1L]
      pos2 <- start
      ln[k, sample] <<- currentlength + 1L
      for (i in seq_len(maxi - j + 1L)) {
        ln[k, pos] <<- currentlength + 1L
        pos <- point[pos]
        ln[k, pos] <<- currentlength + 1L
        if (pos2 > 0L) {
          pos <- list2[pos2, 2L]
          pos2 <- list2[pos2, 1L]
        }
      }
    }
  }
  # DIRInsert: ins into the sorted list after start, returning the new start
  insert <- function(start, ins) {
    for (i in seq_len(maxfunc)) {
      if (point[start] == 0L) {
        point[start] <<- ins
        point[ins] <<- 0L
        return(start)
      }
      if (isTRUE(fv[ins] < fv[point[start]])) {
        help <- point[start]
        point[start] <<- ins
        point[ins] <<- help
        return(start)
      }
      start <- point[start]
    }
    start
  }
  # DIRInsertList: the new pairs and the divided rectangle into the lists
  # of their levels
  insert_list <- function(new, samp) {
    for (j in seq_len(maxi)) {
      pos1 <- new
      pos2 <- point[pos1]
      new <- point[pos2]
      ai <- getlevel(pos1, jones) + 2L
      if (anchor[ai] == 0L) {
        if (isTRUE(fv[pos2] < fv[pos1])) {
          anchor[ai] <<- pos2
          point[pos2] <<- pos1
          point[pos1] <<- 0L
        } else {
          anchor[ai] <<- pos1
          point[pos2] <<- 0L
        }
      } else {
        pos <- anchor[ai]
        if (isTRUE(fv[pos2] < fv[pos1])) {
          if (isTRUE(fv[pos2] < fv[pos])) {
            anchor[ai] <<- pos2
            if (isTRUE(fv[pos1] < fv[pos])) {
              point[pos2] <<- pos1
              point[pos1] <<- pos
            } else {
              point[pos2] <<- pos
              pos <- insert(pos, pos1)
            }
          } else {
            pos <- insert(pos, pos2)
            pos <- insert(pos, pos1)
          }
        } else {
          if (isTRUE(fv[pos1] < fv[pos])) {
            anchor[ai] <<- pos1
            if (isTRUE(fv[pos] < fv[pos2])) {
              point[pos1] <<- pos
              pos <- insert(pos, pos2)
            } else {
              point[pos1] <<- pos2
              point[pos2] <<- pos
            }
          } else {
            pos <- insert(pos, pos1)
            pos <- insert(pos, pos2)
          }
        }
      }
    }
    ai <- getlevel(samp, jones) + 2L
    pos <- anchor[ai]
    if (isTRUE(fv[samp] < fv[pos])) {
      anchor[ai] <<- samp
      point[samp] <<- pos
    } else {
      insert(pos, samp)
    }
  }
  # the result: the center of the best rectangle unscaled (L100), or the
  # last evaluated point when DIRECT stops at its cleanup label
  done <- function(t, numfunc, at_best = TRUE) {
    x <- if (at_best) vapply(seq_len(n), function(i) .fma(cc[i, minpos], xs1[i], xs1[i] * xs2[i]), 0) else xseq
    list(x = x, fun = minf, status = ierror, nfev = numfunc, nit = t)
  }

  # DIRInit: the center and the first division
  if (jones == 0L) {
    for (j in 0:(n - 1L)) w[j + 1L] <- sqrt((n - j) + j / 9) * 0.5
    help2 <- 1
    for (i in seq_len(maxdeep %/% n)) {
      for (j in 0:(n - 1L)) levels[(i - 1L) * n + j + 1L] <- w[j + 1L] / help2
      help2 <- help2 * 3
    }
  } else {
    help2 <- 3
    for (i in seq_len(maxdeep)) {
      levels[i + 1L] <- 1 / help2
      help2 <- help2 * 3
    }
    levels[1L] <- 1
  }
  help2 <- 3
  for (i in seq_len(maxdeep)) {
    thirds[i + 1L] <- 1 / help2
    help2 <- help2 * 3
  }
  thirds[1L] <- 1
  cc[, 1L] <- 0.5
  fv[1L] <- infcn(rep(0.5, n))
  fmax <- fv[1L]
  minf <- fv[1L]
  minpos <- 1L
  point[1L] <- 0L
  ifree <- 2L
  get_i(1L)
  new <- sample_points(1L, thirds[2L])
  if (is.na(new)) {
    ierror <- -4L
    return(done(0L, 1L, FALSE))
  }
  sample_f(new)
  divide(new, 0L, 1L)
  insert_list(new, 1L)
  numfunc <- maxi + 1L + maxi
  actmaxdeep <- 1L
  t <- 2L
  while (t <= maxt - 1L) {
    choose()
    if (jones == 0L) {
      double_insert()
      if (ierror == -6L) return(done(t, numfunc, FALSE))
    }
    for (j in seq_len(maxpos)) {
      help <- S[j, 1L]
      if (help > 0L) {
        actdeep_div <- min(ln[, help])
        delta <- thirds[actdeep_div + 2L]
        actdeep <- S[j, 2L]
        if (actdeep + 1L >= maxdeep) {
          ierror <- -6L
          return(done(t, numfunc))
        }
        actmaxdeep <- max(actdeep, actmaxdeep)
        # take the rectangle out of the list of its level
        if (anchor[actdeep + 2L] != help) {
          pos1 <- anchor[actdeep + 2L]
          while (point[pos1] != help) pos1 <- point[pos1]
          point[pos1] <- point[help]
        } else {
          anchor[actdeep + 2L] <- point[help]
        }
        get_i(help)
        start <- sample_points(help, delta)
        if (is.na(start)) {
          ierror <- -4L
          return(done(t, numfunc, FALSE))
        }
        sample_f(start)
        divide(start, actdeep_div, help)
        insert_list(start, help)
        numfunc <- numfunc + maxi + maxi
      }
    }
    # the stopping rules
    if (isTRUE(thirds[getlevel(minpos, 0L) + 1L] <= volper)) {
      ierror <- 4L
      return(done(t, numfunc))
    }
    if (isTRUE(levels[getlevel(minpos, jones) + 1L] <= sigmaper)) {
      ierror <- 5L
      return(done(t, numfunc))
    }
    if (isTRUE((minf - fglobal) / divfactor <= fglper)) {
      ierror <- 3L
      return(done(t, numfunc))
    }
    if (iepschange == 1L) {
      e <- abs(minf) * 1e-4
      eps <- if (isTRUE(e > epsfix)) e else epsfix
    }
    if (numfunc > maxf) {
      ierror <- 1L
      return(done(t, numfunc))
    }
    if (!is.null(callback)) callback(xseq)
    t <- t + 1L
  }
  ierror <- 2L
  done(t, numfunc)
}

# scipy.optimize.direct: the checks of _direct_py.py and direct_optimize
# (direct_wrap.c). maxfun and maxiter must be whole numbers (SciPy asks for
# a Python int; an R double such as 2000 is accepted). The callback gets the
# last evaluated point after each iteration.
.sp_direct <- function(func, lb, ub, eps = 1e-4, maxfun = NULL, maxiter = 1000L, locally_biased = TRUE,
                       f_min = -Inf, f_min_rtol = 1e-4, vol_tol = 1e-16, len_tol = 1e-6, callback = NULL) {
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (!isTRUE(all(lb < ub))) stop("Bounds are not consistent min < max", call. = FALSE)
  if (any(is.infinite(lb)) || any(is.infinite(ub))) stop("Bounds must not be inf.", call. = FALSE)
  if (vol_tol < 0 || vol_tol > 1) stop("vol_tol must be between 0 and 1.", call. = FALSE)
  if (len_tol < 0 || len_tol > 1) stop("len_tol must be between 0 and 1.", call. = FALSE)
  if (f_min_rtol < 0 || f_min_rtol > 1) stop("f_min_rtol must be between 0 and 1.", call. = FALSE)
  n <- length(lb)
  if (is.null(maxfun)) maxfun <- 1000L * n
  whole <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v)
  if (!whole(maxfun)) stop("maxfun must be of type int.", call. = FALSE)
  if (maxfun < 0) stop("maxfun must be > 0.", call. = FALSE)
  if (!whole(maxiter)) stop("maxiter must be of type int.", call. = FALSE)
  if (maxiter < 0) stop("maxiter must be > 0.", call. = FALSE)
  if (!isTRUE(locally_biased) && !isFALSE(locally_biased)) stop("locally_biased must be True or False.", call. = FALSE)
  if (n < 1L) stop(sprintf("DIRECT requires at least one dimension, got %d", n), call. = FALSE)
  # direct_optimize: tolerances <= 0 are switched off, and f_min_rtol is
  # ignored when f_min is unknown
  volper <- if (vol_tol <= 0) -1 else vol_tol
  sigmaper <- if (len_tol <= 0) -1 else len_tol
  fglper <- if (f_min == -Inf) 0 else f_min_rtol
  r <- .spdr_direct(func, lb, ub, eps, as.integer(maxfun), as.integer(maxiter), if (locally_biased) 1L else 0L,
                    f_min, fglper, volper, sigmaper, callback)
  vals <- c(as.character(as.integer(maxfun)), as.character(as.integer(maxiter)), .spsq_repr(f_min_rtol),
            .spsq_repr(vol_tol), .spsq_repr(len_tol))
  list(x = r$x, fun = r$fun, status = r$status, success = r$status > 2L, message = .spdr_message(r$status, vals),
       nfev = r$nfev, nit = r$nit)
}

# ScipyDirect: the options of OptimizationSciPy.jl, then direct(); the
# solution is Success exactly when SciPy reports success
.scipy_run_direct <- function(f, lb, ub, opts, kw) {
  args <- utils::modifyList(list(eps = 0.0001, maxiter = as.integer(opts$maxiters %||% 1000L), locally_biased = TRUE,
                                 vol_tol = 1e-16, len_tol = 1e-6), kw)
  res <- do.call(.sp_direct, c(list(function(x) as.numeric(f(x)), lb, ub), args))
  res$retcode <- if (res$success) "Success" else "Failure"
  res
}
