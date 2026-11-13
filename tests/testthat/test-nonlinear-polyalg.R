# Polyalgorithms, multistart, Jacobian reuse, tracing, Sobol points, chunk
# sizes, sparsity detection and colorings. Reference values printed by
# tools/nonlinearsolve-polyalg-probe.jl (NonlinearSolve 4.32.0,
# NonlinearSolveFirstOrder 2.10.0, Sobol 1.5.0, SparseConnectivityTracer
# 1.2.4, SparseMatrixColorings 0.4.29).

# long reference comparisons: run on CI (NOT_CRAN = true), not on CRAN
skip_on_cran()

.pa_defs <- function() {
  e <- new.env()
  local({
    rosen <- function(u, p) c(10 * (u[2] - u[1] * u[1]), 1 - u[1])
    rosenj <- function(u, p) matrix(c(-20 * u[1], -1, 10, 0), 2)
    s5 <- sqrt(5); s10 <- sqrt(10)
    powell <- function(u, p) c(u[1] + 10 * u[2], s5 * (u[3] - u[4]), (u[2] - 2 * u[3]) * (u[2] - 2 * u[3]),
                               s10 * (u[1] - u[4]) * (u[1] - u[4]))
    powellj <- function(u, p) rbind(c(1, 10, 0, 0), c(0, 0, s5, -s5), c(0, 2 * (u[2] - 2 * u[3]), -4 * (u[2] - 2 * u[3]), 0),
                                    c(2 * s10 * (u[1] - u[4]), 0, 0, -2 * s10 * (u[1] - u[4])))
    tt <- c(1, 2, 4, 8, 16); yy <- c(0.9, 1.5, 2.1, 2.6, 2.75)
    mm <- function(u, p) u[1] * tt / (u[2] + tt) - yy
    mmj <- function(u, p) cbind(tt / (u[2] + tt), -u[1] * tt / ((u[2] + tt) * (u[2] + tt)))
    circ <- function(u, p) c(u[1] * u[1] + u[2] * u[2] - 4, u[1] * u[2] - 1)
    circj <- function(u, p) rbind(c(2 * u[1], 2 * u[2]), c(u[2], u[1]))
    broyden <- function(u, p) {
      n <- length(u)
      (3 - 2 * u) * u - c(0, u[-n]) - 2 * c(u[-1], 0) + 1
    }
    broydenj <- function(u, p) {
      n <- length(u)
      J <- diag(3 - 4 * u, n)
      for (i in seq_len(n)) {
        if (i > 1) J[i, i - 1] <- -1
        if (i < n) J[i, i + 1] <- -2
      }
      J
    }
    outside <- function(u, p) c(u[1] - 3, u[2] - 0.5)
    outsidej <- function(u, p) diag(2)
    sq <- function(u, p) c(u[1] * u[1] - 1, u[2] - 1)
    sqj <- function(u, p) rbind(c(2 * u[1], 0), c(0, 1))
    nlp <- function(f, j, u0, lb, ub) NonlinearProblem(NonlinearFunction(f, jac = j), u0 = u0, lb = lb, ub = ub)
    nlls <- function(f, j, u0, lb, ub) NonlinearLeastSquaresProblem(NonlinearFunction(f, jac = j), u0 = u0, lb = lb, ub = ub)
    cases <- list(
      rosen_in = nlp(rosen, rosenj, c(-1.2, 1), c(-2, -2), c(2, 2)),
      rosen_bd = nlp(rosen, rosenj, c(-1.2, 1), c(-2, -2), c(0.5, 2)),
      powell = nlp(powell, powellj, c(3, -1, 0, 1), rep(-5, 4), rep(5, 4)),
      circ = nlp(circ, circj, c(3, 0.5), c(0, 0), c(3, 3)),
      mm_ls = nlls(mm, mmj, c(1, 1), c(0, 0), c(2.5, 10)),
      outside = nlp(outside, outsidej, c(1, 1), c(0, 0), c(2, 2)),
      sq = nlp(sq, sqj, c(0, 0), c(-2, -2), c(2, 2)),
      rosen_half = nlp(rosen, rosenj, c(-1.2, 1), c(-2, -Inf), c(2, Inf)))
  }, envir = e)
  e
}

.pa_alg <- function(tag) switch(tag,
  FSBP = FastShortcutBoundedPolyalg(),
  FSBP_post = FastShortcutBoundedPolyalg(must_support_postcondition = TRUE),
  Default = SobolMultistart(FastShortcutBoundedPolyalg()),
  Sobol4 = SobolMultistart(nstarts = 4),
  SobolAll = SobolMultistart(nstarts = 6, early_exit = FALSE),
  SobolNoRest = SobolMultistart(nstarts = 3, restoration = FALSE),
  Poly_DT = NonlinearSolvePolyAlgorithm(list(Dogbox(), TrustRegionReflective())),
  Poly_BLMTR = NonlinearSolvePolyAlgorithm(list(BoundedLevenbergMarquardt(), BoundedTrustRegion()), start_index = 2),
  BTR = BoundedTrustRegion(),
  GFOA_auto = BoundedTrustRegion(),
  GFOA_off = BoundedTrustRegion(jacobian_reuse = FALSE),
  GFOA_on = BoundedTrustRegion(jacobian_reuse = TRUE),
  GFOA_age3 = BoundedTrustRegion(jacobian_reuse = JacobianReuse(max_age = 3)),
  GFOA_ratio = BoundedTrustRegion(jacobian_reuse = JacobianReuse(max_age = 10, max_residual_ratio = 0.5)))

.pa_cases <- list(
  list("rosen_in", "FSBP", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_in", "FSBP_post", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_in", "Default", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_in", "Sobol4", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_in", "SobolAll", "Success", c(1.0, 1.0),
       c(77, 66, 60, 77, 77)),
  list("rosen_in", "SobolNoRest", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_in", "Poly_DT", "Success", c(1.0, 1.0),
       c(17, 18, 16, 16, 16)),
  list("rosen_in", "Poly_BLMTR", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_bd", "FSBP", "Stalled", c(0.5, 0.24870073444322444),
       c(83, 30, 26, 63, 64)),
  list("rosen_bd", "FSBP_post", "Stalled", c(0.5, 0.24870073444322444),
       c(55, 18, 17, 54, 54)),
  list("rosen_bd", "Sobol4", "Stalled", c(0.5, 0.24870073444322444),
       c(321, 153, 133, 276, 280)),
  list("rosen_bd", "SobolAll", "Stalled", c(0.5, 0.24870073444322444),
       c(430, 196, 166, 378, 384)),
  list("rosen_bd", "SobolNoRest", "Stalled", c(0.5, 0.24870073444322444),
       c(221, 101, 89, 189, 192)),
  list("rosen_bd", "Poly_DT", "Stalled", c(0.5, 0.25),
       c(33, 36, 31, 31, 32)),
  list("rosen_bd", "Poly_BLMTR", "Stalled", c(0.5, 0.24870073444322444),
       c(55, 18, 17, 54, 54)),
  list("powell", "FSBP", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("powell", "FSBP_post", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("powell", "Default", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("powell", "Sobol4", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("powell", "SobolAll", "Success", c(0.0, 0.0, 0.0, 0.0),
       c(121, 127, 121, 121, 121)),
  list("powell", "SobolNoRest", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("powell", "Poly_DT", "Success", c(2.974414806415404e-7, -2.9744148064154046e-8, 6.444461612574505e-8, 6.444461612574503e-8),
       c(24, 25, 23, 23, 23)),
  list("powell", "Poly_BLMTR", "Success", c(2.838316417875744e-7, -2.838316417875744e-8, 4.541306268601189e-8, 4.541306268601189e-8),
       c(23, 24, 23, 23, 23)),
  list("circ", "FSBP", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("circ", "FSBP_post", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("circ", "Default", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("circ", "Sobol4", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("circ", "SobolNoRest", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("circ", "Poly_DT", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(6, 7, 5, 5, 5)),
  list("circ", "Poly_BLMTR", "Success", c(1.9318516525781366, 0.5176380902050415),
       c(5, 6, 5, 5, 5)),
  list("mm_ls", "FSBP", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("mm_ls", "FSBP_post", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("mm_ls", "Default", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("mm_ls", "Sobol4", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("mm_ls", "SobolAll", "Success", c(2.5, 1.1174912615968227),
       c(66, 78, 66, 66, 72)),
  list("mm_ls", "SobolNoRest", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("mm_ls", "Poly_DT", "Success", c(2.5, 1.1174912512338095),
       c(11, 13, 10, 10, 11)),
  list("mm_ls", "Poly_BLMTR", "Success", c(2.5, 1.1174912631931133),
       c(11, 13, 11, 11, 12)),
  list("outside", "FSBP", "Stalled", c(2.0, 0.6570028297149824),
       c(37, 8, 4, 35, 36)),
  list("outside", "FSBP_post", "Stalled", c(2.0, 0.6570028297149824),
       c(35, 4, 3, 34, 34)),
  list("outside", "Default", "Stalled", c(2.0, 0.5),
       c(592, 155, 75, 571, 587)),
  list("outside", "Sobol4", "Stalled", c(2.0, 0.5),
       c(149, 39, 19, 143, 147)),
  list("outside", "SobolAll", "Stalled", c(2.0, 0.5),
       c(222, 58, 28, 214, 220)),
  list("outside", "SobolNoRest", "Stalled", c(2.0, 0.5),
       c(109, 23, 11, 104, 107)),
  list("outside", "Poly_DT", "Stalled", c(2.0, 0.5),
       c(11, 14, 9, 9, 10)),
  list("outside", "Poly_BLMTR", "Stalled", c(2.0, 0.6570028297149824),
       c(35, 4, 3, 34, 34)),
  list("sq", "FSBP", "Stalled", c(0.0, 0.0),
       c(36, 7, 3, 34, 35)),
  list("sq", "FSBP_post", "Stalled", c(0.0, 0.0),
       c(34, 3, 2, 33, 33)),
  list("sq", "Default", "Success", c(1.0, 1.0),
       c(74, 17, 8, 70, 72)),
  list("sq", "Sobol4", "Success", c(1.0, 1.0),
       c(74, 17, 8, 70, 72)),
  list("sq", "SobolAll", "Success", c(1.0, 1.0),
       c(85, 31, 19, 81, 83)),
  list("sq", "SobolNoRest", "Success", c(1.0, 1.0),
       c(74, 17, 8, 70, 72)),
  list("sq", "Poly_DT", "Stalled", c(0.0, 1.0),
       c(18, 21, 16, 16, 17)),
  list("sq", "Poly_BLMTR", "Stalled", c(0.0, 0.0),
       c(34, 3, 2, 33, 33)),
  list("rosen_half", "FSBP", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_half", "FSBP_post", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_half", "Default", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_half", "Sobol4", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_half", "SobolAll", "Success", c(1.0, 1.0),
       c(85, 70, 64, 85, 85)),
  list("rosen_half", "SobolNoRest", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("rosen_half", "Poly_DT", "Success", c(1.0, 1.0),
       c(17, 18, 16, 16, 16)),
  list("rosen_half", "Poly_BLMTR", "Success", c(1.0, 1.0),
       c(26, 21, 20, 26, 26)),
  list("broyden8", "BTR", "Success", c(-0.570469836964552, -0.6811405903331907, -0.7004284708169576, -0.7006724537898226, -0.6917363327761061, -0.6657674263517837, -0.5960292391307003, -0.4164109994188799),
       c(5, 6, 5, 5, 5)),
  list("broyden8", "GFOA_auto", "Success", c(-0.570469836964552, -0.6811405903331907, -0.7004284708169576, -0.7006724537898226, -0.6917363327761061, -0.6657674263517837, -0.5960292391307003, -0.4164109994188799),
       c(5, 6, 5, 5, 5)),
  list("broyden8", "GFOA_off", "Success", c(-0.570469836964552, -0.6811405903331907, -0.7004284708169576, -0.7006724537898226, -0.6917363327761061, -0.6657674263517837, -0.5960292391307003, -0.4164109994188799),
       c(5, 6, 5, 5, 5)),
  list("broyden8", "GFOA_on", "Success", c(-0.5704698369645521, -0.6811405903331907, -0.7004284708169577, -0.7006724537898226, -0.6917363327761061, -0.6657674263517838, -0.5960292391307003, -0.41641099941887993),
       c(7, 4, 3, 7, 7)),
  list("broyden8", "GFOA_age3", "Success", c(-0.5704698369645521, -0.6811405903331907, -0.7004284708169577, -0.7006724537898226, -0.6917363327761061, -0.6657674263517837, -0.5960292391307002, -0.4164109994188799),
       c(7, 5, 4, 7, 7)),
  list("broyden8", "GFOA_ratio", "Success", c(-0.5704698369645521, -0.6811405903331909, -0.7004284708169579, -0.7006724537898231, -0.6917363327761075, -0.6657674263517861, -0.5960292391307039, -0.416410999418883),
       c(12, 3, 2, 12, 12)),
  list("broyden16", "BTR", "Success", c(-0.5707610988550857, -0.6819098802488937, -0.7024853557268581, -0.7062587684765307, -0.7069469229013529, -0.7070649519134553, -0.7070648126438835, -0.7070053922882276, -0.7068323068350298, -0.7063576740941611, -0.7050615214754398, -0.701525194231341, -0.6918946287506123, -0.6657975233041955, -0.596035312618988, -0.4164123011651986),
       c(7, 4, 3, 7, 7)),
  list("broyden16", "GFOA_auto", "Success", c(-0.5707610988550857, -0.6819098802488937, -0.7024853557268581, -0.7062587684765307, -0.7069469229013529, -0.7070649519134553, -0.7070648126438835, -0.7070053922882276, -0.7068323068350298, -0.7063576740941611, -0.7050615214754398, -0.701525194231341, -0.6918946287506123, -0.6657975233041955, -0.596035312618988, -0.4164123011651986),
       c(7, 4, 3, 7, 7)),
  list("broyden16", "GFOA_off", "Success", c(-0.5707610988550857, -0.6819098802488937, -0.7024853557268581, -0.7062587684765307, -0.7069469229013529, -0.7070649519134553, -0.7070648126438834, -0.7070053922882275, -0.70683230683503, -0.7063576740941611, -0.7050615214754398, -0.7015251942313411, -0.6918946287506124, -0.6657975233041955, -0.5960353126189879, -0.4164123011651986),
       c(5, 6, 5, 5, 5)),
  list("broyden16", "GFOA_on", "Success", c(-0.5707610988550857, -0.6819098802488937, -0.7024853557268581, -0.7062587684765307, -0.7069469229013529, -0.7070649519134553, -0.7070648126438835, -0.7070053922882276, -0.7068323068350298, -0.7063576740941611, -0.7050615214754398, -0.701525194231341, -0.6918946287506123, -0.6657975233041955, -0.596035312618988, -0.4164123011651986),
       c(7, 4, 3, 7, 7)),
  list("broyden16", "GFOA_age3", "Success", c(-0.5707610988550857, -0.6819098802488937, -0.7024853557268581, -0.7062587684765307, -0.7069469229013529, -0.7070649519134553, -0.7070648126438835, -0.7070053922882276, -0.7068323068350298, -0.7063576740941611, -0.7050615214754398, -0.701525194231341, -0.6918946287506123, -0.6657975233041955, -0.596035312618988, -0.4164123011651986),
       c(7, 4, 3, 7, 7)),
  list("broyden16", "GFOA_ratio", "Success", c(-0.5707610988550859, -0.6819098802488937, -0.702485355726858, -0.7062587684765308, -0.7069469229013529, -0.7070649519134553, -0.7070648126438835, -0.7070053922882275, -0.70683230683503, -0.7063576740941611, -0.7050615214754399, -0.7015251942313415, -0.6918946287506135, -0.665797523304198, -0.5960353126189917, -0.4164123011652017),
       c(12, 3, 2, 12, 12)),
  list("broyden20", "BTR", "Success", c(-0.5707611912831241, -0.6819101243996171, -0.7024860087165649, -0.7062605433175683, -0.7069517656652966, -0.7070781758164413, -0.7071009276079201, -0.7071040253276404, -0.7071016768220528, -0.707093283933818, -0.7070699996740117, -0.7070063419831152, -0.7068324807420124, -0.7063577059533686, -0.705061527318752, -0.7015251953064953, -0.6918946289501835, -0.6657975233421397, -0.5960353126266449, -0.41641230116683975),
       c(7, 4, 3, 7, 7)),
  list("broyden20", "GFOA_auto", "Success", c(-0.5707611912831241, -0.6819101243996171, -0.7024860087165649, -0.7062605433175683, -0.7069517656652966, -0.7070781758164413, -0.7071009276079201, -0.7071040253276404, -0.7071016768220528, -0.707093283933818, -0.7070699996740117, -0.7070063419831152, -0.7068324807420124, -0.7063577059533686, -0.705061527318752, -0.7015251953064953, -0.6918946289501835, -0.6657975233421397, -0.5960353126266449, -0.41641230116683975),
       c(7, 4, 3, 7, 7)),
  list("broyden20", "GFOA_off", "Success", c(-0.5707611912831241, -0.681910124399617, -0.7024860087165647, -0.7062605433175682, -0.7069517656652966, -0.7070781758164414, -0.7071009276079201, -0.7071040253276404, -0.7071016768220528, -0.707093283933818, -0.7070699996740115, -0.7070063419831152, -0.7068324807420124, -0.7063577059533686, -0.7050615273187522, -0.7015251953064953, -0.6918946289501835, -0.6657975233421397, -0.5960353126266449, -0.4164123011668397),
       c(5, 6, 5, 5, 5)),
  list("broyden20", "GFOA_on", "Success", c(-0.5707611912831241, -0.6819101243996171, -0.7024860087165649, -0.7062605433175683, -0.7069517656652966, -0.7070781758164413, -0.7071009276079201, -0.7071040253276404, -0.7071016768220528, -0.707093283933818, -0.7070699996740117, -0.7070063419831152, -0.7068324807420124, -0.7063577059533686, -0.705061527318752, -0.7015251953064953, -0.6918946289501835, -0.6657975233421397, -0.5960353126266449, -0.41641230116683975),
       c(7, 4, 3, 7, 7)),
  list("broyden20", "GFOA_age3", "Success", c(-0.5707611912831241, -0.6819101243996171, -0.7024860087165649, -0.7062605433175683, -0.7069517656652966, -0.7070781758164413, -0.7071009276079201, -0.7071040253276404, -0.7071016768220528, -0.707093283933818, -0.7070699996740117, -0.7070063419831152, -0.7068324807420124, -0.7063577059533686, -0.705061527318752, -0.7015251953064953, -0.6918946289501835, -0.6657975233421397, -0.5960353126266449, -0.41641230116683975),
       c(7, 4, 3, 7, 7)),
  list("broyden20", "GFOA_ratio", "Success", c(-0.5707611912831241, -0.6819101243996171, -0.7024860087165647, -0.7062605433175682, -0.7069517656652967, -0.7070781758164413, -0.70710092760792, -0.7071040253276404, -0.7071016768220528, -0.7070932839338181, -0.7070699996740115, -0.7070063419831151, -0.7068324807420124, -0.7063577059533687, -0.7050615273187524, -0.7015251953064959, -0.6918946289501848, -0.6657975233421422, -0.5960353126266487, -0.4164123011668429),
       c(12, 3, 2, 12, 12)))

test_that("bounded polyalgorithms, multistart and Jacobian reuse match NonlinearSolve", {
  e <- .pa_defs()
  for (cs in .pa_cases) {
    if (startsWith(cs[[1]], "broyden")) {
      n <- as.integer(sub("broyden", "", cs[[1]]))
      prob <- if (cs[[2]] == "BTR") e$nlp(e$broyden, e$broydenj, rep(-1, n), rep(-2, n), rep(2, n))
              else NonlinearProblem(NonlinearFunction(e$broyden, jac = e$broydenj), u0 = rep(-1, n))
    } else prob <- e$cases[[cs[[1]]]]
    # SobolMultistart warns, as NonlinearSolve does, when restoration leaves the box
    sol <- suppressWarnings(solve(prob, .pa_alg(cs[[2]])))
    lab <- paste(cs[[1]], cs[[2]])
    expect_identical(sol$retcode, cs[[3]], label = lab)
    expect_lte(max(abs(unname(sol$u) - cs[[4]])), ref_tol(1e-12) * max(1, abs(cs[[4]])), label = lab)
    st <- sol$stats
    expect_equal(c(st$nf, st$njacs, st$nfactors, st$nsolve, st$nsteps), cs[[5]], label = lab)
  }
})

test_that("trace levels record the NonlinearSolve history", {
  e <- .pa_defs()
  prob <- e$cases$rosen_in
  algs <- list(Dogbox = Dogbox, BTR = BoundedTrustRegion, TRR = TrustRegionReflective)
  lvls <- list(min = TraceMinimal(), cond = TraceWithJacobianConditionNumber(), min2 = TraceMinimal(2))
  ref <- list(
    "Dogbox min" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
      c(4.3999999999999995, 3.5716100288074926, 3.5716100288074926, 1.3476294862879084, 1.3476294862879084, 1.2226294862879084, 0.9726294862879084, 0.7226294862879084, 0.6416343804011332, 0.3294344516343789, 0.3294344516343789, 0.23806284735354943, 0.2130352481392861, 0.07399342922247276, 0.036747308824715, 0.0018865837085746406, 0.0),
      c(0.0, 1.167411767882168, 0.0, 0.25833162369689905, 0.0, 0.12732162345637751, 0.2666515835202896, 0.258889350402022, 0.31926802653015923, 0.3066278877249605, 0.0, 0.12589734939569214, 0.26707630422421036, 0.05937048965627525, 0.12739393636316917, 0.03369463812490503, 0.00018865837085746406),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    "Dogbox cond" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
      c(4.3999999999999995, 3.5716100288074926, 3.5716100288074926, 1.3476294862879084, 1.3476294862879084, 1.2226294862879084, 0.9726294862879084, 0.7226294862879084, 0.6416343804011332, 0.3294344516343789, 0.3294344516343789, 0.23806284735354943, 0.2130352481392861, 0.07399342922247276, 0.036747308824715, 0.0018865837085746406, 0.0),
      c(0.0, 1.167411767882168, 0.0, 0.25833162369689905, 0.0, 0.12732162345637751, 0.2666515835202896, 0.258889350402022, 0.31926802653015923, 0.3066278877249605, 0.0, 0.12589734939569214, 0.26707630422421036, 0.05937048965627525, 0.12739393636316917, 0.03369463812490503, 0.00018865837085746406),
      c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)),
    "Dogbox min2" = list(c(0, 1, 3, 5, 7, 9, 11, 13, 15),
      c(4.3999999999999995, 3.5716100288074926, 1.3476294862879084, 1.2226294862879084, 0.7226294862879084, 0.3294344516343789, 0.23806284735354943, 0.07399342922247276, 0.0018865837085746406),
      c(0.0, 1.167411767882168, 0.25833162369689905, 0.12732162345637751, 0.258889350402022, 0.3066278877249605, 0.12589734939569214, 0.05937048965627525, 0.03369463812490503),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    "BTR min" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26),
      c(4.3999999999999995, 4.3999999999999995, 1.9065030721167389, 1.9065030721167389, 1.7775861739077303, 1.5585326008950808, 1.5585326008950808, 1.4064596332729495, 1.4064596332729495, 1.3200797881602133, 1.1462166280882888, 1.1462166280882888, 1.0489439796890678, 0.8590285964747715, 0.8654940090454308, 0.7765297044166288, 0.5940850864614783, 0.5865299953167833, 0.5069273427568355, 0.4461149395809836, 0.31770073082652783, 0.21018091376957304, 0.21018091376957304, 0.1634015670660871, 0.10051089249740297, 0.041476408399405695, 0.0),
      c(0.0, 0.0, 0.3905124837953327, 0.0, 0.19525624189766633, 0.39051248379533277, 0.0, 0.19525624189766638, 0.0, 0.09762812094883319, 0.19525624189766633, 0.0, 0.09762812094883319, 0.19525624189766633, 0.04881406047441659, 0.09762812094883318, 0.19525624189766635, 0.04881406047441659, 0.09762812094883318, 0.09762812094883319, 0.19525624189766633, 0.19525624189766635, 0.0, 0.09762812094883318, 0.19525624189766633, 0.1455802129527021, 0.0041476408399405695),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    "BTR cond" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26),
      c(4.3999999999999995, 4.3999999999999995, 1.9065030721167389, 1.9065030721167389, 1.7775861739077303, 1.5585326008950808, 1.5585326008950808, 1.4064596332729495, 1.4064596332729495, 1.3200797881602133, 1.1462166280882888, 1.1462166280882888, 1.0489439796890678, 0.8590285964747715, 0.8654940090454308, 0.7765297044166288, 0.5940850864614783, 0.5865299953167833, 0.5069273427568355, 0.4461149395809836, 0.31770073082652783, 0.21018091376957304, 0.21018091376957304, 0.1634015670660871, 0.10051089249740297, 0.041476408399405695, 0.0),
      c(0.0, 0.0, 0.3905124837953327, 0.0, 0.19525624189766633, 0.39051248379533277, 0.0, 0.19525624189766638, 0.0, 0.09762812094883319, 0.19525624189766633, 0.0, 0.09762812094883319, 0.19525624189766633, 0.04881406047441659, 0.09762812094883318, 0.19525624189766635, 0.04881406047441659, 0.09762812094883318, 0.09762812094883319, 0.19525624189766633, 0.19525624189766635, 0.0, 0.09762812094883318, 0.19525624189766633, 0.1455802129527021, 0.0041476408399405695),
      c(Inf, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)),
    "BTR min2" = list(c(0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25),
      c(4.3999999999999995, 4.3999999999999995, 1.9065030721167389, 1.5585326008950808, 1.4064596332729495, 1.3200797881602133, 1.1462166280882888, 0.8590285964747715, 0.7765297044166288, 0.5865299953167833, 0.4461149395809836, 0.21018091376957304, 0.1634015670660871, 0.041476408399405695),
      c(0.0, 0.0, 0.0, 0.39051248379533277, 0.19525624189766638, 0.09762812094883319, 0.0, 0.19525624189766633, 0.09762812094883318, 0.04881406047441659, 0.09762812094883319, 0.19525624189766635, 0.09762812094883318, 0.1455802129527021),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    "TRR min" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27),
      c(4.3999999999999995, 2.0364085830119647, 1.928239743591777, 1.928239743591777, 1.785883908987288, 1.785883908987288, 1.6895539678554514, 1.4882899808754164, 1.4882899808754164, 1.3897132835314348, 1.2090686944463744, 1.2090686944463744, 1.098122789576256, 0.8810568444846869, 0.6802556076331552, 0.5143880647875033, 0.38906047092491924, 0.2806620934262556, 0.20420295324406157, 0.1276873344102325, 0.08613367174569997, 0.03632859378874553, 0.019424882302374824, 0.002936049489061787, 0.000328866542312678, 1.1197851290667415e-6, 5.492473142965082e-11, 0.0),
      c(0.0, 0.16686752794224297, 0.20614965495174012, 0.0, 0.2894296802741775, 0.0, 0.16255673777271137, 0.3359714946362332, 0.0, 0.11401128711287686, 0.22319526560675354, 0.0, 0.11238114294027739, 0.2192562963863475, 0.21967402249906262, 0.22141903548801636, 0.19611252316641917, 0.18438837058111165, 0.14424984389553872, 0.1491927427907025, 0.08889699898017626, 0.10549403054235143, 0.038997308559630885, 0.037267799577973236, 0.0054235568156032095, 0.0007389850013251913, 2.3581585102040983e-6, 1.239417942299626e-10),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    "TRR cond" = list(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27),
      c(4.3999999999999995, 2.0364085830119647, 1.928239743591777, 1.928239743591777, 1.785883908987288, 1.785883908987288, 1.6895539678554514, 1.4882899808754164, 1.4882899808754164, 1.3897132835314348, 1.2090686944463744, 1.2090686944463744, 1.098122789576256, 0.8810568444846869, 0.6802556076331552, 0.5143880647875033, 0.38906047092491924, 0.2806620934262556, 0.20420295324406157, 0.1276873344102325, 0.08613367174569997, 0.03632859378874553, 0.019424882302374824, 0.002936049489061787, 0.000328866542312678, 1.1197851290667415e-6, 5.492473142965082e-11, 0.0),
      c(0.0, 0.16686752794224297, 0.20614965495174012, 0.0, 0.2894296802741775, 0.0, 0.16255673777271137, 0.3359714946362332, 0.0, 0.11401128711287686, 0.22319526560675354, 0.0, 0.11238114294027739, 0.2192562963863475, 0.21967402249906262, 0.22141903548801636, 0.19611252316641917, 0.18438837058111165, 0.14424984389553872, 0.1491927427907025, 0.08889699898017626, 0.10549403054235143, 0.038997308559630885, 0.037267799577973236, 0.0054235568156032095, 0.0007389850013251913, 2.3581585102040983e-6, 1.239417942299626e-10),
      c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)),
    "TRR min2" = list(c(0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27),
      c(4.3999999999999995, 2.0364085830119647, 1.928239743591777, 1.785883908987288, 1.4882899808754164, 1.3897132835314348, 1.2090686944463744, 0.8810568444846869, 0.5143880647875033, 0.2806620934262556, 0.1276873344102325, 0.03632859378874553, 0.002936049489061787, 1.1197851290667415e-6, 0.0),
      c(0.0, 0.16686752794224297, 0.0, 0.0, 0.3359714946362332, 0.11401128711287686, 0.0, 0.2192562963863475, 0.22141903548801636, 0.18438837058111165, 0.1491927427907025, 0.10549403054235143, 0.037267799577973236, 0.0007389850013251913, 1.239417942299626e-10),
      c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
  for (tag in names(algs)) for (lt in names(lvls)) {
    sol <- solve(prob, algs[[tag]](), store_trace = TRUE, trace_level = lvls[[lt]])
    h <- sol$trace$history
    r <- ref[[paste(tag, lt)]]
    lab <- paste(tag, lt)
    # NonlinearSolve traces cond(J) of the Jacobian buffer before the first evaluation
    # (uninitialized memory from safe_similar when the Jacobian is analytic; it printed
    # Inf in the probe); the port reports the condition number of J(u0)
    if (tag == "BTR" && lt == "cond") r[[4]][1] <- kappa(e$rosenj(c(-1.2, 1)), exact = TRUE)
    expect_identical(vapply(h, function(x) as.numeric(x$iteration), 0), r[[1]], label = lab)
    expect_equal(vapply(h, function(x) x$fnorm, 0), r[[2]], tolerance = 1e-12, label = lab)
    expect_equal(vapply(h, function(x) x$stepnorm, 0), r[[3]], tolerance = 1e-12, label = lab)
    expect_equal(vapply(h, function(x) if (is.null(x$condJ)) NA_real_ else x$condJ, 0), r[[4]],
                 tolerance = 1e-10, label = lab)
  }
  sol <- solve(e$cases$mm_ls, Dogbox(), store_trace = TRUE, trace_level = TraceMinimal())
  h <- sol$trace$history
  expect_equal(vapply(h, function(x) x$fnorm, 0),
               c(2.9570573615616556, 1.2629951259005638, 0.7787963211455952, 0.6683954141039944, 0.6683949469096117, 0.66839493400864, 0.6683949336544645, 0.6683949336447317, 0.6683949336444647, 0.6683949336444572, 0.6683949336444568), tolerance = 1e-12)
  expect_equal(vapply(h, function(x) x$stepnorm, 0),
               c(0.0, 1.0114693808982835, 0.5088063084176969, 0.3644658109932226, 0.0009692649107928553, 0.0001609463246499132, 2.6670600479894446e-5, 4.421116839736783e-6, 7.328359392655415e-7, 1.2147462769007689e-7, 2.013556033020336e-8), tolerance = 1e-12)
  expect_true(all(vapply(h, function(x) x$norm_type, "") == "L2"))
})

test_that("the Sobol sequence matches Sobol.jl", {
  s <- .sobol_new(1)
  ref <- list(c(0.5), c(0.75), c(0.25), c(0.375), c(0.875), c(0.625), c(0.125), c(0.1875), c(0.6875), c(0.9375), c(0.4375), c(0.3125))
  for (k in seq_along(ref)) expect_identical(.sobol_next(s), ref[[k]])
  s <- .sobol_new(2)
  ref <- list(c(0.5, 0.5), c(0.75, 0.25), c(0.25, 0.75), c(0.375, 0.375), c(0.875, 0.875), c(0.625, 0.125), c(0.125, 0.625), c(0.1875, 0.3125), c(0.6875, 0.8125), c(0.9375, 0.0625), c(0.4375, 0.5625), c(0.3125, 0.1875))
  for (k in seq_along(ref)) expect_identical(.sobol_next(s), ref[[k]])
  s <- .sobol_new(5)
  ref <- list(c(0.5, 0.5, 0.5, 0.5, 0.5), c(0.75, 0.25, 0.25, 0.25, 0.75), c(0.25, 0.75, 0.75, 0.75, 0.25), c(0.375, 0.375, 0.625, 0.875, 0.375), c(0.875, 0.875, 0.125, 0.375, 0.875), c(0.625, 0.125, 0.875, 0.625, 0.625), c(0.125, 0.625, 0.375, 0.125, 0.125), c(0.1875, 0.3125, 0.9375, 0.4375, 0.5625), c(0.6875, 0.8125, 0.4375, 0.9375, 0.0625), c(0.9375, 0.0625, 0.6875, 0.1875, 0.3125), c(0.4375, 0.5625, 0.1875, 0.6875, 0.8125), c(0.3125, 0.1875, 0.3125, 0.5625, 0.9375))
  for (k in seq_along(ref)) expect_identical(.sobol_next(s), ref[[k]])
  s <- .sobol_new(3); lo <- c(-1, 0, 2); hi <- c(1, 10, 3)
  ref <- list(c(0.0, 5.0, 2.5), c(0.5, 2.5, 2.25), c(-0.5, 7.5, 2.75), c(-0.25, 3.75, 2.625), c(0.75, 8.75, 2.125), c(0.25, 1.25, 2.875))
  for (k in seq_along(ref)) expect_equal(lo + (hi - lo) * .sobol_next(s), ref[[k]])
  s <- .sobol_new(40)
  for (k in 1:40) x <- .sobol_next(s)
  expect_identical(x, c(0.234375, 0.078125, 0.265625, 0.984375, 0.703125, 0.734375, 0.734375, 0.359375, 0.046875, 0.359375, 0.828125, 0.890625, 0.515625, 0.046875, 0.171875, 0.921875, 0.984375, 0.859375, 0.609375, 0.828125, 0.390625, 0.015625, 0.359375, 0.203125, 0.359375, 0.171875, 0.734375, 0.453125, 0.578125, 0.296875, 0.796875, 0.328125, 0.921875, 0.796875, 0.421875, 0.921875, 0.953125, 0.109375, 0.265625, 0.890625))
})

test_that("pickchunksize matches ForwardDiff", {
  expect_identical(as.numeric(sapply(c(1, 5, 12, 13, 24, 25, 30, 100, 1000), pickchunksize)), c(1, 5, 12, 7, 12, 9, 10, 12, 12))
  expect_identical(as.numeric(pickchunksize(numeric(30))), 10)
})

.pa_banded12 <- function() {
  n <- 12
  S <- diag(TRUE, n)
  for (i in 1:n) {
    if (i < n) S[i, i + 1] <- S[i + 1, i] <- TRUE
    if (i + 3 <= n) S[i, i + 3] <- S[i + 3, i] <- TRUE
  }
  S
}
.pa_arrow12 <- function() {
  A <- diag(TRUE, 12); A[1, ] <- TRUE; A[, 1] <- TRUE
  A
}
test_that("sparsity detection and colorings match SparseConnectivityTracer and SparseMatrixColorings", {
  spf1 <- function(x) c(x[1] * x[2], sin(x[3]), x[1] + x[4] * x[4], abs(x[2]) + floor(x[3]), x[4] / x[2])
  spf2 <- function(x) c(sum(x[1:3]), prod(x[2:4]), max(x[1], x[4]), x[3]^2, exp(x[1]) * x[4])
  hsf1 <- function(x) x[1] * x[2] + exp(x[3]) + x[4] / x[1] + x[2]^3
  hsf2 <- function(x) sum(x * x) + x[1] * x[3] * x[5] + sin(x[2] + x[4])
  pat <- function(S) {
    S <- as(S, "TsparseMatrix")
    o <- order(S@j, S@i)
    list(dim(S), S@i[o] + 1L, S@j[o] + 1L)
  }
  expect_equal(pat(jacobian_sparsity(spf1, c(0.3, 0.7, 0.2, 0.9))),
               list(c(5L, 4L), c(1L, 3L, 1L, 4L, 5L, 2L, 3L, 5L),
                    c(1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L)))
  expect_equal(pat(jacobian_sparsity(spf2, c(0.3, 0.7, 0.2, 0.9))),
               list(c(5L, 4L), c(1L, 3L, 5L, 1L, 2L, 1L, 2L, 4L, 2L, 3L, 5L),
                    c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L)))
  expect_equal(pat(hessian_sparsity(hsf1, c(0.3, 0.7, 0.2, 0.9))),
               list(c(4L, 4L), c(1L, 2L, 4L, 1L, 2L, 3L, 1L),
                    c(1L, 1L, 1L, 2L, 2L, 3L, 4L)))
  expect_equal(pat(hessian_sparsity(hsf2, c(0.3, 0.7, 0.2, 0.9, 0.5))),
               list(c(5L, 5L), c(1L, 3L, 5L, 2L, 4L, 1L, 3L, 5L, 2L, 4L, 1L, 3L, 5L),
                    c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 5L, 5L, 5L)))
  pt <- function(m, n, a, b, c) outer(seq_len(m), seq_len(n), function(i, j) (a * i + b * j) %% c == 0 | i == j)
  expect_equal(column_coloring(pt(6, 6, 3, 5, 7)), c(1, 1, 2, 1, 2, 2))
  expect_equal(row_coloring(pt(6, 6, 3, 5, 7)), c(1, 1, 2, 1, 2, 2))
  expect_equal(symmetric_coloring(pt(6, 6, 3, 5, 7) | t(pt(6, 6, 3, 5, 7))), c(1, 1, 2, 1, 3, 4))
  expect_equal(column_coloring(pt(8, 5, 2, 3, 5)), c(1, 1, 1, 1, 1))
  expect_equal(row_coloring(pt(8, 5, 2, 3, 5)), c(1, 1, 1, 1, 1, 2, 2, 2))
  expect_equal(column_coloring(pt(5, 9, 4, 1, 6)), c(1, 2, 1, 1, 2, 2, 1, 3, 1))
  expect_equal(row_coloring(pt(5, 9, 4, 1, 6)), c(1, 2, 1, 3, 1))
  expect_equal(column_coloring(pt(10, 10, 1, 1, 4)), c(1, 1, 2, 1, 3, 2, 4, 2, 5, 3))
  expect_equal(row_coloring(pt(10, 10, 1, 1, 4)), c(1, 1, 2, 1, 3, 2, 4, 2, 5, 3))
  expect_equal(symmetric_coloring(pt(10, 10, 1, 1, 4) | t(pt(10, 10, 1, 1, 4))), c(1, 1, 2, 1, 1, 2, 3, 2, 1, 3))
  expect_equal(symmetric_coloring(.pa_banded12()), c(1, 2, 1, 3, 1, 4, 1, 5, 1, 2, 1, 3))
  expect_equal(symmetric_coloring(.pa_arrow12()), c(1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
})

test_that("ADTypes constructors return their backends and replacements", {
  exact <- c("AutoChainRules", "AutoDiffractor", "AutoEnzyme", "AutoFastDifferentiation", "AutoGTPSA",
             "AutoHyperHessians", "AutoMooncake", "AutoMooncakeForward", "AutoPolyesterForwardDiff",
             "AutoReactant", "AutoReverseDiff", "AutoSymbolics", "AutoTapir", "AutoTaylorDiff", "AutoTracker",
             "AutoZygote")
  for (a in exact) {
    o <- get(a)()
    expect_identical(class(o), c(a, "AutoForwardDiff", "ADType"), label = a)
  }
  expect_identical(class(AutoFiniteDifferences()), c("AutoFiniteDifferences", "AutoFiniteDiff", "ADType"))
  sparse <- c(AutoSparseFastDifferentiation = "AutoFastDifferentiation", AutoSparseFiniteDiff = "AutoFiniteDiff",
              AutoSparseForwardDiff = "AutoForwardDiff", AutoSparsePolyesterForwardDiff = "AutoPolyesterForwardDiff",
              AutoSparseReverseDiff = "AutoReverseDiff", AutoSparseZygote = "AutoZygote")
  for (a in names(sparse)) {
    o <- get(a)()
    expect_s3_class(o, "AutoSparse")
    expect_identical(class(o$dense_ad)[1L], sparse[[a]], label = a)
  }
  expect_identical(class(AutoModelingToolkit())[1L], "AutoSymbolics")
  expect_s3_class(AutoModelingToolkit(obj_sparse = TRUE), "AutoSparse")
  expect_identical(class(NoAutoDiff()), c("NoAutoDiff", "ADType"))
  f <- function(u, p) c(u[1] * u[1] - 2, u[2] - 1)
  prob <- NonlinearProblem(f, u0 = c(1, 0), lb = c(0, 0), ub = c(2, 2))
  expect_error(solve(prob, TrustRegionReflective(autodiff = NoAutoDiff())), class = "NoAutoDiffSelectedError")
  expect_equal(solve(prob, TrustRegionReflective(autodiff = AutoZygote()))$u, c(sqrt(2), 1), tolerance = 1e-10)
  e <- NoAutoDiffSelectedError()
  expect_s3_class(e, "error")
  expect_identical(conditionMessage(e), "Automatic differentiation can not be used with NoAutoDiff()")
})

test_that("NonlinearVerbosity resolves presets, groups and toggles", {
  v <- NonlinearVerbosity()
  expect_s3_class(v, "NonlinearVerbosity")
  expect_identical(v$toggles[["forcing"]], "InfoLevel")
  expect_identical(v$toggles[["unsupported_postcondition"]], "ErrorLevel")
  expect_identical(v$linear_verbosity, "None")
  expect_true(all(NonlinearVerbosity("None")$toggles == "Silent"))
  expect_identical(NonlinearVerbosity(FALSE), NonlinearVerbosity("None"))
  expect_identical(NonlinearVerbosity(TRUE), NonlinearVerbosity())
  w <- NonlinearVerbosity(numerical = "Silent", termination_condition = "ErrorLevel")
  expect_identical(w$toggles[c("threshold_state", "forcing", "termination_condition")],
                   c(threshold_state = "Silent", forcing = "Silent", termination_condition = "ErrorLevel"))
  expect_error(NonlinearVerbosity("Loud"), "unknown verbosity preset")
  expect_error(NonlinearVerbosity(bogus = "Silent"), "unknown verbosity toggle")
  expect_error(NonlinearVerbosity(forcing = "Loud"), "unknown verbosity level")
})

test_that("PostconditionSpecifier wraps a corrector", {
  H <- function(u, u_prev, p, cache) pmax(u, 0)
  post <- PostconditionSpecifier(H, space = PostconditionSpace$Transformed)
  expect_identical(attr(post, "space"), "Transformed")
  expect_identical(attr(PostconditionSpecifier(H), "space"), "Original")
  expect_identical(post(c(-1, 2), c(0, 0), NULL), c(0, 2))
  expect_error(PostconditionSpecifier(H, space = "Physical"))
  sol <- solve(NonlinearProblem(function(u, p) u * u - 2, u0 = 1), BoundedTrustRegion(), postcondition = post)
  expect_equal(sol$u, sqrt(2), tolerance = 1e-12)
})

test_that("EisenstatWalkerForcing2 reproduces the forcing terms of NonlinearSolve", {
  # tools/nonlinearsolve-forcing-probe.jl: NewtonRaphson(linsolve = KrylovJL_GMRES(),
  # forcing = ...). Julia calls the post step hook before the step is applied
  # (NonlinearSolveFirstOrder src/solve.jl), so it records the residual norm at
  # the start of each step; replaying those norms gives the same forcing terms.
  cases <- list(
    circ = list(EisenstatWalkerForcing2(),
      c(5.058665832015394, 0.8790280286177874, 0.2810531689830291, 0.12359120748106263, 0.038689858912543534, 0.00047059843517766704, 3.0348254974231134e-8),
      c(0.5, 0.9, 0.7290000000000001, 0.47829690000000014, 0.2058911320946491, 0.08819859174652676, 0.00013315245399621694)),
    broyden = list(EisenstatWalkerForcing2(),
      c(4.123105625617661, 1.3687401654612674, 0.5397781154138267, 0.23765115071544923, 0.1108657467134924, 0.008786987938037404, 0.00030636634706792364, 1.0494612981218367e-6),
      c(0.5, 0.9, 0.7290000000000001, 0.47829690000000014, 0.2058911320946491, 0.19586527169395776, 0.005653635392209312, 0.0010940686295993528)),
    circ_nosg = list(EisenstatWalkerForcing2(eta_0 = 0.3, eta_max = 0.6, gamma = 0.8, alpha = 1.5, safeguard = FALSE),
      c(5.058665832015394, 0.8790280286177874, 0.2810531689830291, 0.022148136135603173, 5.772009598188011e-5, 6.538476422956662e-10),
      c(0.3, 0.6, 0.057948217898269476, 0.14463341211622524, 0.017697545341327312, 0.00010643262597460791)),
    circ_sg = list(EisenstatWalkerForcing2(eta_0 = 0.8, eta_max = 0.95, gamma = 0.9, alpha = 1.2, safeguard_threshold = 0.05),
      c(5.058665832015394, 0.8790280286177874, 0.2810531689830291, 0.12359120748106263, 0.038689858912543534, 0.018945273230970124, 0.006482533665872703, 1.1868179684279175e-5, 2.446779188684926e-6),
      c(0.8, 0.9, 0.7931101735123912, 0.6814637116711687, 0.5680332927108219, 0.4565528500681578, 0.3820559626056429, 0.2836579028788634, 0.19842499350743414)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    cc <- .ew_init(cs[[1]], 1)
    eta <- numeric(0)
    for (k in seq_along(cs[[3]])) {
      eta[k] <- .ew_pre(cc, cs[[2]][k], k - 1L)
      .ew_post(cc, cs[[2]][k])
    }
    expect_equal(eta, cs[[3]], tolerance = 1e-15, label = nm)
  }
  expect_s3_class(EisenstatWalkerForcing2(), "EisenstatWalkerForcing2")
})

test_that("TraceAll stores the state of each step", {
  # tools/nonlinearsolve-forcing-probe.jl; the problem has no Jacobian (ForwardDiff
  # in Julia, differences in R), hence the tolerance
  f <- function(u, p) c(u[1] * u[1] - 2, u[2] - 1)
  prob <- NonlinearProblem(f, u0 = c(1, 0), lb = c(0, 0), ub = c(2, 2))
  sol <- solve(prob, Dogbox(), store_trace = TRUE, trace_level = TraceAll())
  h <- sol$trace$history
  expect_identical(vapply(h, function(e) e$iteration, 0L), 0:5)
  expect_true(all(vapply(h, function(e) is.null(e$storage$J) && e$condJ == -1, TRUE)))
  expect_equal(h[[2]]$storage$u, c(1.5, 1), tolerance = 1e-9)
  expect_equal(h[[2]]$storage$fu, c(0.25, 0), tolerance = 1e-9)
  expect_equal(h[[2]]$storage$du, c(0.5, 1), tolerance = 1e-9)
  expect_equal(h[[3]]$storage$u, c(1.4166666666666667, 1), tolerance = 1e-9)
  sol <- solve(prob, BoundedTrustRegion(), store_trace = TRUE, trace_level = TraceAll())
  h <- sol$trace$history
  expect_identical(vapply(h, function(e) e$iteration, 0L), 0:5)
  expect_equal(h[[1]]$storage$J, diag(c(2, 1)), tolerance = 1e-9)
  expect_equal(h[[1]]$condJ, 2, tolerance = 1e-9)
  expect_true(all(vapply(h[-1], function(e) is.null(e$storage$J), TRUE)))
  expect_equal(h[[2]]$storage$u, c(1.5180837643025114, 0.8553298855799086), tolerance = 1e-9)
  expect_equal(h[[6]]$storage$u, c(sqrt(2), 1), tolerance = 1e-12)
})

