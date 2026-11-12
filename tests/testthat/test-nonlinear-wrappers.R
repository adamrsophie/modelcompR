# Wrapped solvers. Reference values printed by
# tools/nonlinearsolve-wrappers-probe.jl (NonlinearSolve 4.32.0 with NLsolve,
# LeastSquaresOptim, FastLevenbergMarquardt, NLSolvers, SpeedMapping,
# FixedPointAcceleration and SIAMFANLEquations) and, for CMINPACK, by
# tools/nonlinearsolve-cminpack-reference.cpp (the MINPACK port of Eigen: the
# cminpack_jll binary Julia loads computes every norm as garbage).

.wx_TT <- c(0, 0.5, 1, 1.5, 2, 2.5)
.wx_YY <- c(1, 1.6, 2.8, 4.4, 7.6, 12)
.wx_fun <- list(
  circ = list(function(u, p) c(u[1]^2 + u[2]^2 - 4, u[1] * u[2] - 1),
              function(u, p) rbind(c(2 * u[1], 2 * u[2]), c(u[2], u[1])), c(3, 0.1)),
  rosen = list(function(u, p) c(10 * (u[2] - u[1]^2), 1 - u[1]),
               function(u, p) rbind(c(-20 * u[1], 10), c(-1, 0)), c(-1.2, 1)),
  fixp = list(function(u, p) c(0.5 * cos(u[2]) - u[1] + 0.1, 0.5 * sin(u[1]) - u[2] + 0.2),
              function(u, p) rbind(c(-1, -0.5 * sin(u[2])), c(0.5 * cos(u[1]), -1)), c(0, 0)),
  burden = list(function(u, p) c(3 * u[1] - cos(u[2] * u[3]) - 0.5,
                                 u[1]^2 - 81 * (u[2] + 0.1)^2 + sin(u[3]) + 1.06,
                                 exp(-u[1] * u[2]) + 20 * u[3] + (10 * pi - 3) / 3),
                function(u, p) rbind(c(3, u[3] * sin(u[2] * u[3]), u[2] * sin(u[2] * u[3])),
                                     c(2 * u[1], -162 * (u[2] + 0.1), cos(u[3])),
                                     c(-u[2] * exp(-u[1] * u[2]), -u[1] * exp(-u[1] * u[2]), 20)),
                c(0.1, 0.1, -0.1)),
  expfit = list(function(u, p) u[1] * exp(u[2] * .wx_TT) - .wx_YY,
                function(u, p) cbind(exp(u[2] * .wx_TT), u[1] * .wx_TT * exp(u[2] * .wx_TT)), c(1, 0.5)))

.wx_prob <- function(pn) {
  if (pn %in% c("cubic", "cubicfar"))
    return(NonlinearProblem(function(u, p) u * u * u - 2 * u - 5, u0 = if (pn == "cubic") 2 else 10))
  base <- sub("0$", "", pn)
  fj <- .wx_fun[[base]]
  jac <- if (endsWith(pn, "0")) NULL else fj[[2]]
  if (base == "expfit")
    return(NonlinearLeastSquaresProblem(NonlinearFunction(fj[[1]], jac = jac, resid_prototype = numeric(6)), u0 = fj[[3]]))
  NonlinearProblem(NonlinearFunction(fj[[1]], jac = jac), u0 = fj[[3]])
}

.wx_alg <- function(pn, tag) switch(tag,
  CMINPACK = list(CMINPACK(), list(store_trace = startsWith(pn, "expfit"))),
  CMINPACKTrace = list(CMINPACK(), list(store_trace = TRUE)),
  CMINPACKlm = list(CMINPACK(method = "lm"), list(store_trace = TRUE)),
  NLsolveTR = list(NLsolveJL(), list()),
  NLsolveTRNoScale = list(NLsolveJL(autoscale = FALSE), list()),
  NLsolveNewton = list(NLsolveJL(method = "newton"), list()),
  NLsolveNewtonBT = list(NLsolveJL(method = "newton", linesearch = "BackTracking"), list()),
  NLsolveAnderson = list(NLsolveJL(method = "anderson"), list()),
  NLsolveAnderson2 = list(NLsolveJL(method = "anderson", m = 2, beta = 0.5), list()),
  NLsolveBroyden = list(NLsolveJL(method = "broyden"), list()),
  NLSolvers = list(NLSolversJL(), list()),
  LSOlm = , LSOlmJ = list(LeastSquaresOptimJL("lm"), list()),
  LSOdogleg = , LSOdoglegJ = list(LeastSquaresOptimJL("dogleg"), list()),
  LSOlmChol = , LSOlmCholJ = list(LeastSquaresOptimJL("lm", linsolve = "cholesky"), list()),
  FastLM = list(FastLevenbergMarquardtJL(), list()),
  FastLMqr = list(FastLevenbergMarquardtJL("qr"), list()),
  FastLMnielsen = list(FastLevenbergMarquardtJL(factorupdate = "nielson"), list()),
  SIAMnewton = list(SIAMFANLEquationsJL(), list()),
  SIAMptc = list(SIAMFANLEquationsJL(method = "pseudotransient"), list()),
  SIAMnewtonGMRES = list(SIAMFANLEquationsJL(linsolve = "gmres"), list()),
  SIAMnewtonBiCGSTAB = list(SIAMFANLEquationsJL(linsolve = "bicgstab"), list()),
  SIAMptcGMRES = list(SIAMFANLEquationsJL(method = "pseudotransient", linsolve = "gmres", delta = 1e-6), list()),
  SIAMptcBiCGSTAB = list(SIAMFANLEquationsJL(method = "pseudotransient", linsolve = "bicgstab", delta = 1e-6), list()),
  SIAMaa0 = list(SIAMFANLEquationsJL(method = "anderson"), list()),
  SIAMaa2 = list(SIAMFANLEquationsJL(method = "anderson", m = 2), list()),
  SIAMsecant = list(SIAMFANLEquationsJL(method = "secant"), list()),
  SpeedMapping = list(SpeedMappingJL(), list()),
  FPAAnderson = list(FixedPointAccelerationJL(), list()),
  FPAMPE = list(FixedPointAccelerationJL(algorithm = "MPE"), list()),
  FPARRE = list(FixedPointAccelerationJL(algorithm = "RRE"), list()),
  FPAVEA = list(FixedPointAccelerationJL(algorithm = "VEA"), list()),
  FPASEA = list(FixedPointAccelerationJL(algorithm = "SEA"), list()),
  FPASimple = list(FixedPointAccelerationJL(algorithm = "Simple"), list()),
  FPAAitken = list(FixedPointAccelerationJL(algorithm = "Aitken"), list()),
  FPANewton = list(FixedPointAccelerationJL(algorithm = "Newton"), list()))

.wx_solve <- function(pn, tag) {
  a <- .wx_alg(pn, tag)
  do.call(solve, c(list(.wx_prob(pn), a[[1]]), a[[2]]))
}

# problem, algorithm, retcode, u, stats (nf, njacs, nfactors, nsolve, nsteps),
# check: "exact" (u to 1e-12 and stats), "loose" (u to 1e-8: finite
# differences in R against ForwardDiff, or the stopping test of the nonzero
# residual fit, which reacts to rounding) or "diverge" (retcode only)
.wx_ref <- list(
  list("circ", "NLsolveTR", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, 6, 6, 5), "exact"),
  list("circ", "NLsolveTRNoScale", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, 6, 6, 5), "exact"),
  list("circ", "NLsolveNewton", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, 6, 6, 5), "exact"),
  list("circ", "NLsolveNewtonBT", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 5, 5, 5, 5), "exact"),
  list("circ", "NLsolveAnderson", "Success", c(1.9318516525781366, 0.5176380902050415), c(13, 0, 0, 0, 13), "exact"),
  list("circ", "NLsolveAnderson2", "Success", c(1.9318516525781217, 0.5176380902050564), c(12, 0, 0, 0, 12), "exact"),
  list("circ", "NLsolveBroyden", "Success", c(-1.931851652578137, -0.5176380902050413), c(93, 0, 0, 0, 16), "exact"),
  list("circ", "NLSolvers", "Success", c(1.9318516525781388, 0.5176380902050394), c(-1, -1, -1, -1, 5), "exact"),
  list("circ", "LSOlmJ", "Success", c(1.9318516525781366, 0.5176380902050415), c(9, 8, -1, -1, 8), "exact"),
  list("circ", "LSOdoglegJ", "Success", c(1.9318516525781366, 0.5176380902050415), c(7, 6, -1, -1, 6), "exact"),
  list("circ", "LSOlmCholJ", "Success", c(1.9318516525781366, 0.5176380902050415), c(9, 8, -1, -1, 8), "exact"),
  list("circ", "FastLM", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, -1, -1, 5), "exact"),
  list("circ", "FastLMqr", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, -1, -1, 5), "exact"),
  list("circ", "FastLMnielsen", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 6, -1, -1, 5), "exact"),
  list("circ", "SIAMnewton", "Success", c(1.9318516525781388, 0.5176380902050394), c(6, 5, 0, 0, 0), "exact"),
  list("circ", "SIAMptc", "Success", c(1.93185165257829, 0.5176380902047196), NULL, "exact"),
  list("circ", "SIAMnewtonGMRES", "Success", c(1.931851652578341, 0.5176380902048321), c(6, 9, 0, 0, 0), "exact"),
  list("circ", "SIAMnewtonBiCGSTAB", "Success", c(1.9318516525781366, 0.5176380902050415), c(7, 20, 0, 0, 0), "exact"),
  list("circ", "SIAMptcGMRES", "MaxIters", c(2.994990207059891, 0.10070126421321252), NULL, "exact"),
  list("circ", "SIAMptcBiCGSTAB", "MaxIters", c(2.994990207415286, 0.1007012667521591), NULL, "exact"),
  list("circ", "SIAMaa0", "Unstable", c(4801.94154201, -446.4098206), NULL, "exact"),
  list("circ", "SIAMaa2", "Success", c(1.9318516525781366, 0.5176380902050415), NULL, "exact"),
  list("circ", "SpeedMapping", "Failure", c(110.90018063783646, 14.367664963006973), c(1001, 0, 0, 0, 1001), "exact"),
  list("circ", "FPAAnderson", "Success", c(1.9318516525781377, 0.5176380902050415), c(14, 0, 0, 0, 14), "exact"),
  list("circ", "FPAMPE", "Failure", c(3.0, 0.1), c(13, 0, 0, 0, 13), "exact"),
  list("circ", "FPARRE", "Failure", c(3.0, 0.1), c(1000, 0, 0, 0, 1000), "exact"),
  list("circ", "FPAVEA", "Failure", c(3.0, 0.1), c(29, 0, 0, 0, 29), "exact"),
  list("circ", "FPASEA", "Failure", c(3.0, 0.1), c(24, 0, 0, 0, 24), "exact"),
  list("circ", "FPASimple", "Failure", c(3.0, 0.1), c(9, 0, 0, 0, 9), "exact"),
  list("circ", "FPAAitken", "Failure", c(3.0, 0.1), c(23, 0, 0, 0, 23), "exact"),
  list("circ", "FPANewton", "Failure", c(1.9318516525779426, 0.5176380902051382), c(80, 0, 0, 0, 80), "exact"),
  list("rosen", "NLsolveTR", "Success", c(1.0, 1.0), c(30, 17, 17, 17, 29), "exact"),
  list("rosen", "NLsolveTRNoScale", "Success", c(1.0, 1.0), c(24, 13, 13, 13, 23), "exact"),
  list("rosen", "NLsolveNewton", "Success", c(0.9999999999999998, 1.0), c(3, 3, 3, 3, 2), "exact"),
  list("rosen", "NLsolveNewtonBT", "Success", c(1.0, 1.0), c(27, 14, 14, 14, 14), "exact"),
  list("rosen", "NLsolveAnderson2", "Success", c(1.0, 1.0000000000000002), c(28, 0, 0, 0, 28), "exact"),
  list("rosen", "NLsolveBroyden", "Success", c(1.0, 1.0000000000000002), c(847, 0, 0, 0, 52), "exact"),
  list("rosen", "NLSolvers", "Success", c(0.9999999999999998, 1.0), c(-1, -1, -1, -1, 2), "exact"),
  list("rosen", "LSOlmJ", "Success", c(1.0, 1.0), c(38, 29, -1, -1, 37), "exact"),
  list("rosen", "LSOdoglegJ", "Success", c(1.0, 1.0), c(22, 13, -1, -1, 21), "exact"),
  list("rosen", "LSOlmCholJ", "Success", c(1.0, 1.0), c(38, 29, -1, -1, 37), "exact"),
  list("rosen", "FastLM", "Success", c(0.9999999999999917, 0.9999999999999833), c(30, 15, -1, -1, 29), "exact"),
  list("rosen", "FastLMqr", "Success", c(0.9999999999999917, 0.9999999999999833), c(30, 15, -1, -1, 29), "exact"),
  list("rosen", "FastLMnielsen", "Success", c(1.0, 1.0), c(37, 20, -1, -1, 36), "exact"),
  list("rosen", "SIAMnewton", "Success", c(1.0, 1.0), c(94, 26, 0, 0, 67), "exact"),
  list("rosen", "SIAMptc", "MaxIters", c(3.1570642333819015, -0.326625485366076), NULL, "exact"),
  list("rosen", "SIAMnewtonGMRES", "Success", c(0.9999999999998858, 0.9999999999997708), c(100, 60, 0, 0, 67), "exact"),
  list("rosen", "SIAMnewtonBiCGSTAB", "Success", c(1.0, 1.0), c(97, 118, 0, 0, 64), "exact"),
  list("rosen", "SIAMptcGMRES", "MaxIters", c(-1.1956076191224736, 0.9977850191381877), NULL, "exact"),
  list("rosen", "SIAMptcBiCGSTAB", "MaxIters", c(-1.1956076342520667, 0.9977849891355544), NULL, "exact"),
  list("rosen", "SIAMaa0", "Unstable", c(-287.2, 9.8), NULL, "exact"),
  list("rosen", "SIAMaa2", "Success", c(1.000000000000006, 1.0000000000000018), NULL, "exact"),
  list("rosen", "SpeedMapping", "Failure", c(1.0087463000469965, 1.0240654253608712), c(1001, 0, 0, 0, 1001), "exact"),
  list("rosen", "FPAAnderson", "Success", c(1.0, 1.0), c(9, 0, 0, 0, 9), "exact"),
  list("rosen", "FPAMPE", "Failure", c(-1.2, 1.0), c(10, 0, 0, 0, 10), "exact"),
  list("rosen", "FPARRE", "Failure", c(-1.2, 1.0), c(1000, 0, 0, 0, 1000), "exact"),
  list("rosen", "FPAVEA", "Failure", c(-1.2, 1.0), c(19, 0, 0, 0, 19), "exact"),
  list("rosen", "FPASEA", "Failure", c(-1.2, 1.0), c(19, 0, 0, 0, 19), "exact"),
  list("rosen", "FPASimple", "Failure", c(-1.2, 1.0), c(8, 0, 0, 0, 8), "exact"),
  list("rosen", "FPAAitken", "Failure", c(-1.2, 1.0), c(20, 0, 0, 0, 20), "exact"),
  list("rosen", "FPANewton", "Failure", c(-1.2, 1.0), c(1000, 0, 0, 0, 1000), "exact"),
  list("fixp", "NLsolveTR", "Success", c(0.5479244644090591, 0.4604583296926159), c(5, 5, 5, 5, 4), "exact"),
  list("fixp", "NLsolveTRNoScale", "Success", c(0.5479244644090591, 0.4604583296926159), c(5, 5, 5, 5, 4), "exact"),
  list("fixp", "NLsolveNewton", "Success", c(0.5479244644090591, 0.4604583296926159), c(5, 5, 5, 5, 4), "exact"),
  list("fixp", "NLsolveNewtonBT", "Success", c(0.5479244644090591, 0.4604583296926159), c(5, 4, 4, 4, 4), "exact"),
  list("fixp", "NLsolveAnderson", "Success", c(0.5479244644090588, 0.4604583296926071), c(8, 0, 0, 0, 8), "exact"),
  list("fixp", "NLsolveAnderson2", "Success", c(0.5479244644090617, 0.460458329692568), c(8, 0, 0, 0, 8), "exact"),
  list("fixp", "NLsolveBroyden", "Success", c(0.5479244644090593, 0.46045832969261524), c(93, 0, 0, 0, 15), "exact"),
  list("fixp", "NLSolvers", "Success", c(0.5479244644090591, 0.4604583296926159), c(-1, -1, -1, -1, 4), "exact"),
  list("fixp", "LSOlmJ", "Success", c(0.547924464409059, 0.4604583296926156), c(8, 7, -1, -1, 7), "exact"),
  list("fixp", "LSOdoglegJ", "Success", c(0.547924464409059, 0.4604583296926156), c(6, 5, -1, -1, 5), "exact"),
  list("fixp", "LSOlmCholJ", "Success", c(0.547924464409059, 0.4604583296926156), c(8, 7, -1, -1, 7), "exact"),
  list("fixp", "FastLM", "Success", c(0.5479244644090591, 0.46045832969261585), c(5, 5, -1, -1, 4), "exact"),
  list("fixp", "FastLMqr", "Success", c(0.5479244644090591, 0.46045832969261585), c(5, 5, -1, -1, 4), "exact"),
  list("fixp", "FastLMnielsen", "Success", c(0.5479244644090591, 0.46045832969261585), c(5, 5, -1, -1, 4), "exact"),
  list("fixp", "SIAMnewton", "Success", c(0.5479244644090591, 0.4604583296926159), c(5, 4, 0, 0, 0), "exact"),
  list("fixp", "SIAMptc", "MaxIters", c(-0.6240462154941342, -0.08686544254893445), NULL, "exact"),
  list("fixp", "SIAMnewtonGMRES", "Success", c(0.5479244644090587, 0.4604583296926152), c(5, 8, 0, 0, 0), "exact"),
  list("fixp", "SIAMnewtonBiCGSTAB", "Success", c(0.5479244644089127, 0.4604583296923063), c(8, 16, 0, 0, 0), "exact"),
  list("fixp", "SIAMptcGMRES", "MaxIters", c(-0.0006000454360974666, -0.00019986532376708793), NULL, "exact"),
  list("fixp", "SIAMptcBiCGSTAB", "MaxIters", c(-0.0006000455260630193, -0.00019986505366850004), NULL, "exact"),
  list("fixp", "SIAMaa0", "Success", c(0.5479244644088233, 0.4604583296924911), NULL, "exact"),
  list("fixp", "SIAMaa2", "Success", c(0.5479244644090588, 0.4604583296926071), NULL, "exact"),
  list("fixp", "SpeedMapping", "Success", c(0.5479244644091313, 0.46045832969264), c(26, 0, 0, 0, 26), "exact"),
  list("fixp", "FPAAnderson", "Success", c(0.5479244644090608, 0.4604583296926155), c(8, 0, 0, 0, 8), "exact"),
  list("fixp", "FPAMPE", "Success", c(0.5479244644090864, 0.4604583296926283), c(16, 0, 0, 0, 16), "exact"),
  list("fixp", "FPARRE", "Success", c(0.5479244644090592, 0.4604583296926153), c(19, 0, 0, 0, 19), "exact"),
  list("fixp", "FPAVEA", "Success", c(0.5479244644090034, 0.4604583296926101), c(14, 0, 0, 0, 14), "exact"),
  list("fixp", "FPASEA", "Success", c(0.5479244644090591, 0.46045832969257733), c(15, 0, 0, 0, 15), "exact"),
  list("fixp", "FPASimple", "Success", c(0.5479244644090866, 0.460458329692515), c(25, 0, 0, 0, 25), "exact"),
  list("fixp", "FPAAitken", "Success", c(0.5479244644090902, 0.46045832969263206), c(38, 0, 0, 0, 38), "exact"),
  list("fixp", "FPANewton", "Success", c(0.5479244644091129, 0.4604583296926333), c(49, 0, 0, 0, 49), "exact"),
  list("burden", "NLsolveTR", "Success", c(0.5, -1.1667064928554629e-17, -0.5235987755982988), c(7, 7, 7, 7, 6), "exact"),
  list("burden", "NLsolveTRNoScale", "Success", c(0.5, 8.189882571253331e-18, -0.5235987755982989), c(7, 7, 7, 7, 6), "exact"),
  list("burden", "NLsolveNewton", "Success", c(0.5, 3.3552188928291758e-18, -0.5235987755982989), c(6, 6, 6, 6, 5), "exact"),
  list("burden", "NLsolveNewtonBT", "Success", c(0.5, 3.3552188928291758e-18, -0.5235987755982989), c(6, 5, 5, 5, 5), "exact"),
  list("burden", "NLsolveAnderson", "Success", c(0.49999999999999994, -6.345197237738916e-16, -0.5235987755982989), c(21, 0, 0, 0, 21), "exact"),
  list("burden", "NLsolveBroyden", "Success", c(0.5, -2.240779183650091e-16, -0.5235987755982989), c(108, 0, 0, 0, 18), "exact"),
  list("burden", "NLSolvers", "Success", c(0.5, 3.3552188928291758e-18, -0.5235987755982989), c(-1, -1, -1, -1, 5), "exact"),
  list("burden", "LSOlmJ", "Success", c(0.5, -6.359067014693155e-18, -0.5235987755982989), c(9, 8, -1, -1, 8), "exact"),
  list("burden", "LSOdoglegJ", "Success", c(0.5, 9.169324755464393e-18, -0.5235987755982988), c(8, 7, -1, -1, 7), "exact"),
  list("burden", "LSOlmCholJ", "Success", c(0.5, -8.988468421023769e-18, -0.5235987755982989), c(9, 8, -1, -1, 8), "exact"),
  list("burden", "FastLM", "Success", c(0.5, 1.2189383964598782e-17, -0.5235987755982989), c(6, 6, -1, -1, 5), "exact"),
  list("burden", "FastLMqr", "Success", c(0.5, -2.2040485915416845e-18, -0.5235987755982989), c(6, 6, -1, -1, 5), "exact"),
  list("burden", "FastLMnielsen", "Success", c(0.5, 1.2189383964598782e-17, -0.5235987755982989), c(6, 6, -1, -1, 5), "exact"),
  list("burden", "SIAMnewton", "Success", c(0.5, 3.3552188928291758e-18, -0.5235987755982989), c(6, 5, 0, 0, 0), "exact"),
  list("burden", "SIAMptc", "MaxIters", c(0.1661089873399955, 8.88411818014502, -0.40450270304794433), NULL, "exact"),
  list("burden", "SIAMnewtonGMRES", "Success", c(0.49999999999992667, -4.5709144953415136e-15, -0.5235987755982984), c(11, 19, 0, 0, 0), "exact"),
  list("burden", "SIAMnewtonBiCGSTAB", "Success", c(0.49999999999998385, -1.3213747508803257e-15, -0.5235987755982989), c(7, 28, 0, 0, 0), "exact"),
  list("burden", "SIAMptcGMRES", "MaxIters", c(0.10120772976763456, 0.10233006113773, -0.10844430217513744), NULL, "exact"),
  list("burden", "SIAMptcBiCGSTAB", "MaxIters", c(0.10120774522885219, 0.10233018483965574, -0.1084442658170176), NULL, "exact"),
  list("burden", "SIAMaa0", "Unstable", c(-5.661182169446718, -346.0472857468996, 185.16643957323782), NULL, "exact"),
  list("burden", "SIAMaa2", "Unstable", c(-22.00536017792209, 0.641744826879492, -1.9419426284000758), NULL, "exact"),
  list("burden", "SpeedMapping", "Failure", c(-3.636162754105527, 1.8189894035458565e-12, 48.406873857696155), NULL, "diverge"),
  list("burden", "FPAAnderson", "Success", c(0.5, -1.3335706582591802e-14, -0.523598775598297), c(19, 0, 0, 0, 19), "exact"),
  list("burden", "FPAMPE", "Failure", c(0.1, 0.1, -0.1), c(6, 0, 0, 0, 6), "exact"),
  list("burden", "FPARRE", "Failure", c(0.1, 0.1, -0.1), c(1000, 0, 0, 0, 1000), "exact"),
  list("burden", "FPAVEA", "Failure", c(0.1, 0.1, -0.1), c(22, 0, 0, 0, 22), "exact"),
  list("burden", "FPASEA", "Failure", c(0.1, 0.1, -0.1), c(7, 0, 0, 0, 7), "exact"),
  list("burden", "FPASimple", "Failure", c(0.1, 0.1, -0.1), c(8, 0, 0, 0, 8), "exact"),
  list("burden", "FPAAitken", "Failure", c(0.1, 0.1, -0.1), c(5, 0, 0, 0, 5), "exact"),
  list("burden", "FPANewton", "Success", c(0.4981446845894977, -0.19960589554379177, -0.5288259775733912), c(358, 0, 0, 0, 358), "exact"),
  list("expfit", "LSOlmJ", "Success", c(1.0341927082552274, 0.9831024240402416), NULL, "loose"),
  list("expfit", "LSOdoglegJ", "Success", c(1.0341927083210898, 0.9831024240114542), c(9, 7, -1, -1, 8), "exact"),
  list("expfit", "FastLM", "Success", c(1.0341927082270266, 0.9831024240526219), c(23, 10, -1, -1, 22), "exact"),
  list("expfit", "FastLMqr", "Success", c(1.0341927082133302, 0.9831024240586161), c(23, 8, -1, -1, 22), "exact"),
  list("expfit", "FastLMnielsen", "Success", c(1.0341927080516755, 0.9831024241293667), c(16, 9, -1, -1, 15), "exact"),
  list("expfit0", "LSOlm", "Success", c(1.0341927093675318, 0.9831024235724787), NULL, "loose"),
  list("expfit0", "LSOdogleg", "Success", c(1.034192708256635, 0.9831024240428937), NULL, "loose"),
  list("expfit0", "LSOlmChol", "Success", c(1.0341927084654046, 0.9831024239652082), NULL, "loose"),
  list("expfit0", "LSOlmJ", "Success", c(1.0341927082552274, 0.9831024240402415), NULL, "loose"),
  list("expfit0", "LSOdoglegJ", "Success", c(1.0341927083210898, 0.9831024240114542), NULL, "loose"),
  list("expfit0", "FastLM", "Success", c(1.0341927082270261, 0.9831024240526219), NULL, "loose"),
  list("expfit0", "FastLMqr", "Success", c(1.0341927082270266, 0.9831024240526219), NULL, "loose"),
  list("expfit0", "FastLMnielsen", "Success", c(1.0341927080516755, 0.9831024241293667), NULL, "loose"),
  list("cubic", "SIAMsecant", "Success", c(2.0945514815423265), c(7, 0, 0, 0, 0), "exact"),
  list("cubicfar", "SIAMsecant", "Success", c(2.0945514815423265), c(14, 0, 0, 0, 0), "exact"),
  list("circ", "CMINPACK", "Success", c(1.9318516525781366, 0.51763809020504148), c(0, 0, 0, 0, -1), "exact"),
  list("circ", "CMINPACKTrace", "Success", c(1.9318516525781366, 0.51763809020504148), c(11, 1, 1, 1, -1), "exact"),
  list("circ", "CMINPACKlm", "Success", c(1.9318516525781366, 0.51763809020504148), c(7, 6, 6, 6, -1), "exact"),
  list("circ0", "CMINPACK", "Success", c(1.9318516525781364, 0.51763809020504148), c(0, 0, 0, 0, -1), "exact"),
  list("circ0", "CMINPACKTrace", "Success", c(1.9318516525781364, 0.51763809020504148), c(11, 1, 1, 1, -1), "exact"),
  list("circ0", "CMINPACKlm", "Success", c(1.9318516525781364, 0.51763809020504148), c(7, 6, 6, 6, -1), "exact"),
  list("rosen", "CMINPACK", "Success", c(1, 1), c(0, 0, 0, 0, -1), "exact"),
  list("rosen", "CMINPACKTrace", "Success", c(1, 1), c(16, 3, 3, 3, -1), "exact"),
  list("rosen", "CMINPACKlm", "Success", c(1, 1), c(21, 16, 16, 16, -1), "exact"),
  list("rosen0", "CMINPACK", "Success", c(1, 1), c(0, 0, 0, 0, -1), "exact"),
  list("rosen0", "CMINPACKTrace", "Success", c(1, 1), c(16, 3, 3, 3, -1), "exact"),
  list("rosen0", "CMINPACKlm", "Success", c(1, 1), c(22, 17, 17, 17, -1), "exact"),
  list("fixp", "CMINPACK", "Success", c(0.54792446440905895, 0.46045832969261558), c(0, 0, 0, 0, -1), "exact"),
  list("fixp", "CMINPACKTrace", "Success", c(0.54792446440905895, 0.46045832969261558), c(9, 1, 1, 1, -1), "exact"),
  list("fixp", "CMINPACKlm", "Success", c(0.54792446440905895, 0.46045832969261558), c(6, 5, 5, 5, -1), "exact"),
  list("fixp0", "CMINPACK", "Success", c(0.54792446440905895, 0.46045832969261558), c(0, 0, 0, 0, -1), "exact"),
  list("fixp0", "CMINPACKTrace", "Success", c(0.54792446440905895, 0.46045832969261558), c(9, 1, 1, 1, -1), "exact"),
  list("fixp0", "CMINPACKlm", "Success", c(0.54792446440905906, 0.46045832969261558), c(6, 5, 5, 5, -1), "exact"),
  list("burden", "CMINPACK", "Success", c(0.5, 6.1817552610436517e-18, -0.52359877559829893), c(0, 0, 0, 0, -1), "exact"),
  list("burden", "CMINPACKTrace", "Success", c(0.5, 6.1817552610436517e-18, -0.52359877559829893), c(10, 1, 1, 1, -1), "exact"),
  list("burden", "CMINPACKlm", "Success", c(0.5, 1.1885082001623693e-18, -0.52359877559829893), c(7, 6, 6, 6, -1), "exact"),
  list("burden0", "CMINPACK", "Success", c(0.5, -3.4315066862089832e-19, -0.52359877559829893), c(0, 0, 0, 0, -1), "exact"),
  list("burden0", "CMINPACKTrace", "Success", c(0.5, -3.4315066862089832e-19, -0.52359877559829893), c(10, 1, 1, 1, -1), "exact"),
  list("burden0", "CMINPACKlm", "Failure", c(0.49999999981031212, -1.1744584226992683e-10, -0.52359877560493118), NULL, "loose"),
  list("expfit", "CMINPACK", "Failure", c(1.0341927080596549, 0.98310242412587356), c(10, 8, 8, 8, -1), "exact"),
  list("expfit0", "CMINPACK", "Failure", c(1.0341927088849527, 0.98310242374559598), NULL, "loose"))

test_that("wrapped solvers follow the reference runs", {
  skip_on_cran()
  for (rf in .wx_ref) {
    sol <- .wx_solve(rf[[1]], rf[[2]])
    lab <- paste(rf[[1]], rf[[2]])
    expect_identical(sol$retcode, rf[[3]], label = lab)
    if (rf[[6]] == "diverge") next
    tol <- if (rf[[6]] == "exact") 1e-12 else 1e-8
    expect_lte(max(abs(unname(sol$u) - rf[[4]]) / pmax(1, abs(rf[[4]]))), tol, label = lab)
    if (!is.null(rf[[5]])) {
      s <- sol$stats
      expect_equal(c(s$nf, s$njacs, s$nfactors, s$nsolve, s$nsteps), rf[[5]], label = lab)
    }
  }
})

test_that("wrapped solvers report non finite residuals as NLsolve does", {
  expect_error(.wx_solve("rosen", "NLsolveAnderson"), "non-finite number")
  expect_error(.wx_solve("burden", "NLsolveAnderson2"), "non-finite number")
})

test_that("the fused multiply add emulation rounds once", {
  expect_identical(.wr_fma(1 + 2^-30, 1 - 2^-30, -1), -2^-60)
  expect_identical(.wr_fma(3, 5, 7), 22)
  expect_identical(.wr_dot(c(1, 2, 3), c(4, 5, 6)), 32)
})

test_that("wrapper constructors validate their options", {
  expect_error(SpeedMappingJL(orders = c(3L, 3L)), "at least one 2")
  expect_error(SpeedMappingJL(orders = c(2L, 4L)), "2 or 3")
  expect_error(NLsolveJL(method = "secant"))
  expect_error(NLsolveJL(linesearch = "HagerZhang"))
  expect_error(LeastSquaresOptimJL("newton"))
  expect_error(LeastSquaresOptimJL(linsolve = "lu"))
  expect_error(FastLevenbergMarquardtJL("lu"))
  expect_identical(FastLevenbergMarquardtJL(factorupdate = "nielsen")$factorupdate, "nielson")
  expect_error(FastLevenbergMarquardtJL(factorupdate = "levenberg"))
  expect_error(FixedPointAccelerationJL(algorithm = "Broyden"))
  expect_error(FixedPointAccelerationJL(algorithm = "MPE", m = 3), "only valid for Anderson")
  expect_error(FixedPointAccelerationJL(algorithm = "SEA", extrapolation_period = 3), "multiples of 2")
  expect_error(SIAMFANLEquationsJL(method = "broyden"))
  expect_error(SIAMFANLEquationsJL(linsolve = "cg"))
  expect_error(CMINPACK(method = "dogleg"))
  expect_error(NLSolversJL(method = "TrustRegion"))
  expect_error(solve(.wx_prob("circ"), CMINPACK(method = "hybrd")), "unknown method hybrd")
  expect_error(solve(.wx_prob("circ"), CMINPACK(method = "lmdif")), "unknown method lmdif")
  expect_identical(solve(.wx_prob("circ0"), CMINPACK(method = "hybrd"))$retcode, "Success")
  expect_error(solve(.wx_prob("circ"), SIAMFANLEquationsJL(method = "secant")), "scalar equations")
})

test_that("PETScSNES runs Newton with backtracking", {
  sol <- solve(.wx_prob("circ"), PETScSNES())
  expect_identical(sol$retcode, "Success")
  expect_lte(max(abs(sol$resid)), 1e-10)
  expect_identical(sol$original$reason, "CONVERGED_FNORM_ABS")
  sol0 <- solve(.wx_prob("burden0"), PETScSNES(snes_linesearch_order = 2L))
  expect_identical(sol0$retcode, "Success")
  expect_lte(max(abs(sol0$resid)), 1e-10)
})

