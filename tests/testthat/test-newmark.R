# NewmarkBeta and GeneralizedAlpha (R/solvers-newmark.R) against
# OrdinaryDiffEqNewmark 2.4.0: the reference values come from
# tools/newmark-probe.jl (out of place functions, state c(v, u)).

newmark_test_probs <- list(
  L1 = function() DynamicalODEProblem(function(v, u, p, t) -4 * u, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)),
  P1 = function() DynamicalODEProblem(function(v, u, p, t) -u - u * u * u, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)),
  P5 = function() DynamicalODEProblem(function(v, u, p, t) -u + t, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)))

newmark_refs <- list(
  list(alg = "NewmarkBeta()", case = "L1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-2.0291263860455326, -0.2417805275718581, -0.18210306796824965, -0.21882765849354713),
       tq = 0.55, uq = c(-1.5506948276259125, -0.4928777847063105, 0.6792368419031105, -0.042064466779610135)),
  list(alg = "NewmarkBeta()", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2081059578085984, -0.33098387732371526, 0.5098918713184697, -0.14504424077852568),
       tq = 0.55, uq = c(-0.6635492331511565, -0.3614469782262172, 0.9500423896218934, 0.013345105716840201)),
  list(alg = "NewmarkBeta()", case = "P5", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.06947120731902955, 0.170543849648136, 1.143441621965683, 0.03782289363794345),
       tq = 0.55, uq = c(0.07743529359597555, -0.18693820857372598, 1.1489773754793364, 0.04892799187740188)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "L1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.9939360826708608, -0.23441669729431777, -0.16948245206079507, -0.21274378377426323),
       tq = 0.55, uq = c(-1.5381818938699545, -0.4849299809948878, 0.6807303301370811, -0.03984264421441626)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.194760019627542, -0.32801174870883765, 0.5132508894148597, -0.14342536710691875),
       tq = 0.55, uq = c(-0.6607243109337165, -0.3596032729509321, 0.9495775074277552, 0.013859394035646028)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "P5", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.07043555494324834, 0.17172748210359898, 1.142384758266111, 0.03890791546307848),
       tq = 0.55, uq = c(0.07606044591784014, -0.185521380882261, 1.148481005574743, 0.04935819594007701)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "L1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-2.001125907840663, -0.23366928222009017, -0.18078557712552026, -0.21661988488904246),
       tq = 0.55, uq = c(-1.5427086834013175, -0.4869665760928246, 0.678028210719486, -0.042119453672924495)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.1975916634390693, -0.329052254335963, 0.5100669482318778, -0.14460383619479458),
       tq = 0.55, uq = c(-0.6619710119603385, -0.360300834204774, 0.9493721273806611, 0.01334319407414489)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "P5", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.09384726383143266, 0.147715147796202, 1.1309444927492849, 0.026169706538774953),
       tq = 0.55, uq = c(0.06198207540653892, -0.20068541794981082, 1.1452553734388227, 0.045268325586893934)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "L1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-2.015062409158462, -0.2388519602601632, -0.17695014785483712, -0.21635991800041388),
       tq = 0.55, uq = c(-1.5456885621980092, -0.4896997016447683, 0.6798634014822913, -0.04115415574268554)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2022555051731292, -0.3297941011757831, 0.5114175669728889, -0.14438245284895512),
       tq = 0.55, uq = c(-0.6623212587719987, -0.3607000694146051, 0.9498804158554205, 0.013558577355484854)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "P5", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.07489172298407268, 0.16598224149939964, 1.1402632549775575, 0.03550812113604603),
       tq = 0.55, uq = c(0.07375537646859354, -0.18950140507161287, 1.1478914433296203, 0.04821941202296933)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "L1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.990842674618328, -0.23295448146184397, -0.17157173193864395, -0.21317612543655098),
       tq = 0.55, uq = c(-1.537793779498605, -0.4843804330874942, 0.6799814102143429, -0.04031564606423612)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.1931456238344553, -0.32790810842088547, 0.5127766944866977, -0.14358895888424777),
       tq = 0.55, uq = c(-0.6606270979174573, -0.3595494377083224, 0.9494452773001087, 0.013756498508319271)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "P5", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.0831204563649578, 0.15917194075499635, 1.1357014987300043, 0.03209686762445264),
       tq = 0.55, uq = c(0.06819445297667587, -0.19333018795173926, 1.1464356755346967, 0.047145810534630087)),
  list(alg = "NewmarkBeta()", case = "L1", nf = 27, nf2 = 0, naccept = 22, nreject = 2, nt = 23,
       t2 = 0.0015284805965253485, u = c(-2.027481915144206, -0.23982626224928177, -0.18662601892767938, -0.21936463020251168),
       tq = 0.55, uq = c(-1.5542078963060926, -0.4926516188835249, 0.6772089047850118, -0.04270526362639512)),
  list(alg = "NewmarkBeta()", case = "P1", nf = 19, nf2 = 0, naccept = 16, nreject = 0, nt = 17,
       t2 = 0.002058009239399836, u = c(-1.2086018301833605, -0.3309502869054824, 0.5091834493354043, -0.1450810506088163),
       tq = 0.55, uq = c(-0.6646581620305582, -0.3614349593384299, 0.9496610212486507, 0.013293725599288486)),
  list(alg = "NewmarkBeta()", case = "P1", tol = 1.0e-6, nf = 123, nf2 = 0, naccept = 116, nreject = 4, nt = 117,
       t2 = 0.00010855061738442037, u = c(-1.2083521247863098, -0.33085262829731166, 0.5068631919415164, -0.1453058085502234)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "L1", nf = 42, nf2 = 0, naccept = 30, nreject = 9, nt = 31,
       t2 = 2.3362529339544852e-6, u = c(-2.0128425858976122, -0.23637342179424367, -0.18265433504344558, -0.2171000009370736),
       tq = 0.55, uq = c(-1.5499947849759503, -0.4895374758101681, 0.6774247817755646, -0.0419667733357531)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "P1", nf = 27, nf2 = 0, naccept = 22, nreject = 2, nt = 23,
       t2 = 4.23540202945509e-6, u = c(-1.200976411408695, -0.32927710488127354, 0.5096937662693607, -0.1443642066229708),
       tq = 0.55, uq = c(-0.6641476029368545, -0.3604497842936053, 0.9490301091746575, 0.013500911385427661)),
  list(alg = "NewmarkBeta(0.3, 0.6)", case = "P1", tol = 1.0e-6, nf = 324, nf2 = 0, naccept = 312, nreject = 9, nt = 313,
       t2 = 1.1783236534538827e-8, u = c(-1.207944762895436, -0.3307514139100102, 0.5069462184830403, -0.14524987425214048)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "L1", nf = 35, nf2 = 0, naccept = 25, nreject = 7, nt = 26,
       t2 = 2.3362529339544852e-6, u = c(-2.014699933171807, -0.23640560469877822, -0.1848962135755979, -0.2180245301120337),
       tq = 0.55, uq = c(-1.5505130731837062, -0.49003717205625535, 0.6770558951230594, -0.04246061933393357)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "P1", nf = 25, nf2 = 0, naccept = 18, nreject = 4, nt = 19,
       t2 = 4.23540202945509e-6, u = c(-1.2011922196072666, -0.32951082666276493, 0.5093627442018913, -0.1447279552475454),
       tq = 0.55, uq = c(-0.6633903774899634, -0.3605489011358985, 0.9494356677408927, 0.013323932026035757)),
  list(alg = "GeneralizedAlpha(rho_inf = 0.5)", case = "P1", tol = 1.0e-6, nf = 264, nf2 = 0, naccept = 241, nreject = 20, nt = 242,
       t2 = 1.1783236534538827e-8, u = c(-1.2080507033398158, -0.33077819915032664, 0.5069169471632697, -0.1452664955620319)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "L1", nf = 33, nf2 = 0, naccept = 24, nreject = 6, nt = 25,
       t2 = 2.3362529339544852e-6, u = c(-2.019768857964906, -0.2382063622061065, -0.18372359853948472, -0.2179949425898787),
       tq = 0.55, uq = c(-1.5515799193493647, -0.49099357173808156, 0.6776408083920702, -0.04217721138056477)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "P1", nf = 22, nf2 = 0, naccept = 16, nreject = 3, nt = 17,
       t2 = 4.23540202945509e-6, u = c(-1.2031574157946099, -0.3299482984858097, 0.5109277832032965, -0.14450104604319372),
       tq = 0.55, uq = c(-0.66249142287234, -0.36076713162199636, 0.9499479775405615, 0.013528351825374466)),
  list(alg = "GeneralizedAlpha(alpha_hht = -0.1)", case = "P1", tol = 1.0e-6, nf = 230, nf2 = 0, naccept = 204, nreject = 23, nt = 205,
       t2 = 1.1783236534538827e-8, u = c(-1.208094590018169, -0.33078944294265694, 0.506912434722382, -0.14527170614317023)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "L1", nf = 42, nf2 = 0, naccept = 30, nreject = 9, nt = 31,
       t2 = 2.3362529339544852e-6, u = c(-2.0123624548562487, -0.23611711988556153, -0.1829376720887888, -0.2171553769743366),
       tq = 0.55, uq = c(-1.5499029884843247, -0.4894362283198922, 0.6773159021674003, -0.04203967847302518)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "P1", nf = 27, nf2 = 0, naccept = 22, nreject = 2, nt = 23,
       t2 = 4.23540202945509e-6, u = c(-1.2004918592229419, -0.3292488273713366, 0.5095436086122791, -0.14441553283810038),
       tq = 0.55, uq = c(-0.6641161202199058, -0.36043329830424237, 0.9489909644969631, 0.013469345602138834)),
  list(alg = "GeneralizedAlpha(-0.2, 0.2, 0.36, 0.9)", case = "P1", tol = 1.0e-6, nf = 324, nf2 = 0, naccept = 312, nreject = 9, nt = 313,
       t2 = 1.1783236534538827e-8, u = c(-1.2079431874756972, -0.33075126641621655, 0.506945787824246, -0.14525006232262253)),
  list(alg = "NewmarkBeta(1/6, 1/2)", case = "P1", dt = 0.1, nf = 11, nf2 = 0, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2074692866690102, -0.3308453874411636, 0.5087715757904385, -0.14529359421133353))
)

newmark_check <- function(r) {
  alg <- eval(parse(text = r$alg))
  prob <- newmark_test_probs[[r$case]]()
  sol <- if (!is.null(r$dt)) solve(prob, alg, dt = r$dt, adaptive = FALSE)
    else if (!is.null(r$tol)) solve(prob, alg, reltol = r$tol, abstol = r$tol) else solve(prob, alg)
  lab <- paste(r$alg, r$case, if (!is.null(r$dt)) "fixed" else "adaptive", r$tol %||% "")
  st <- sol$stats
  expect_identical(c(st$nf, st$nf2, st$naccept, st$nreject, length(sol$t)),
                   as.integer(c(r$nf, r$nf2, r$naccept, r$nreject, r$nt)), label = lab)
  expect_equal(sol$t[2], r$t2, tolerance = 1e-12, label = lab)
  expect_equal(unname(sol$u[[length(sol$t)]]), r$u, tolerance = 1e-12, label = lab)
  if (!is.null(r$tq)) expect_equal(unname(sol(r$tq)), r$uq, tolerance = 1e-12, label = lab)
}

test_that("NewmarkBeta and GeneralizedAlpha match OrdinaryDiffEqNewmark", {
  for (r in newmark_refs) newmark_check(r)
})

test_that("the parameterizations and orders match Julia", {
  # tools/newmark-probe.jl: alg_order and the fields of the constructors
  expect_equal(NewmarkBeta()$order, 2L)
  expect_equal(NewmarkBeta(0.3, 0.6)$order, 1L)
  ga <- GeneralizedAlpha(rho_inf = 0.5)
  expect_equal(unname(ga$newmark), c(0.0, 0.3333333333333333, 0.17361111111111108, 0.8333333333333333))
  expect_equal(ga$order, 1L)
  expect_equal(unname(GeneralizedAlpha(alpha_hht = -0.1)$newmark), c(0.0, 0.1, 0.30250000000000005, 0.6))
  expect_equal(GeneralizedAlpha(rho_inf = 1)$order, 2L)
  expect_error(NewmarkBeta(0.6), "Beta outside admissible range")
  expect_error(NewmarkBeta(0.25, 1.5), "Gamma outside admissible range")
  expect_error(GeneralizedAlpha(), "Must provide either rho_inf")
  expect_error(GeneralizedAlpha(rho_inf = 2), "rho_inf must be in")
  expect_error(GeneralizedAlpha(alpha_hht = 0.5), "HHT alpha must be in")
  expect_error(GeneralizedAlpha(0.3, 0.2, 0.36, 0.9), "requires alpha_m <= alpha_f")
  expect_error(GeneralizedAlpha(0, 0.2, 0.1, 0.9), "requires beta >=")
})

test_that("the Zienkiewicz-Xie estimate refuses blind parameters", {
  # Julia throws an ArgumentError for NewmarkBeta(1/6, 1/2) with adaptive steps
  expect_error(solve(newmark_test_probs$P1(), NewmarkBeta(1 / 6, 1 / 2)), "Zienkiewicz-Xie")
})

test_that("saveat interpolates with the derivatives at both ends of the step", {
  # tools/newmark-probe.jl: saveat = c(0.15, 0.55) with dt = 0.1
  sol <- solve(newmark_test_probs$P1(), NewmarkBeta(), dt = 0.1, adaptive = FALSE, saveat = c(0.15, 0.55))
  expect_equal(sol$t, c(0.15, 0.55))
  expect_equal(unname(sol$u[[1]]), c(0.18271001546967117, -0.3273502648031634, 1.0514315840887474, 0.15289514175692762),
               tolerance = 1e-12)
  expect_equal(unname(sol$u[[2]]), c(-0.6635492331511565, -0.3614469782262172, 0.9500423896218934, 0.013345105716840205),
               tolerance = 1e-12)
})

test_that("a start at rest does not freeze the step size", {
  # tools/newmark-probe.jl, "at rest": the first steps leave the acceleration
  # unchanged to rounding and the estimate is 0 (taken as 1); as in
  # NonlinearSolve, the nonlinear solve takes one Newton step and keeps it
  # when it lowers the residual, so dt grows once the correction moves the
  # acceleration (steps 2 to 11 hold dt, step 12 grows it).
  prob <- SecondOrderODEProblem(function(v, u, p, t) -u, 0, 1, c(0, 1))
  sol <- solve(prob, NewmarkBeta(), reltol = 1e-6)
  expect_identical(c(sol$stats$nf, sol$stats$naccept, sol$stats$nreject), c(69L, 63L, 3L))
  expect_equal(sol$t[2:3], c(0.00011892071150027212, 0.011491469800730342), tolerance = 1e-12)
  sol <- solve(prob, GeneralizedAlpha(rho_inf = 0.8), reltol = 1e-6, maxiters = 500)
  expect_identical(c(sol$stats$nf, sol$stats$naccept, sol$stats$nreject), c(103L, 83L, 17L))
  expect_equal(sol$t[c(2, 3, 12, 13)],
               c(1.4142135623730952e-8, 1.6970562748477143e-8, 6.788225099390861e-8, 3.5072496346852787e-7),
               tolerance = 1e-12)
})

test_that("Newmark methods need a partitioned problem", {
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), NewmarkBeta(), dt = 0.1),
               "designed for partitioned ODE")
})
