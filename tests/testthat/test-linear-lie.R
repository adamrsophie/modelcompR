# Reference values from tools/linear-probe.jl (Julia, OrdinaryDiffEqLinear 2.3.1):
# final states and statistics (nf, naccept, nreject). For LieEuler, MagnusGauss4
# and CayleyEuler with dt = 0.1 upstream steps past tf (they cannot change dt)
# and appends tf; the values are those of its step before, and the counts
# leave that extra step out.
.LIN_REF <- list(
  `p1_MagnusMidpoint` = list(u = c(0.55510091148288254, -0.41021880944654437, -0.22476123164128201),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusLeapfrog` = list(u = c(0.55580688150521917, -0.41071747025155042, -0.22485606869283131),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_LieEuler` = list(u = c(0.55863979181579015, -0.40500880271017514, -0.22958584279663119),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusGauss4` = list(u = c(0.55512487874661676, -0.41018437311049499, -0.22464249397736577),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusNC6` = list(u = c(0.55508209669488451, -0.41023326544843813, -0.22475633465848951),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusGL6` = list(u = c(0.55508209669488451, -0.41023326544843813, -0.22475633465848951),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusGL8` = list(u = c(0.55508209669723263, -0.4102332654820417, -0.22475633466380252),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusNC8` = list(u = c(0.55508209669723263, -0.4102332654820417, -0.22475633466380252),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_MagnusGL4` = list(u = c(0.55508210073375908, -0.41023321457073358, -0.22475630899781501),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_RKMK2` = list(u = c(0.55926545590715515, -0.4055097340017757, -0.22973503986063731),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_RKMK4` = list(u = c(0.55930233693958908, -0.40550595463259059, -0.22971285674704545),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_LieRK4` = list(u = c(0.55930232350854248, -0.40550590954772336, -0.22971267123235764),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_CG2` = list(u = c(0.55924578202156294, -0.40552391924712566, -0.22973580333085725),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_CG3` = list(u = c(0.55574655451708965, -0.41074101295447135, -0.22491825947202509),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1_CG4a` = list(u = c(0.55574513728767982, -0.41073848827453135, -0.2249144097549359),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusMidpoint` = list(u = c(0.31396918897390114, 0.2443938402120534, 0.26226818965700399, 0.21612072058611273,
    0.28011864264817898, 0.097136616785022944),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusLeapfrog` = list(u = c(0.31336187601983834, 0.24432796803594783, 0.26239099662358195, 0.21613614552061869,
    0.27995565315406146, 0.097223509459610044),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_LieEuler` = list(u = c(0.31715534177241589, 0.24685198529553176, 0.26489998878754567, 0.21829232117597111,
    0.28293893317373769, 0.098110436227664494),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusGauss4` = list(u = c(0.31396918897390125, 0.24439384021205343, 0.26226818965700394, 0.21612072058611259,
    0.28011864264817898, 0.097136616785022944),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusNC6` = list(u = c(0.31396918897390114, 0.24439384021205349, 0.26226818965700388, 0.21612072058611256,
    0.28011864264817898, 0.097136616785022958),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusGL6` = list(u = c(0.31396918897390114, 0.24439384021205349, 0.26226818965700394, 0.21612072058611265,
    0.28011864264817904, 0.097136616785022958),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusGL8` = list(u = c(0.31396918897390114, 0.24439384021205349, 0.26226818965700388, 0.21612072058611256,
    0.28011864264817904, 0.097136616785022958),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusNC8` = list(u = c(0.31396918897390097, 0.24439384021205329, 0.26226818965700383, 0.21612072058611248,
    0.28011864264817898, 0.097136616785022917),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_MagnusGL4` = list(u = c(0.31396918897390119, 0.24439384021205349, 0.26226818965700388, 0.21612072058611259,
    0.28011864264817904, 0.097136616785022917),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_RKMK2` = list(u = c(0.31659579767658352, 0.2467958829698165, 0.26500437402544674, 0.21830774952681839,
    0.28279765177499777, 0.09818441040013616),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_RKMK4` = list(u = c(0.31657836437953529, 0.24679598553161664, 0.26500511681065109, 0.21830754622599616,
    0.28279759953211864, 0.098183919934412445),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2k_LieRK4` = list(u = c(0.31657788375119189, 0.24679592737549824, 0.2650048385254764, 0.21830798367582724,
    0.28279714672035688, 0.098184143524085402),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1dt_MagnusGL8` = list(u = c(0.55396120484848643, -0.40952011496905366, -0.22464670952502297),
    stats = c(5, 4, 0), nt = 5, tend = 1),
  `p1dt_CG3` = list(u = c(0.55577999619161467, -0.41079862503030601, -0.22504841718662391),
    stats = c(5, 4, 0), nt = 5, tend = 1),
  `p1dt_LieEuler` = list(u = c(0.56286401041591227, -0.39615100899354067, -0.23623273787698945),
    stats = c(5, 4, 0), nt = 5, tend = 1),
  `p1k_MagnusMidpoint` = list(u = c(0.59149236718379872, -0.4324938901941266, -0.23799255711423797),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1k_RKMK4` = list(u = c(0.59527393026695274, -0.42712715371386251, -0.24303985505239872),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p1k_MagnusGL6` = list(u = c(0.59147352934198782, -0.43250912961757348, -0.23798480438067071),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2kd_MagnusNC8` = list(u = c(0.31456836248904363, 0.2444475038981154, 0.26216478468260596, 0.21610700930997939,
    0.28025451238527943, 0.097067199822220326),
    stats = c(6, 5, 0), nt = 6, tend = 1),
  `p2kd_LieRK4` = list(u = c(0.31978734971062456, 0.24927802980361083, 0.26766447729919984, 0.22050163958640795,
    0.28564386343359949, 0.099168738829830044),
    stats = c(6, 5, 0), nt = 6, tend = 1),
  `p3k_MagnusMidpoint` = list(u = c(0.30236845273899221, 0.2906191242221029, 0.27540919090941102, 0.18530216313591097,
    0.22588751541313865, 0.1517387171439071),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p3k_LieEuler` = list(u = c(0.30542613078741226, 0.2935368135392431, 0.27817908323502871, 0.18716118413895216,
    0.22815588116714677, 0.15330925272010248),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p3k_MagnusGL6` = list(u = c(0.30236845273899232, 0.29061912422210257, 0.27540919090940991, 0.18530216313591036,
    0.22588751541313912, 0.1517387171439076),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p3k_MagnusGauss4` = list(u = c(0.30236845273899221, 0.29061912422210312, 0.27540919090941118, 0.18530216313591055,
    0.22588751541313876, 0.15173871714390724),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2_MagnusGL4` = list(u = c(0.31396921378699971, 0.24439384308980197, 0.26226777871227142, 0.21612116712958679,
    0.28011832554432059, 0.097136764841022016),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `p2_RKMK4` = list(u = c(0.31657839056830039, 0.2467959903227617, 0.26500470135147591, 0.21830799273580248,
    0.28279728254047104, 0.098184068306703801),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `a1_default` = list(u = c(-0.36463266317631138, 0.0082961505333659959, 0.027611466679208514),
    stats = c(16, 13, 0), nt = 14, tend = 3),
  `a1_tight` = list(u = c(-0.36463116175459992, 0.0083013862428154683, 0.027609306685045795),
    stats = c(87, 84, 0), nt = 85, tend = 3),
  `a2_default` = list(u = c(0.066445622099362292, 0.055136487224351077, 0.076357879930546935, 0.061556877044677161,
    0.09770179009634096, 0.002283452549691476),
    stats = c(15, 12, 0), nt = 13, tend = 2),
  `a1_fixed` = list(u = c(0.55573929979753234, -0.41073087065165753, -0.22489284281833233),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `le_nodt` = list(u = c(0.097669936109443523, 0.08163583424185944, 0.11451364254168242, 0.092054842994620439,
    0.1443384776447619, 0.0044922312802591624),
    stats = c(2, 1, 0), nt = 2, tend = 2),
  `le_off` = list(u = c(0.097669936109443536, 0.081635834241859426, 0.11451364254168242, 0.092054842994620439,
    0.1443384776447619, 0.0044922312802591312),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `le_simple` = list(u = c(0.097669936109443398, 0.081635834241859648, 0.11451364254168242, 0.092054842994620481,
    0.14433847764476193, 0.0044922312802589499),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `le_simple4` = list(u = c(0.097667488804277836, 0.081621643500488616, 0.11456592007681804, 0.092027906206475898,
    0.14436461800546072, 0.0044825330522629981),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `le_adapt` = list(u = c(0.097669936109443398, 0.081635834241859648, 0.11451364254168242, 0.092054842994620481,
    0.14433847764476193, 0.0044922312802589499),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `le_adapt3` = list(u = c(0.097669936752089126, 0.081635831782976565, 0.11451364731955577, 0.09205484032069608,
    0.14433847717738946, 0.0044922322680845551),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `le_adapt3n` = list(u = c(0.097589874774124313, 0.08156437331293713, 0.11489187216913711, 0.091941123955189025,
    0.14446744672077241, 0.0043316776598298322),
    stats = c(2, 1, 0), nt = 2, tend = 2),
  `le_sym` = list(u = c(0.10992538481117832, 0.1160698735448427, 0.11558613971626652, 0.054052718326913642,
    0.1111079215144067, 0.047523828938768317),
    stats = c(5, 4, 0), nt = 5, tend = 2),
  `cay` = list(u = c(1.0824144617719129, -0.14100935748475887, -0.15288963501196423, -0.050844382580941129,
    0.86367294348461854, -0.053158366674597932, 0.24444205978615779, 0.20857763609485505,
    0.94514925998646582),
    stats = c(11, 10, 0), nt = 11, tend = 1),
  `cays` = list(u = c(1.0847169561821894, -0.13800958269924024, -0.15265690638713511, -0.048575136481556834,
    0.86089081261771194, -0.055406648083618115, 0.24109557898982636, 0.21105631762074811,
    0.94359778951539897),
    stats = c(5, 4, 0), nt = 5, tend = 1)
)

lin_upd1 <- function(A, u, p, t) matrix(c(-0.5, 1 + 0.2 * t, 0.1 * u[1] * u[2],
                                          -1, -0.3, 0.5 * t,
                                          0.2 * u[3] * u[3], -0.4, -0.1 - 0.3 * t * t), 3, byrow = TRUE)
lin_L1 <- MatrixOperator(lin_upd1(NULL, c(1, 0.5, -0.3), NULL, 0), lin_upd1)
lin_prob1 <- ODEProblem(lin_L1, c(1, 0.5, -0.3), c(0, 1))
lin_base6 <- outer(1:6, 1:6, function(i, j) ((3 * i + 5 * j) %% 7 - 3) / 10) - diag(6)
lin_sym6 <- (lin_base6 + t(lin_base6)) / 2
lin_upd2 <- function(A, u, p, t) {
  A <- lin_base6
  for (i in 1:6) A[i, i] <- A[i, i] - 0.2 * t
  A[1, 6] <- A[1, 6] + 0.1 * u[1] * u[1]
  A
}
lin_upd3 <- function(A, u, p, t) {
  A <- lin_sym6
  for (i in 1:6) A[i, i] <- A[i, i] - 0.2 * t
  A[1, 6] <- A[1, 6] + 0.05 * u[1] * u[1]
  A[6, 1] <- A[6, 1] + 0.05 * u[1] * u[1]
  A
}
lin_u06 <- 1 - 0.1 * (0:5)
lin_prob2 <- ODEProblem(MatrixOperator(lin_base6, lin_upd2), lin_u06, c(0, 1))
lin_prob3 <- ODEProblem(MatrixOperator(lin_sym6, lin_upd3), lin_u06, c(0, 1))

# the final state, the statistics (nf, naccept, nreject), the number of
# saved points and the final time against the reference
lin_check <- function(sol, nm, tol = 1e-12) {
  ref <- .LIN_REF[[nm]]
  expect_length(sol$u, ref$nt)
  expect_equal(sol$t[[ref$nt]], ref$tend, tolerance = 1e-14, label = nm)
  expect_equal(as.numeric(sol$u[[ref$nt]]), ref$u, tolerance = tol, label = nm)
  expect_equal(c(sol$stats$nf, sol$stats$naccept, sol$stats$nreject), ref$stats, label = nm)
}

lin_algs <- c("MagnusMidpoint", "MagnusLeapfrog", "LieEuler", "MagnusGauss4", "MagnusNC6", "MagnusGL6",
              "MagnusGL8", "MagnusNC8", "MagnusGL4", "RKMK2", "RKMK4", "LieRK4", "CG2", "CG3", "CG4a")

test_that("Lie group and Magnus methods match OrdinaryDiffEqLinear", {
  for (nm in lin_algs) lin_check(solve(lin_prob1, get(nm)(), dt = 0.1), paste0("p1_", nm))
  lin_check(solve(lin_prob1, MagnusGL8(), dt = 0.25), "p1dt_MagnusGL8")
  lin_check(solve(lin_prob1, CG3(), dt = 0.25), "p1dt_CG3")
  lin_check(solve(lin_prob1, LieEuler(), dt = 0.25), "p1dt_LieEuler")
  lin_check(solve(lin_prob2, MagnusGL4(), dt = 0.1), "p2_MagnusGL4")
  lin_check(solve(lin_prob2, RKMK4(), dt = 0.1), "p2_RKMK4")
})

test_that("Krylov variants match OrdinaryDiffEqLinear", {
  for (nm in lin_algs[1:12])
    lin_check(solve(lin_prob2, get(nm)(krylov = TRUE, m = 4, iop = 2), dt = 0.1), paste0("p2k_", nm))
  lin_check(solve(lin_prob1, MagnusMidpoint(krylov = TRUE, m = 2, iop = 1), dt = 0.1), "p1k_MagnusMidpoint")
  lin_check(solve(lin_prob1, RKMK4(krylov = TRUE, m = 2, iop = 1), dt = 0.1), "p1k_RKMK4")
  lin_check(solve(lin_prob1, MagnusGL6(krylov = TRUE, m = 2, iop = 1), dt = 0.1), "p1k_MagnusGL6")
  lin_check(solve(lin_prob2, MagnusNC8(krylov = TRUE), dt = 0.2), "p2kd_MagnusNC8")
  lin_check(solve(lin_prob2, LieRK4(krylov = TRUE), dt = 0.2), "p2kd_LieRK4")
  # symmetric operators use Lanczos
  for (nm in c("MagnusMidpoint", "LieEuler", "MagnusGL6", "MagnusGauss4"))
    lin_check(solve(lin_prob3, get(nm)(krylov = TRUE, m = 4), dt = 0.1), paste0("p3k_", nm))
})

test_that("MagnusAdapt4 matches OrdinaryDiffEqLinear", {
  lp <- ODEProblem(lin_L1, c(1, 0.5, -0.3), c(0, 3))
  lin_check(solve(lp, MagnusAdapt4()), "a1_default")
  lin_check(solve(lin_prob1, MagnusAdapt4(), dt = 0.1, adaptive = FALSE), "a1_fixed")
  # the embedded error estimate is the difference of two nearly equal
  # exponentials, near rounding level for these smooth problems: the step
  # counts agree, the steps (and states) only to the last bits of exp
  lin_check(solve(lp, MagnusAdapt4(), reltol = 1e-7, abstol = 1e-9), "a1_tight", tol = 1e-7)
  lin_check(solve(ODEProblem(MatrixOperator(lin_base6, lin_upd2), lin_u06, c(0, 2)), MagnusAdapt4(),
                  reltol = 1e-5, abstol = 1e-7), "a2_default", tol = 1e-5)
})

test_that("LinearExponential matches OrdinaryDiffEqLinear", {
  pc <- ODEProblem(MatrixOperator(lin_base6), lin_u06, c(0, 2))
  lin_check(solve(pc, LinearExponential()), "le_nodt")
  lin_check(solve(pc, LinearExponential(), dt = 0.5), "le_off")
  lin_check(solve(pc, LinearExponential(krylov = "simple"), dt = 0.5), "le_simple")
  lin_check(solve(pc, LinearExponential(krylov = "simple", m = 4, iop = 2), dt = 0.5), "le_simple4")
  lin_check(solve(pc, LinearExponential(krylov = "adaptive"), dt = 0.5), "le_adapt")
  lin_check(solve(pc, LinearExponential(krylov = "adaptive", m = 3), dt = 0.5, reltol = 1e-8), "le_adapt3")
  lin_check(solve(pc, LinearExponential(krylov = "adaptive", m = 3, iop = 1)), "le_adapt3n")
  ps <- ODEProblem(MatrixOperator(lin_sym6), lin_u06, c(0, 2))
  lin_check(solve(ps, LinearExponential(krylov = "simple", m = 4), dt = 0.5), "le_sym")
  expect_equal(solve(pc, LinearExponential())$u[[2]],
               as.numeric(.eu_exp_base(2 * lin_base6) %*% lin_u06), tolerance = 1e-12)
})

test_that("CayleyEuler matches OrdinaryDiffEqLinear on a matrix state", {
  updc <- function(A, U, p, t) matrix(c(0, 1 + 0.1 * t, 0.2 * U[1, 1] * U[2, 2],
                                        -1 - 0.1 * t, 0, 0.3,
                                        -0.2 * U[1, 1] * U[2, 2], -0.3, -0.05 * U[3, 3] * U[3, 3]),
                                      3, byrow = TRUE)
  U0 <- matrix(c(1, 0.2, 0, 0.1, 0.9, 0.3, 0, -0.2, 1.1), 3, byrow = TRUE)
  L <- MatrixOperator(updc(NULL, U0, NULL, 0), updc)
  lin_check(solve(ODEProblem(L, U0, c(0, 1)), CayleyEuler(), dt = 0.1), "cay")
  lin_check(solve(SplitODEProblem(L, function(u, p, t) 0.1 * u * u, U0, c(0, 1)), CayleyEuler(), dt = 0.25),
            "cays")
})

test_that("linear methods need an operator problem", {
  prob <- ODEProblem(function(u, p, t) -u, c(1, 2), c(0, 1))
  expect_error(solve(prob, MagnusGL4(), dt = 0.1), "SciML operator")
  expect_error(solve(lin_prob1, MagnusGL4()), "give dt")
  expect_error(solve(lin_prob1, CayleyEuler(), dt = 0.1), "matrix state")
  expect_s3_class(LinearExponential(krylov = "adaptive"), "ODEAlgorithm")
  expect_error(LinearExponential(krylov = "other"))
})
