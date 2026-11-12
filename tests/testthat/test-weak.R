# Reference values: one fixed noise path (NoiseGrid) through StochasticDiffEq
# 6 (StochasticDiffEqWeak), two steps of size 0.1.

weak_f <- function(u, p, t) c(-0.5 * u[1] + 0.2 * u[2] * u[1], 0.3 * u[1] - u[2])
weak_gG <- function(u, p, t) matrix(c(0.3 * u[1], 0.05 * u[1] * u[2], 0.1 * u[2], 0.2 * u[2]), 2, 2)
weak_gD <- function(u, p, t) c(0.3 * u[1] + 0.1, 0.2 * u[2]^2)
weak_ts <- c(0, 0.1, 0.2)

test_that("weak methods match StochasticDiffEq on a fixed noise path with two channels", {
  W <- list(c(0, 0), c(0.4, -0.1), c(0.38, 0.35))
  Z <- list(c(0, 0), c(0.2, -0.3), c(-0.1, -0.5))
  ref <- list(
    DRI1_diag = list("diag", c(1.13054794717022777, 0.49477598203620754)),
    DRI1_gen = list("gen", c(1.10520676601850365, 0.53214142192833092)),
    RI1_diag = list("diag", c(1.13059409086667961, 0.49469335489997079)),
    RI1_gen = list("gen", c(1.10522935522842736, 0.53212870810923241)),
    RI3_diag = list("diag", c(1.13049357225845082, 0.49489381709122737)),
    RI3_gen = list("gen", c(1.10519379354484393, 0.53215290228192069)),
    RI5_diag = list("diag", c(0.72573514361091451, 0.43100831608555873)),
    RI5_gen = list("gen", c(0.75128457637260038, 0.39765884543699842)),
    RI6_diag = list("diag", c(1.13050841331807006, 0.49484945077035591)),
    RI6_gen = list("gen", c(1.10518704336568985, 0.53215105685499386)),
    RDI1WM_diag = list("diag", c(1.12898736125595556, 0.49486758608200299)),
    RDI1WM_gen = list("gen", c(1.10474199016504993, 0.53207256215075849)),
    RDI2WM_diag = list("diag", c(1.13050826135661731, 0.49484773046900871)),
    RDI2WM_gen = list("gen", c(1.10518765399752295, 0.53215365981992568)),
    RDI3WM_diag = list("diag", c(1.13051326086509030, 0.49483219262073819)),
    RDI3WM_gen = list("gen", c(1.10519109879596766, 0.53214779871694684)),
    RDI4WM_diag = list("diag", c(1.13054793353321070, 0.49477582736947900)),
    RDI4WM_gen = list("gen", c(1.10520515898194516, 0.53214110216911059)),
    PL1WM_gen = list("gen", c(1.10518704490698472, 0.53223343339306706)),
    NON_diag = list("diag", c(1.14371599504026755, 0.49608524326280917)),
    NON_gen = list("gen", c(1.11718475417476393, 0.53559092421016707)),
    NON2_diag = list("diag", c(1.14371599504026755, 0.49608524326280917)),
    NON2_gen = list("gen", c(1.11403496907357358, 0.53516670859802651)),
    COM_diag = list("diag", c(1.07183869187013769, 0.48394101981449494)),
    COM_gen = list("gen", c(1.05134493679312802, 0.50898511845488237)),
    SIEA_diag = list("diag", c(1.05935845530199302, 0.48284866500895707)),
    SMEB_diag = list("diag", c(1.05962419573705913, 0.48125328981799831)),
    W2Ito1_diag = list("diag", c(1.13349007842532501, 0.49413803422610481)),
    W2Ito1_gen = list("gen", c(1.10769655380299503, 0.53185949803405674)),
    IRI1_diag = list("diag", c(1.13346050228391482, 0.49514563546959833)),
    IRI1_gen = list("gen", c(1.10790828241859707, 0.53314817682060023)),
    IRI1h_diag = list("diag", c(1.13351309531831190, 0.49512454225356350)),
    IRI1h_gen = list("gen", c(1.10795101168206700, 0.53315702538995580)),
    DRI1NM_diag = list("diag", c(1.13054794717022800, 0.49477598203620760))
  )
  for (nm in names(ref)) {
    alg_name <- sub("_(diag|gen)$", "", nm)
    alg <- if (alg_name == "IRI1h") IRI1(theta = 0.5) else get(alg_name)()
    nz <- NoiseGrid(weak_ts, W, Z)
    prob <- if (ref[[nm]][[1]] == "diag") SDEProblem(weak_f, weak_gD, c(1, 0.5), c(0, 0.2), noise = nz) else
      SDEProblem(weak_f, weak_gG, c(1, 0.5), c(0, 0.2), noise = nz, noise_rate_prototype = matrix(0, 2, 2))
    sol <- solve(prob, alg, dt = 0.1, adaptive = FALSE)
    # IRI1 stops its Newton iterations at a different point than Julia
    tol <- if (startsWith(alg_name, "IRI1")) 1e-6 else 1e-12
    expect_equal(sol$u[[3]], ref[[nm]][[2]], tolerance = tol, info = nm)
  }
})

test_that("weak methods match StochasticDiffEq for scalar noise and one Wiener process", {
  ref <- list(
    RS1_scal = c(0.92824329963204355),
    RS1_one = c(1.1154588648590296, 0.4826009845053858),
    RS2_scal = c(0.92826996347302237),
    RS2_one = c(1.11545450605234797, 0.48260065834337168),
    PL1WM_scal = c(0.91059616337325722),
    PL1WM_one = c(1.10388845096418486, 0.48121162909262943),
    PL1WMA_scal = c(0.90913305553199275),
    PL1WMA_one = c(1.10332201702675725, 0.48118153203253383),
    NON_scal = c(0.92835622113066363),
    NON_one = c(1.11568697539665118, 0.48246315925066136),
    NON2_scal = c(0.92835622113066363),
    NON2_one = c(1.11568697539665118, 0.48246315925066136),
    COM_scal = c(0.86753819182795022),
    COM_one = c(1.05299261500892216, 0.47580287417663047),
    DRI1_scal = c(0.91062240476615364),
    DRI1_one = c(1.1039238019770627, 0.4811439115863837),
    W2Ito1_scal = c(0.91344211813094023),
    W2Ito1_one = c(1.1072651467947836, 0.4816408448347429),
    SIEA_scal = c(0.85193698639542725),
    SMEA_scal = c(0.85238595677318341)
  )
  fs <- function(u, p, t) -0.5 * u + 0.2 * u^2
  gs <- function(u, p, t) 0.3 * u + 0.1 * u^2
  g1 <- function(u, p, t) matrix(c(0.3 * u[1] + 0.1 * u[2], 0.05 * u[1] * u[2]), 2, 1)
  for (nm in names(ref)) {
    alg <- get(sub("_(scal|one)$", "", nm))()
    prob <- if (endsWith(nm, "_scal"))
      SDEProblem(fs, gs, 0.8, c(0, 0.2), noise = NoiseGrid(weak_ts, c(0, 0.4, 0.38), c(0, 0.2, -0.1))) else
      SDEProblem(weak_f, g1, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(weak_ts, list(0, 0.4, 0.38), list(0, 0.2, -0.1)),
                 noise_rate_prototype = matrix(0, 2, 1))
    sol <- solve(prob, alg, dt = 0.1, adaptive = FALSE)
    expect_equal(as.numeric(sol$u[[3]]), ref[[nm]], tolerance = 1e-12, info = nm)
  }
})

test_that("diagonal noise is the general noise of the diagonal matrix", {
  W <- list(c(0, 0), c(0.4, -0.1), c(0.38, 0.35))
  Z <- list(c(0, 0), c(0.2, -0.3), c(-0.1, -0.5))
  gM <- function(u, p, t) diag(weak_gD(u, p, t))
  for (alg in list(RS1(), RS2(), PL1WM(), NON(), DRI1(), W2Ito1())) {
    a <- solve(SDEProblem(weak_f, weak_gD, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(weak_ts, W, Z)), alg,
               dt = 0.1, adaptive = FALSE)
    b <- solve(SDEProblem(weak_f, gM, c(1, 0.5), c(0, 0.2), noise = NoiseGrid(weak_ts, W, Z),
                          noise_rate_prototype = matrix(0, 2, 2)), alg, dt = 0.1, adaptive = FALSE)
    expect_equal(a$u[[3]], b$u[[3]], tolerance = 1e-14)
  }
})

test_that("the shipped weak tableaus satisfy their order conditions", {
  for (nm in c("DRI1", "RI1", "RI3", "RI5", "RI6", "RDI2WM", "RDI3WM", "RDI4WM"))
    expect_true(all(checkRIOrder(get(paste0("construct", nm))())), info = nm)
  expect_true(all(checkRIOrder(constructRDI1WM(), ps = 1)))
  expect_length(checkRIOrder(constructDRI1()), 59)
  expect_length(checkRIOrder(constructDRI1(), ps = 1), 9)
  for (nm in c("RS1", "RS2")) expect_true(all(checkRSOrder(get(paste0("construct", nm))())), info = nm)
  expect_length(checkRSOrder(constructRS1()), 55)
  expect_true(all(checkNONOrder(constructNON())))
  expect_true(all(checkNONOrder(constructNON2())))
  expect_s3_class(constructRS2(), "RoesslerRS")
  expect_error(checkRIOrder(constructDRI1(), ps = 3), "first and second order")
  tab <- constructRI1()
  tab$alpha <- tab$alpha * 1.1
  expect_false(all(checkRIOrder(tab)))
})

test_that("Ihat2 follows each family", {
  dW <- c(0.3, -0.2, 0.1); dZ <- c(0.5, -0.5, 0.5)
  expect_equal(Ihat2(DRI1(), dW, dZ, 0.5, 1, 2), (dW[1] * dW[2] - 0.5 * dZ[1]) / 2)
  expect_equal(Ihat2(DRI1(), dW, dZ, 0.5, 2, 1), (dW[2] * dW[1] + 0.5 * dZ[1]) / 2)
  expect_equal(Ihat2("NON", dW, dZ, 0.5, 1, 2), dZ[1])
  expect_equal(Ihat2("NON", dW, dZ, 0.5, 2, 1), dW[1])
  expect_equal(Ihat2("PL1WM", dW, dZ, 0.5, 2, 2), -0.25)
  expect_error(Ihat2("SIEA", dW, dZ, 0.5, 1, 2), "no Ihat2")
})

test_that("weak methods reproduce the moments of geometric Brownian motion", {
  a <- 1.01; b <- 0.87; u0 <- 0.5
  prob <- SDEProblem(function(u, p, t) a * u, function(u, p, t) b * u, u0, c(0, 1))
  for (alg in list(DRI1(), RI1(), W2Ito1(), SIEA(), PL1WM())) {
    sim <- solve(EnsembleProblem(prob), alg, trajectories = 4000, dt = 1 / 8, adaptive = FALSE, seed = 11)
    uT <- vapply(sim$u, function(s) s$u[[length(s$u)]], 0)
    expect_lt(abs(mean(uT) - u0 * exp(a)), 4 * sd(uT) / sqrt(length(uT)))
  }
  # Stratonovich: E u(1) = u0 exp(a + b^2 / 2)
  for (alg in list(RS1(), NON2(), COM())) {
    expect_equal(alg$interpretation, "Stratonovich")
    sim <- solve(EnsembleProblem(prob), alg, trajectories = 4000, dt = 1 / 8, adaptive = FALSE, seed = 12)
    uT <- vapply(sim$u, function(s) s$u[[length(s$u)]], 0)
    expect_lt(abs(mean(uT) - u0 * exp(a + b^2 / 2)), 4 * sd(uT) / sqrt(length(uT)))
  }
})

test_that("DRI1 adapts its step for diagonal noise and DRI1NM refuses general noise", {
  prob <- SDEProblem(weak_f, weak_gD, c(1, 0.5), c(0, 1))
  sol <- solve(prob, DRI1(), seed = 5)
  expect_equal(sol$retcode, "Success")
  expect_gt(length(sol$t), 3)
  prob2 <- SDEProblem(weak_f, weak_gG, c(1, 0.5), c(0, 1), noise_rate_prototype = matrix(0, 2, 2))
  expect_error(solve(prob2, DRI1NM(), dt = 0.1, adaptive = FALSE, seed = 1), "non mixing")
})
