# Reference values from tools/globaldiffeq-probe.jl (GlobalDiffEq 1.7.3,
# OrdinaryDiffEqCore 4.19.1, Richardson 1.4).

glee_lv <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1] * u[2], -p[3] * u[2] + p[4] * u[1] * u[2])
glee_prob <- function() ODEProblem(glee_lv, c(1, 1), c(0, 5), c(1.5, 1, 3, 1))
glee_last <- function(sol) sol$u[[length(sol$t)]]

glee_ref <- list(
  GLEE23 = list(
    stats = c(500, 1503, 499, 1), t2 = 6.338317814055903e-5,
    u = c(6.098532177271675, 0.6281538372294413), ge = c(-3.527693513049049e-5, -1.4501661371011646e-5),
    ge3 = c(5.7857413327643676e-9, 7.610220820629187e-9),
    tq = c(1.2491649391473465, 0.5668861439183221, 6.190103474556511, 0.6635994249230687,
           2.1230998278662687, 0.2927115724388852),
    fixed_nf = 121, fixed_u = c(6.102742616247266, 0.646652297679813),
    fixed_ge = c(0.003927860581505887, -0.011254546253699132),
    save_u = c(1.53119270770476, 0.41012572415744797, 2.178460795544732, 4.368604667013049),
    save_ge = c(4.123910434218749e-5, 3.0442221856893907e-5, -0.00016887729355400183, 0.0005597953250991538)),
  GLEE24 = list(
    stats = c(631, 3178, 630, 5), t2 = 6.338317814055903e-5,
    u = c(6.098471744744195, 0.6281380442401729), ge = c(2.4161966307960977e-5, 3.7448331408812324e-6),
    ge3 = c(-5.355749355112079e-9, -7.353083965875359e-9),
    tq = c(1.2491808487449707, 0.5669023657823127, 6.190165007715442, 0.6635637067886428,
           2.123101018423962, 0.29274162546532667),
    fixed_nf = 201, fixed_u = c(6.119492219195065, 0.6408355489315587),
    fixed_ge = c(-0.014562423643726972, -0.015130133142662473),
    save_u = c(1.5312853732292737, 0.4102027948439936, 2.1781889856409498, 4.369793298783503),
    save_ge = c(-5.128446866002132e-5, -4.521419668100974e-5, 0.00014109151011413567, -0.0005964569411830037)),
  GLEE35 = list(
    stats = c(421, 2523, 420, 0), t2 = 0.0015897030511203302,
    u = c(6.098476570315634, 0.6280843911165669), ge = c(1.6774412733692566e-5, 4.955501093142941e-5),
    ge3 = c(1.192577453890814e-9, 2.7734131740307644e-9),
    tq = c(1.2491704914900312, 0.566890794755631, 6.190180564066755, 0.6635702499832633,
           2.1230415153625093, 0.29270917553810105),
    fixed_nf = 241, fixed_u = c(6.098425761072058, 0.6254240902832877),
    fixed_ge = c(-0.0010149549187692315, 0.0025497450350950167),
    save_u = c(1.5312185839356685, 0.41014219508620114, 2.178113049148977, 4.369630352051817),
    save_ge = c(1.1663435137122444e-5, 1.3105223008765816e-5, 0.00021996978218867468, -0.00036917178073861526)),
  MM5GEE = list(
    stats = c(85, 777, 84, 2), t2 = 0.020930982021551248,
    u = c(6.098500128247746, 0.6281442812963458), ge = c(-4.909067643397926e-6, -6.701832667716531e-6),
    ge3 = c(-2.5270035669024138e-9, 7.116740391444054e-10),
    tq = c(1.2491703481656447, 0.566890249680637, 6.190149873107869, 0.6635779298233053,
           2.1230944957789943, 0.2927210261249001),
    fixed_nf = 361, fixed_u = c(6.0985424048175565, 0.6281451329478511),
    fixed_ge = c(-0.000100354657695233, -0.00014122616628414029),
    save_u = c(1.531228403469777, 0.4101505430010504, 2.1783467046493685, 4.369188060091519),
    save_ge = c(-2.2300872341714513e-7, 8.620414984344388e-8, 1.399445585645544e-5, 6.935133973589903e-5)))

for (nm in names(glee_ref)) {
  test_that(paste(nm, "matches GlobalDiffEq (adaptive, fixed step, saveat)"), {
    r <- glee_ref[[nm]]
    alg <- get(nm)()
    sol <- solve(glee_prob(), alg, reltol = 1e-6, abstol = 1e-8)
    expect_identical(sol$retcode, "Success")
    expect_equal(c(length(sol$t), sol$stats$nf, sol$stats$naccept, sol$stats$nreject), r$stats)
    expect_equal(sol$t[2], r$t2, tolerance = 1e-13)
    expect_equal(glee_last(sol), r$u, tolerance = 1e-12)
    expect_equal(global_error_estimate(sol)[[length(sol$t)]], r$ge, tolerance = 1e-10)
    expect_equal(global_error_estimate(sol, 3), r$ge3, tolerance = 1e-10)
    expect_equal(as.numeric(t(sol(c(0.3, 1.7, 4.1)))), r$tq, tolerance = 1e-12)

    solf <- solve(glee_prob(), alg, dt = 0.125, adaptive = FALSE)
    expect_equal(c(length(solf$t), solf$stats$nf), c(41, r$fixed_nf))
    expect_equal(glee_last(solf), r$fixed_u, tolerance = 1e-13)
    expect_equal(solf$global_error[[41]], r$fixed_ge, tolerance = 1e-12)

    sols <- solve(glee_prob(), alg, reltol = 1e-5, abstol = 1e-7, saveat = c(0.5, 2.5))
    expect_equal(sols$t, c(0.5, 2.5))
    expect_equal(unlist(sols$u), r$save_u, tolerance = 1e-12)
    expect_equal(unlist(sols$global_error), r$save_ge, tolerance = 1e-10)
  })
}

test_that("GLEE24 and MM5GEE match GlobalDiffEq on a length one state", {
  dec <- function(u, p, t) -p[1] * u[1] * u[1]
  probd <- ODEProblem(dec, 2, c(0, 3), 0.7)
  sol <- solve(probd, GLEE24(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(c(length(sol$t), sol$stats$nf, sol$stats$nreject), c(647, 3808, 115))
  expect_equal(glee_last(sol), 0.38461607528409103, tolerance = 1e-12)
  expect_equal(sol$global_error[[647]], -6.890807211633085e-7, tolerance = 1e-9)
  sol <- solve(probd, MM5GEE(), reltol = 1e-8, abstol = 1e-10)
  expect_equal(c(length(sol$t), sol$stats$nf, sol$stats$nreject), c(31, 273, 0))
  expect_equal(glee_last(sol), 0.3846153896716246, tolerance = 1e-12)
  expect_equal(sol$global_error[[31]], -5.082657609893626e-9, tolerance = 1e-9)
  # the estimate tracks the true global error 1 / (1 / 2 + 0.7 t) - y
  expect_equal(sol$global_error[[31]], 1 / (0.5 + 0.7 * 3) - glee_last(sol), tolerance = 0.05)
})

test_that("global_error_estimate needs a GLEE solution", {
  sol <- solve(glee_prob(), Tsit5())
  expect_error(global_error_estimate(sol), "GLEE method")
  # has_global_error_glee of tools/globaldiffeq-probe.jl
  expect_true(has_global_error(GLEE24()))
  expect_false(has_global_error(Tsit5()))
})

test_that("GlobalRichardson matches GlobalDiffEq", {
  cases <- list(list(Tsit5(), 0.5, c(161, 961), 0.03125, c(6.098494682229677, 0.6281379037829629)),
                list(Euler(), 0.25, c(5121, 5121), 0.0009765625, c(6.0912407993082285, 0.6078230977114201)),
                list(RK4(), 0.25, c(641, 2561), 0.0078125, c(6.098494664428515, 0.6281378942651421)))
  for (cs in cases) {
    sol <- solve(glee_prob(), GlobalRichardson(cs[[1]]), dt = cs[[2]], reltol = 1e-6, abstol = 1e-8)
    expect_equal(c(length(sol$t), sol$stats$nf), cs[[3]])
    expect_equal(sol$t[2], cs[[4]])
    expect_equal(glee_last(sol), cs[[5]], tolerance = 1e-13)
    expect_equal(dim(sol$extrapolated), c(2L, 11L + 10L * (cs[[2]] == 0.25)))
  }
  sol <- solve(glee_prob(), GlobalRichardson(Tsit5()), dt = 0.5)
  expect_equal(length(sol$t), 41)
  expect_equal(glee_last(sol), c(6.098526634267035, 0.6281475266705442), tolerance = 1e-13)
  expect_error(solve(glee_prob(), GlobalRichardson(Tsit5())), "dt")
})
