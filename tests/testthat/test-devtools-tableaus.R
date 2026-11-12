# Tableau tools of DiffEqDevTools.jl (src/tableau_info.jl). Reference values
# from tools/DiffEqDevTools-probe.jl (section "Tableaus").

rk4_tableau <- function() {
  A <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0), 4, byrow = TRUE)
  ExplicitRKTableau(A, c(0, 0.5, 0.5, 1), c(1 / 6, 1 / 3, 1 / 3, 1 / 6), 4)
}

test_that("residual_order_condition matches DiffEqDevTools", {
  tab <- constructTsit5ExplicitRK()
  rk4 <- rk4_tableau()
  expect_equal(residual_order_condition(tab, 3), c(-9.992007221626409e-16, -2.1371793224034263e-15),
    tolerance = 1e-14)
  expect_equal(residual_order_condition(tab, 4),
    c(-1.249000902703301e-16, -6.245004513516506e-17, -2.1094237467877974e-15, -1.0362081563168126e-15),
    tolerance = 1e-14)
  res <- vapply(1:6, function(i) residual_order_condition(tab, i, `+`, abs), 0)
  expect_equal(res[5:6], c(0.002301776475535945, 0.006074815959888372), tolerance = 1e-12)
  expect_true(all(res[1:4] < 1e-14))
  emb <- vapply(1:5, function(i) residual_order_condition(tab, i, `+`, abs, embedded = TRUE), 0)
  expect_true(all(emb < 1e-14))
  expect_equal(vapply(1:5, function(i) residual_order_condition(rk4, i, `+`, abs), 0),
    c(1.1102230246251565e-16, 0, 0, 0, 0.03506944444444444), tolerance = 1e-14)
  expect_equal(vapply(1:5, function(i) residual_order_condition(rk4, i, `+`), 0)[5],
    0.00021038290895061725, tolerance = 1e-14)
})

test_that("check_tableau matches DiffEqDevTools", {
  expect_error(check_tableau(constructTsit5ExplicitRK()), "Tableau's order is not correct.", fixed = TRUE)
  expect_true(check_tableau(rk4_tableau()))
  # Deliberate difference: Julia iterates over the integer adaptiveorder and
  # checks only the conditions of that order (probe line check_embedded
  # returns true although sum(bhat) = 1/2); R checks every order up to it.
  rk4 <- rk4_tableau()
  bad <- ExplicitRKTableau(rk4$A, rk4$c, rk4$alpha, 4, alphaEEst = c(0, 0, 0, 0.5), adaptiveorder = 2)
  expect_error(check_tableau(bad), "Tableau's embedded order is not correct.", fixed = TRUE)
})

test_that("stability_region matches DiffEqDevTools", {
  tab <- constructTsit5ExplicitRK()
  rk4 <- rk4_tableau()
  expect_equal(stability_region(-1, tab), 0.3674947481323657, tolerance = 1e-13)
  expect_equal(stability_region(-1, tab, embedded = TRUE), 0.36809887799147256, tolerance = 1e-13)
  expect_equal(stability_region(-2.5, rk4), 0.6484375000000001, tolerance = 1e-14)
  expect_equal(stability_region(0.3, tab), 1.3498597214754864, tolerance = 1e-14)
  expect_equal(stability_region(-1 + 2i, tab), -0.16879264946414585 + 0.2836711935050898i, tolerance = 1e-12)
  expect_equal(stability_region(0.5i, rk4), 0.8776041666666666 + 0.47916666666666663i, tolerance = 1e-14)
  expect_equal(stability_region(tab), -4.056036047197509, tolerance = 1e-12)
  expect_equal(stability_region(rk4), -2.785293563405282, tolerance = 1e-12)
  algs <- list(Euler = c(-0.5, -2), RK4 = c(0.27343750000000006, -2.785293563405282),
    Tsit5 = c(0.22647003212163064, -3.506846993804954), SSPRK33 = c(0.0625, -2.5127453266183286),
    Midpoint = c(0.625, -2), Heun = c(0.625, -2), BS3 = c(0.0625, -2.5127453266183286),
    DP5 = c(0.22914062499999943, -3.3065678926351065))
  for (nm in names(algs)) {
    alg <- get(nm)()
    expect_equal(c(stability_region(-1.5, alg), stability_region(alg)), algs[[nm]], tolerance = 1e-11, info = nm)
  }
  expect_equal(stability_region(-0.5 + 1i, Tsit5()), 0.32787867799483106 + 0.5099528547141972i, tolerance = 1e-12)
  expect_equal(stability_region(-2, ImplicitEuler()), 0.33333333333333337, tolerance = 1e-14)
})

test_that("imaginary_stability_interval matches DiffEqDevTools", {
  expect_equal(imaginary_stability_interval(rk4_tableau()), 2.8284271247461903, tolerance = 1e-14)
  expect_equal(imaginary_stability_interval(RK4()), 2.8284271247461903, tolerance = 1e-14)
  expect_equal(imaginary_stability_interval(SSPRK33()), 1.732050807568905, tolerance = 1e-14)
  # |R(iy)| - 1 is flat near the root of Tsit5 (slope about 2e-9), so both
  # roots are within the solver tolerance of each other.
  expect_equal(imaginary_stability_interval(Tsit5()), 0.4779788742764177, tolerance = 1e-6)
  # Julia throws an ArgumentError for the tableau of Tsit5 as well.
  expect_error(imaginary_stability_interval(constructTsit5ExplicitRK()))
})

test_that("deduce_Butcher_tableau matches DiffEqDevTools", {
  r <- deduce_Butcher_tableau(RK4())
  A <- matrix(0, 5, 5)
  A[2, 1] <- 0.5; A[3, 2] <- 0.5; A[4, 3] <- 1; A[5, 1:4] <- c(1 / 6, 1 / 3, 1 / 3, 1 / 6)
  expect_equal(r$A, A, tolerance = 1e-15, ignore_attr = TRUE)
  expect_equal(r$b, c(1 / 6, 1 / 3, 1 / 3, 1 / 6, 0), tolerance = 1e-15)
  expect_equal(r$c, c(0, 0.5, 0.5, 1, 1), tolerance = 1e-15)
  r <- deduce_Butcher_tableau(SSPRK33())
  expect_equal(r$A, matrix(c(0, 0, 0, 1, 0, 0, 0.25, 0.25, 0), 3, byrow = TRUE), tolerance = 1e-15,
    ignore_attr = TRUE)
  expect_equal(r$c, c(0, 1, 0.5), tolerance = 1e-15)
  r <- deduce_Butcher_tableau(BS3())
  expect_equal(r$b, c(2 / 9, 1 / 3, 4 / 9, 0), tolerance = 1e-15)
  expect_equal(r$c, c(0, 0.5, 0.75, 1), tolerance = 1e-15)
  r <- deduce_Butcher_tableau(ExplicitRK(constructTsit5ExplicitRK()))
  expect_equal(r$A[7, 1:3], c(0.09646076681806523, 0.01, 0.4798896504144996), tolerance = 1e-14)
  expect_equal(r$A[3, 1:2], c(-0.008480655492356989, 0.335480655492357), tolerance = 1e-14)
})
