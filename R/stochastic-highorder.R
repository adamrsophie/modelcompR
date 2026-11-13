# The configurable Roessler methods of StochasticDiffEqHighOrder:
# tableaus.jl (RosslerSRI, RosslerSRA, the construct* functions,
# checkSRIOrder, checkSRAOrder), perform_step/sri.jl and sra.jl (SRI, SRA
# and the four stage SRIW2 and SOSRI2). The tableau constants are generated
# into R/stochastic-tableaus.R by tools/sde-gen-tableaus.jl.

#' Roessler SRI and SRA methods
#'
#' `SRI(tableau, error_terms)` is the strong order 1.5 stochastic Runge-Kutta
#' method of Roessler for diagonal or scalar Ito noise with any
#' `RosslerSRI()` tableau (`constructSRIW1()`, `constructSRIW2()`,
#' `constructSRIOpt1()`, `constructSRIOpt2()`); `error_terms` drift stages
#' enter the error estimate. `SRIW2()` and `SOSRI2()` are the four stage
#' methods with the SRIW2 and the second stability optimized tableau.
#' `SRA(tableau)` is the additive noise method with any `RosslerSRA()`
#' tableau (`constructSRA1()`, `constructSRA2()`, `constructSRA3()`,
#' `constructSOSRA()`, `constructSOSRA2()`; `constructSKenCarp()` and
#' `constructExplicitSKenCarp()` are the KenCarp pair as SRA tableaus).
#' `checkSRIOrder(tab)` and `checkSRAOrder(tab)` evaluate the strong order
#' 1.5 conditions (25 and 8 of them) and return which hold.
#'
#' @param tableau A `RosslerSRI` or `RosslerSRA` tableau.
#' @param error_terms Number of drift evaluations in the error estimate.
#' @param c0,c1,A0,A1,B0,B1,alpha,beta1,beta2,beta3,beta4 Tableau entries.
#' @param order Strong order.
#' @param tol Tolerance of the order conditions.
#' @return An algorithm, a tableau or a logical vector.
#' @examples
#' all(checkSRIOrder(constructSRIW1()))
#' prob <- SDEProblem(function(u, p, t) 1.01 * u, function(u, p, t) 0.87 * u, 0.5, c(0, 1))
#' solve(prob, SRI(constructSRIW2()), seed = 1)$retcode
#' @name rossler_methods
NULL

#' @rdname rossler_methods
#' @export
RosslerSRI <- function(c0, c1, A0, A1, B0, B1, alpha, beta1, beta2, beta3, beta4, order = 1.5)
  structure(list(c0 = c0, c1 = c1, A0 = as.matrix(A0), A1 = as.matrix(A1), B0 = as.matrix(B0), B1 = as.matrix(B1),
                 alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4, order = order),
            class = c("RosslerSRI", "Tableau"))

#' @rdname rossler_methods
#' @export
RosslerSRA <- function(c0, c1, A0, B0, alpha, beta1, beta2, order = 1.5)
  structure(list(c0 = c0, c1 = c1, A0 = as.matrix(A0), B0 = as.matrix(B0), alpha = alpha, beta1 = beta1,
                 beta2 = beta2, order = order), class = c("RosslerSRA", "Tableau"))

.sde_tab <- function(name, cls) {
  tb <- .SDE_TABLEAUS[[name]]
  if (cls == "RosslerSRI")
    return(RosslerSRI(tb$c0, tb$c1, tb$A0, tb$A1, tb$B0, tb$B1, tb$alpha, tb$beta1, tb$beta2, tb$beta3, tb$beta4,
                      tb$order %||% 1.5))
  RosslerSRA(tb$c0, tb$c1, tb$A0, tb$B0, tb$alpha, tb$beta1, tb$beta2, tb$order %||% 1.5)
}

#' @rdname rossler_methods
#' @export
constructSRIW1 <- function() .sde_tab("SRIW1", "RosslerSRI")
#' @rdname rossler_methods
#' @export
constructSRIW2 <- function() .sde_tab("SRIW2", "RosslerSRI")
#' @rdname rossler_methods
#' @export
constructSRIOpt1 <- function() .sde_tab("SRIOpt1", "RosslerSRI")
#' @rdname rossler_methods
#' @export
constructSRIOpt2 <- function() .sde_tab("SRIOpt2", "RosslerSRI")
#' @rdname rossler_methods
#' @export
constructSRA1 <- function() .sde_tab("SRA1", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructSRA2 <- function() .sde_tab("SRA2", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructSRA3 <- function() .sde_tab("SRA3", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructSOSRA <- function() .sde_tab("SOSRA", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructSOSRA2 <- function() .sde_tab("SOSRA2", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructSKenCarp <- function() .sde_tab("SKenCarp", "RosslerSRA")
#' @rdname rossler_methods
#' @export
constructExplicitSKenCarp <- function() .sde_tab("ExplicitSKenCarp", "RosslerSRA")

#' @rdname rossler_methods
#' @export
checkSRIOrder <- function(tableau, tol = 1e-6) {
  t <- tableau
  e <- rep(1, length(t$alpha))
  s <- function(b, M) sum(b %*% M)
  sq <- function(M) as.numeric(M %*% e)^2
  B1e <- as.numeric(t$B1 %*% e)
  c(abs(sum(t$alpha) - 1) < tol, abs(sum(t$beta1) - 1) < tol, abs(sum(t$beta2)) < tol, abs(sum(t$beta3)) < tol,
    abs(sum(t$beta4)) < tol,
    abs(s(t$beta1, t$B1)) < tol, abs(s(t$beta2, t$B1) - 1) < tol, abs(s(t$beta3, t$B1)) < tol,
    abs(s(t$beta4, t$B1)) < tol,
    abs(s(t$alpha, t$A0) - 0.5) < tol, abs(s(t$alpha, t$B0) - 1) < tol, abs(sum(t$alpha * sq(t$B0)) - 3 / 2) < tol,
    abs(s(t$beta1, t$A1) - 1) < tol, abs(s(t$beta2, t$A1)) < tol, abs(s(t$beta3, t$A1) + 1) < tol,
    abs(s(t$beta4, t$A1)) < tol,
    abs(sum(t$beta1 * sq(t$B1)) - 1) < tol, abs(sum(t$beta2 * sq(t$B1))) < tol, abs(sum(t$beta3 * sq(t$B1)) + 1) < tol,
    abs(sum(t$beta4 * sq(t$B1)) - 2) < tol,
    abs(sum(t$beta1 %*% t$B1 %*% B1e)) < tol, abs(sum(t$beta2 %*% t$B1 %*% B1e)) < tol,
    abs(sum(t$beta3 %*% t$B1 %*% B1e)) < tol, abs(sum(t$beta4 %*% t$B1 %*% B1e) - 1) < tol,
    abs(0.5 * sum(t$beta1 * (t$A1 %*% (t$B0 %*% e))) + (1 / 3) * sum(t$beta3 * (t$A1 %*% (t$B0 %*% e)))) < tol)
}

#' @rdname rossler_methods
#' @export
checkSRAOrder <- function(tableau, tol = 1e-6) {
  t <- tableau
  e <- rep(1, length(t$alpha))
  c(abs(sum(t$alpha) - 1) < tol, abs(sum(t$beta1) - 1) < tol, abs(sum(t$beta2)) < tol,
    abs(sum(t$alpha * (t$B0 %*% e)) - 1) < tol, abs(sum(t$alpha * (t$A0 %*% e)) - 1 / 2) < tol,
    abs(sum(t$alpha * (t$B0 %*% e)^2) - 3 / 2) < tol, abs(sum(t$beta1 * t$c1) - 1) < tol,
    abs(sum(t$beta2 * t$c1) + 1) < tol)
}

# perform_step! for SRIConstantCache
.sri_generic_step <- function(tab, error_terms) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  s <- length(tab$alpha)
  sq <- .sqdt(dt)
  chi1 <- (dW^2 - abs(dt)) / (2 * sq)
  chi2 <- (dW + dZ / sqrt(3)) / 2
  chi3 <- (dW^3 - 3 * dW * dt) / (6 * dt)
  H0 <- vector("list", s); H1 <- vector("list", s)
  fk <- vector("list", s); gk <- vector("list", s)
  for (i in seq_len(s)) {
    A0t <- 0 * u; A1t <- 0 * u; B0t <- 0 * u; B1t <- 0 * u
    for (j in seq_len(i - 1L)) {
      A0t <- A0t + tab$A0[i, j] * fk[[j]]; A1t <- A1t + tab$A1[i, j] * fk[[j]]
      B0t <- B0t + tab$B0[i, j] * gk[[j]]; B1t <- B1t + tab$B1[i, j] * gk[[j]]
    }
    H0[[i]] <- u + dt * A0t + chi2 * B0t
    H1[[i]] <- u + dt * A1t + sq * B1t
    fk[[i]] <- f(H0[[i]], p, t + tab$c0[i] * dt)
    gk[[i]] <- g(H1[[i]], p, t + tab$c1[i] * dt)
  }
  at <- 0 * u; bt <- 0 * u; E2 <- 0 * u; E1 <- 0 * u
  for (i in seq_len(s)) {
    at <- at + tab$alpha[i] * fk[[i]]
    bt <- bt + (tab$beta1[i] * dW + tab$beta2[i] * chi1) * gk[[i]]
    E2 <- E2 + (tab$beta3[i] * chi2 + tab$beta4[i] * chi3) * gk[[i]]
    if (i <= error_terms) E1 <- E1 + fk[[i]]
  }
  list(u = u + dt * at + bt + E2, err = delta * dt * E1 + E2, nf = s)
}

# perform_step! for SRAConstantCache
.sra_generic_step <- function(tab) function(f, g, u, p, t, dt, dW, dZ, type, delta) {
  s <- length(tab$alpha)
  chi2 <- (dW + dZ / sqrt(3)) / 2
  H0 <- vector("list", s); fk <- vector("list", s); gk <- vector("list", s)
  for (i in seq_len(s)) {
    A0t <- 0 * u; B0t <- 0 * u
    for (j in seq_len(i - 1L)) {
      A0t <- A0t + tab$A0[i, j] * fk[[j]]
      B0t <- B0t + .gdw(tab$B0[i, j] * gk[[j]], chi2, type)
    }
    H0[[i]] <- u + A0t * dt + B0t
    fk[[i]] <- f(H0[[i]], p, t + tab$c0[i] * dt)
    gk[[i]] <- g(H0[[i]], p, t + tab$c1[i] * dt)
  }
  at <- 0 * u; bt <- 0 * u; E2 <- 0 * u; E1 <- 0 * u
  for (i in seq_len(s)) {
    E1 <- E1 + fk[[i]]
    at <- at + tab$alpha[i] * fk[[i]]
    bt <- bt + .gdw(tab$beta1[i] * gk[[i]], dW, type)
    E2 <- E2 + .gdw(tab$beta2[i] * gk[[i]], chi2, type)
  }
  list(u = u + dt * at + bt + E2, err = delta * dt * E1 + E2, nf = s)
}

#' @rdname rossler_methods
#' @export
SRI <- function(tableau = constructSRIW1(), error_terms = 4L)
  .new_sde_alg("SRI", tableau$order %||% 1.5, TRUE, .sri_generic_step(tableau, as.integer(error_terms)),
               needs_dZ = TRUE, noise = c("diagonal", "scalar"), tableau = tableau)

#' @rdname rossler_methods
#' @export
SRIW2 <- function()
  .new_sde_alg("SRIW2", 1.5, TRUE, .sri_generic_step(constructSRIW2(), 4L), needs_dZ = TRUE,
               noise = c("diagonal", "scalar"))

#' @rdname rossler_methods
#' @export
SOSRI2 <- function()
  .new_sde_alg("SOSRI2", 1.5, TRUE, .four_stage_sri(.SOSRI2), needs_dZ = TRUE, noise = c("diagonal", "scalar"),
               stability_size = 10.6)

# SOSRI2ConstantCache of caches/rossler_caches.jl
.SOSRI2 <- list(
  a021 = 0.13804532298278663, a031 = 0.5818361298250374, a032 = 0.4181638701749618,
  a041 = 0.4670018408674211, a042 = 0.8046204792187386, a043 = -0.27162232008616016,
  a121 = 0.45605532163856893, a131 = 0.7555807846451692, a132 = 0x1.f4920fc3b4257p-3,
  a141 = 0.6981181143266059, a142 = 0.3453277086024727, a143 = -0.04344582292908241,
  b021 = 0.08852381537667678, b031 = 0x1.08226c29f34f5p+0, b032 = 0.4563552922077882,
  b041 = 1.73078280444124, b042 = -0.46089678470929774, b043 = -0.9637509618944188,
  b121 = 0.6753186815412179, b131 = -0.07452812525785148, b132 = -0.49783736486149366,
  b141 = -0.5591906709928903, b142 = 0.022696571806569924, b143 = -0.8984927888368557,
  a1 = -0.15036858140642623, a2 = 0.7545275856696072, a3 = 0.686995463807979, a4 = -0.2911544680711602,
  c02 = 0.13804532298278663, c03 = 0.9999999999999992, c04 = 0.9999999999999994,
  c11 = 0.0, c12 = 0.45605532163856893, c13 = 0.999999999999996, c14 = 0.9999999999999962,
  beta11 = -0x1.d0085c96faa21p-2, beta12 = 0.8330937231303951, beta13 = 0.3792843195533544,
  beta14 = 0.24077885458934192, beta21 = -0.4994383733810986, beta22 = 0x1.d61b8205ade46p-1,
  beta23 = -0.25613778661003145, beta24 = -0.16260245862427797, beta31 = 1.4531568972730915,
  beta32 = -0.8330937231303933, beta33 = -0.3792843195533583, beta34 = -0.24077885458934023,
  beta41 = -0.4976090683622265, beta42 = 0.9148155835648892, beta43 = -1.4102107084476505,
  beta44 = 0x1.fc6b0bad8e615p-1)

#' @rdname rossler_methods
#' @export
SRA <- function(tableau = constructSRA1())
  .new_sde_alg("SRA", tableau$order %||% 1.5, TRUE, .sra_generic_step(tableau), needs_dZ = TRUE, tableau = tableau)
