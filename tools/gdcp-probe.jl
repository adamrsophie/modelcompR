# Reference values for tests/testthat/test-gdcp.R: the geodesic (gDCP)
# analysis of SymbolicAnalysis.jl 0.5.8 on the SymmetricPositiveDefinite and
# Lorentz manifolds. Run from the repository root:
#   julia tools/gdcp-probe.jl > gdcp.tsv
# One line per case: name, curvature, sign, gcurvature and the warnings
# logged (joined by " | "); a case that throws prints name, ERROR and the
# message.
using Pkg
Pkg.activate(; temp = true)
Pkg.add([
    PackageSpec(name = "SymbolicAnalysis", version = "0.5.8"),
    PackageSpec(name = "Symbolics"),
    PackageSpec(name = "Manifolds"),
]; io = devnull)
Pkg.status(; io = stderr)
using SymbolicAnalysis, Symbolics, Manifolds, LinearAlgebra, Logging
using Symbolics: unwrap
using Test: TestLogger
const SA = SymbolicAnalysis

function case(name, build, M = nothing)
    logger = TestLogger(; min_level = Logging.Warn)
    line = try
        res = with_logger(logger) do
            ex = unwrap(build())
            M === nothing ? analyze(ex) : analyze(ex, M)
        end
        g = res.gcurvature === nothing ? "nothing" : string(res.gcurvature)
        join([name, string(res.curvature), string(res.sign), g,
            join([string(l.message) for l in logger.logs], " | ")], "\t")
    catch e
        e isa InterruptException && rethrow()
        join([name, "ERROR", replace(sprint(showerror, e), "\n" => " ", "\t" => " ")], "\t")
    end
    println(line)
end

@variables X[1:3, 1:3] Y[1:3, 1:3] p[1:3] x[1:3]
S = SymmetricPositiveDefinite(3)
L = Lorentz(3)
A = [2.0 0.5 0.0; 0.5 3.0 0.25; 0.0 0.25 1.5]
B = [1.0 0.2 0.1; 0.2 2.0 0.3; 0.1 0.3 1.0]
C = [1.5 0.0 0.3; 0.0 1.0 0.0; 0.3 0.0 2.5]
ys = [[1.0, 0.5, 0.25], [0.0, 1.0, 2.0]]
y1 = [1.0, 0.5, 0.25]

# the upstream test set (test/dgp.jl and test/lorentz.jl)
case("logdet_conj_inv", () -> logdet(SA.conjugation(inv(X), A)), S)
case("logdet_tr_inv", () -> logdet(tr(inv(X))), S)
case("lqf_sum", () -> SA.log_quad_form(ys[1], inv(X)) + SA.log_quad_form(ys[2], inv(X)) + 1 / 5 * logdet(X), S)
case("brascamp", () -> logdet(SA.conjugation(X, A)) - logdet(X), S)
case("tr_conj", () -> tr(SA.conjugation(X, A)), S)
case("sdiv", () -> SA.sdivergence(X, A), S)
case("sdiv_sum", () -> SA.sdivergence(X, A) + SA.sdivergence(X, B) + SA.sdivergence(X, C), S)
case("dist_sq", () -> Manifolds.distance(S, A, X)^2, S)
case("dist_sq_sum", () -> Manifolds.distance(S, A, X)^2 + Manifolds.distance(S, B, X)^2, S)
case("sqrt_XY", () -> sqrt(X * Y), S)
case("matsqrt", () -> SA.sdivergence(X, A) + SA.sdivergence(X, Matrix{Float64}(I(3))), S)
case("tr_inv_logdet", () -> tr(inv(X)) + logdet(X), S)
case("tr_shift_sq", () -> (tr(X + 0.5 * I(3)))^2, S)
case("lqf_inv", () -> SA.log_quad_form(x, inv(X)), S)
case("lqf_ys", () -> SA.log_quad_form(ys, X), S)
case("lqf_ys_inv", () -> SA.log_quad_form(ys, inv(X)), S)
case("eigsummax_log", () -> sum(SA.eigsummax(log(X), 2)), S)
case("schatten_log", () -> sum(SA.schatten_norm(log(X), 3)), S)
case("exp_eigsummax_log", () -> exp(SA.eigsummax(log(X), 2)), S)
case("sle", () -> SA.sum_log_eigmax(X, 2), S)
case("sle_exp", () -> SA.sum_log_eigmax(exp, X, 2), S)
case("tr_affine_conj", () -> tr(SA.affine_map(SA.conjugation, X, B, C)), S)
case("hadamard", () -> SA.hadamard_product(X, B), S)
case("logdet_affine_hadamard", () -> logdet(SA.affine_map(SA.hadamard_product, X, A, B)), S)
case("tr_sign", () -> tr(X), S)
case("eigmax_sign", () -> eigmax(X), S)
case("tr_sign_euclid", () -> tr(X))
case("eigmin", () -> eigmin(X), S)
case("sqrt", () -> sqrt(X), S)
case("norm", () -> norm(X), S)
case("tr", () -> tr(X), S)
case("eigmax", () -> eigmax(X), S)
case("logdet", () -> logdet(X), S)
case("log_tr", () -> log(tr(X)), S)
case("abs_logdet_tr", () -> abs(logdet(X) + tr(X)), S)
case("lqf_sign", () -> SA.log_quad_form([1.0, 0.0, 0.0], X), S)
case("sle_sign", () -> SA.sum_log_eigmax(X, 2), S)
case("lor_dist", () -> Manifolds.distance(L, [0.0, 0.0, 1.0], p), L)
case("lor_barrier", () -> SA.lorentz_log_barrier(p), L)
case("lor_homq", () -> SA.lorentz_homogeneous_quadratic([2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0], p), L)
case("lor_homd", () -> SA.lorentz_homogeneous_diagonal([2.0, 2.0, 1.0], p), L)
case("lor_ls", () -> SA.lorentz_least_squares([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [0.0, 0.0, -1.0], p), L)
case("lor_nhq", () -> SA.lorentz_nonhomogeneous_quadratic(A, [0.0, 0.0, 2.0], 1.0, p), L)
case("lor_comp", () -> 2.0 * Manifolds.distance(L, [0.0, 0.0, 1.0], p) + SA.lorentz_log_barrier(p), L)
case("lor_homq_neg", () -> SA.lorentz_homogeneous_quadratic(Matrix(Diagonal([1.0, 1.0, -0.5])), p), L)
case("lor_homd_neg", () -> SA.lorentz_homogeneous_diagonal([1.0, 1.0, -0.5], p), L)

# branches and quirks of find_gcurvature
case("min_logdet_tr", () -> min(logdet(X), tr(X)), S)
case("min_tr_logdet", () -> min(tr(X), logdet(X)), S)
case("negeig_pow", () -> (-eigmax(X))^0.3, S)
case("sqrt_eigmax", () -> sqrt(eigmax(X)), S)
case("sle_exp_inv", () -> SA.sum_log_eigmax(exp, inv(X), 2), S)
case("sle_inv", () -> SA.sum_log_eigmax(inv(X), 2), S)
case("affine_tr_inv", () -> SA.affine_map(tr, inv(X), B), S)
case("affine_diag", () -> SA.affine_map(diag, X, B), S)
case("eigmax_inv", () -> eigmax(inv(X)), S)
case("tr_adjoint", () -> tr(X'), S)
case("logdet_diag", () -> logdet(diag(X)), S)
case("logdet_sum_shift", () -> logdet(X + B), S)
case("logdet_XY", () -> logdet(X + Y), S)
case("tr_XY", () -> tr(X + Y), S)
case("scalar_mat", () -> SA.scalar_mat(X, 3), S)
case("quad_canon", () -> p' * X * p, S)
case("quad_canon_euclid", () -> p' * X * p)
case("conj_euclid", () -> SA.conjugation(X, A))
case("abs_eigmax", () -> abs(eigmax(X)), S)
case("abs_eigmax_euclid", () -> abs(eigmax(X)))
case("neg_tr", () -> -tr(X), S)
case("exp_neg_tr", () -> exp(-tr(X)), S)
case("log_neg_tr", () -> log(-tr(X)), S)
case("XY_mul", () -> tr(X * Y), S)
case("dist_spd_on_lorentz", () -> Manifolds.distance(S, A, X), L)
case("barrier_on_spd", () -> SA.lorentz_log_barrier(p), S)
case("euclid_manifold", () -> tr(X), Manifolds.Euclidean(3))
case("eigsummax_euclid", () -> SA.eigsummax(X, 2))
case("inv_euclid", () -> tr(inv(X)))
