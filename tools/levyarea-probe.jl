# Reference values of StochasticDiffEqLevyArea and the iterated integral
# evaluators of StochasticDiffEqCore for tests/testthat/test-sde-levyarea.R.
# Run: julia --project=tools/julia --startup-file=no tools/levyarea-probe.jl
# The random numbers are fixed: z(k) = mod(37k, 101) / 50 - 1, filled column
# by column into X, Y and the tail, as generate_coefficients draws them.
using StochasticDiffEq, StochasticDiffEqLevyArea, Random
const trigamma = StochasticDiffEqLevyArea.trigamma
const LA = StochasticDiffEqLevyArea

z(k) = mod(37 * k, 101) / 50 - 1
function coeffs(m, n, alg)
    k = 0
    nxt(len) = [z(k += 1) for _ in 1:len]
    if alg isa MronRoe
        X = reshape(nxt(m * n), m, n); Y = reshape(nxt(n * m), n, m)
        tail = vcat(nxt(m), nxt((m^2 - m) ÷ 2))
    else
        X = reshape(nxt(n * m), n, m); Y = reshape(nxt(m * n), m, n)
        tail = alg isa Fourier ? Float64[] : alg isa Milstein ? nxt(m) : nxt((m^2 - m) ÷ 2)
    end
    return LevyAreaCoefficients{Float64}(X, Y, tail, m, n)
end

W = [0.3, -0.7, 1.1]
m, n = 3, 5
println("trigamma ", repr([trigamma(k + 1) for k in 1:8]))
for alg in (Fourier(), Milstein(), Wiktorsson(), MronRoe())
    c = coeffs(m, n, alg)
    println(nameof(typeof(alg)), " area ", repr(vec(levyarea(W, n, alg, c))))
    println(nameof(typeof(alg)), " length ", coefficient_length(m, n, alg))
    # the RNG form draws as generate_coefficients
    Random.seed!(3); A = levyarea(W, n, alg)
    Random.seed!(3); B = levyarea(W, n, alg, generate_coefficients(m, n, alg))
    println(nameof(typeof(alg)), " rng_diff ", maximum(abs.(A - B)))
end
h = 0.04
dW = W * sqrt(h)
cM = coeffs(m, n, MronRoe())
println("ii_mronroe_ito ", repr(vec(iterated_integrals(dW, h, cM))))
println("ii_mronroe_strat ", repr(vec(iterated_integrals(dW, h, cM; ito_correction = false))))
cF = coeffs(m, n, Fourier())
println("ii_fourier_ito ", repr(vec(iterated_integrals(dW, h, cF; alg = Fourier()))))
println("ii_scalar ", repr([iterated_integrals(0.3, 0.04), iterated_integrals(0.3, 0.04; ito_correction = false),
                            iterated_integrals(0.3, 0.5, 0.04, 0.01)]))
for (dim, hh) in ((2, 0.1), (3, 0.01), (5, 0.001), (10, 0.5), (4, 1.0e-4))
    eps = hh^(3 / 2)
    tn = [terms_needed(dim, hh, eps, a, MaxL2()) for a in LA.ITER_INT_ALGS]
    tf = [terms_needed(dim, hh, eps, a, FrobeniusL2()) for a in LA.ITER_INT_ALGS]
    println("terms ", dim, " ", hh, " ", repr(tn), " ", repr(tf), " ", nameof(typeof(optimal_algorithm(dim, hh))))
end
q = [1.0, 0.5, 0.25]
println("terms_q ", repr([terms_needed(3, q, 0.01, 0.001, a, FrobeniusL2()) for a in LA.ITER_INT_ALGS]),
        " ", repr([terms_needed(3, q, 0.01, 0.001, a, MaxL2()) for a in LA.ITER_INT_ALGS]),
        " ", nameof(typeof(optimal_algorithm(3, q, 0.01, 0.001))))
# the Q-Wiener iterated_integrals passes the algorithm as a keyword levyarea
# does not take; the intended value, assembled by hand from coefficients
qerr = try
    iterated_integrals(dW, q, h, 0.01; alg = MronRoe()); "none"
catch e
    string(typeof(e))
end
println("qwiener_error ", qerr)
nq = terms_needed(3, q, h, 0.01, MronRoe(), FrobeniusL2())
cq = coeffs(m, nq, MronRoe())
Iq = levyarea(dW ./ q ./ √h, nq, MronRoe(), cq)
LA.ito_correction!(Iq)
Iq .= 0.5 .* dW .* dW' .+ h .* q' .* Iq .* q
println("qwiener n ", nq, " ", repr(vec(Iq)))
pts = [0.0, 0.01, 0.025, 0.04]
println("path_mronroe ", repr(reduce(vcat, reconstruct_path(dW, h, cM, pts))))
println("path_fourier ", repr(reduce(vcat, reconstruct_path(dW, h, cF, pts))))
println("sub_ito ", repr(vec(iterated_integrals_subinterval(dW, h, cM, 0.01, 0.03))))
println("sub_strat ", repr(vec(iterated_integrals_subinterval(dW, h, cM, 0.0, 0.04; n_quadrature = 16, ito_correction = false))))
# iterated integral evaluators
d2 = [0.1, -0.2]
println("jdiag ", repr(get_iterated_I(0.01, d2, nothing, JDiagonal_oop())))
println("jcomm ", repr(vec(get_iterated_I(0.01, d2, nothing, JCommute_oop()))))
for alg in (Fourier(), MronRoe())
    p = terms_needed(2, 0.01, 1 * 0.01^(1 + 1 // 2), alg, MaxL2())
    c2 = coeffs(2, p, alg)
    Jc = 1 // 2 * d2 .* d2' .+ 0.01 .* levyarea(d2 / √0.01, p, alg, c2)
    Random.seed!(5); J1 = get_iterated_I(0.01, d2, nothing, alg)
    Random.seed!(5); c3 = generate_coefficients(2, p, alg)
    J2 = 1 // 2 * d2 .* d2' .+ 0.01 .* levyarea(d2 / √0.01, p, alg, c3)
    println("jlevy ", nameof(typeof(alg)), " p ", p, " ", repr(vec(Jc)), " rng_diff ", maximum(abs.(J1 - J2)))
end
