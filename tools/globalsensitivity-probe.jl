# Reference values for tests/testthat/test-global-sensitivity.R: the point
# sets of QuasiMonteCarlo.jl and the analyses of GlobalSensitivity.jl on given
# samples, written as raw Float64 to
# tests/testthat/fixtures/julia/globalsensitivity.bin (the test reads them in
# the same order).
#
#   julia --project=tools/julia tools/globalsensitivity-probe.jl
using GlobalSensitivity, QuasiMonteCarlo, Random, Statistics
import KernelDensity, ComplexityMeasures
const GS = GlobalSensitivity
io = open(joinpath(@__DIR__, "..", "tests", "testthat", "fixtures", "julia", "globalsensitivity.bin"), "w")
W(x) = write(io, Float64.(vec(collect(x))))
# QuasiMonteCarlo
W(QuasiMonteCarlo.sample(37, 5, SobolSample()))
W(QuasiMonteCarlo.sample(50, 4, HaltonSample()))
W(QuasiMonteCarlo.sample(50, 5, FaureSample()))
W(QuasiMonteCarlo.sample(18, 3, FaureSample()))
W(QuasiMonteCarlo.sample(40, 6, KroneckerSample()))
W(QuasiMonteCarlo.sample(33, 7, LatticeRuleSample()))
W(QuasiMonteCarlo.sample(5, 300, LatticeRuleSample()))
W(QuasiMonteCarlo.sample(30, 1, VanDerCorputSample(base=3)))
W(QuasiMonteCarlo.sample(125, 1, VanDerCorputSample(base=5)))
W(QuasiMonteCarlo.sample(10, [0.0, -1.0, 2.0], [1.0, 3.0, 2.5], HaltonSample()))
W(QuasiMonteCarlo.unif2bits(permutedims(QuasiMonteCarlo.sample(27, 3, FaureSample())), 3; pad = 20))
W(reduce(hcat, QuasiMonteCarlo.generate_design_matrices(16, [0.0,0.0], [1.0,2.0], SobolSample(), 3)))
# Julia ranges
for (a,b,n) in [(0.0, 10.0, 7), (-pi, pi, 1000), (0.1, 0.7, 11), (-1.2345, 3.3, 2048), (2.5e-3, 17.25, 33), (-0.37, 1e5 + 0.1, 101)]
  r = range(a, b, length=n); W(r); W(Float64(step(r)))
end
ishi(X) = sin(X[1]) + 7 * sin(X[2])^2 + 0.1 * X[3]^4 * sin(X[1])
lb = -ones(3) * π; ub = ones(3) * π
A, B = QuasiMonteCarlo.generate_design_matrices(512, lb, ub, SobolSample(), 2)
for est in (:Jansen1999, :Homma1996, :Sobol2007, :Janon2014)
  r = gsa(ishi, Sobol(order=[0,1,2]), A, B; Ei_estimator=est)
  W(r.S1); W(r.ST); W(r.S2)
end
r = gsa(ishi, Sobol(order=[0,1,2], nboot=4), A, B)
W(r.S1); W(r.ST); W(r.S2); W(r.S1_Conf_Int); W(r.ST_Conf_Int); W(r.S2_Conf_Int)
mo(X) = [ishi(X), X[1]*X[2] + X[3]]
r = gsa(mo, Sobol(order=[0,1,2]), A, B)
W(r.S1); W(r.ST); W(r.S2)
# regression / EASI / delta / mi / rsa on Julia samples
X = QuasiMonteCarlo.sample(600, lb, ub, SobolSample())
Y = [ishi(X[:, j]) for j in 1:size(X, 2)]
W(X); W(Y)
r = gsa(X, reshape(Y, 1, :), RegressionGSA(rank=true))
W(r.pearson); W(r.standard_regression); W(r.partial_correlation)
r = gsa(X, Y, EASI()); W(r.S1); W(r.S1_Corr)
r = gsa(X, Y, EASI(dct_method=true)); W(r.S1); W(r.S1_Corr)
grid = range(minimum(Y), maximum(Y), length=2048)
W(KernelDensity.pdf(KernelDensity.kde(Y), grid)); W(KernelDensity.pdf(KernelDensity.kde(Y[1:137]), grid))
M = 20; cuts = range(0, 600, length=M+1)
W([GS._calc_delta(X[i, :], Y, grid, cuts) for i in 1:3])
nb = Int(round(sqrt(600)))
est = ComplexityMeasures.ValueHistogram(nb)
W([ComplexityMeasures.entropy(est, X[i, :]) for i in 1:3]); W(ComplexityMeasures.entropy(est, Y))
W([ComplexityMeasures.entropy(est, ComplexityMeasures.StateSpaceSet(X[i, :], Y)) for i in 1:3])
# RSA analysis with 2 dummies appended as extra rows
Xr = vcat(X, QuasiMonteCarlo.sample(600, zeros(2), ones(2), SobolSample()))
r = GS._compute_rsa(Xr, Y, RSA(n_dummy_parameters=2)); W(r.S); W([r.Sd[1]])
W(Xr[4:5, :])
# RBD-FAST analysis on given Y and s
s0 = range(-π, stop = π, length = 301)
s = [s0[randperm(Xoshiro(i), 301)] for i in 1:3]
x = hcat([0.5 .+ asin.(sin.(s[i])) ./ pi for i in 1:3]...)
Yb = [ishi(-π .+ 2π .* x[j, :]) for j in 1:301]
for si in s; W(si); end; W(Yb)
W(GS._rbdfast_analysis(Yb, s, RBDFAST()))
# eFAST with injected uniforms
rng = Xoshiro(7); us = rand(Xoshiro(7), 3)
r = gsa(ishi, eFAST(), [[lb[i], ub[i]] for i in 1:3]; samples=400, rng=rng)
W(us); W(r.S1); W(r.ST)
# fractional factorial
ff(x) = 2x[1] - 3x[2] + 0.5x[3]*x[1] + x[4]
me, me2 = gsa(ff, FractionalFactorial(); num_params=3, p_range=[[0,1],[-1,2],[1,3]])
W(me); W(me2)
close(io)
