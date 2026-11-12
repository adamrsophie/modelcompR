# The low storage Runge-Kutta methods (R/solvers-lowstorage.R) against
# OrdinaryDiffEqLowStorageRK: tools/lowstorage-probe.jl runs every method with
# a fixed step, a fixed step with saveat and, for the adaptive ones, with
# tolerances, with step rejections and with saveat. Steps, function
# evaluations, accepted and rejected steps, the final state and two of the
# saved times agree to the bit on the reference platform (helper-platform.R;
# elsewhere the final state is compared to a tolerance). The runs of the R+
# methods with rejected steps differ on purpose (see PORTING.md) and are
# checked against a reference solution.

lsrk_f <- function(u, p, t) c(u[2], -u[1] + 0.3 * (1 - t / 3) * u[1]^2)
lsrk_g <- function(u, p, t) c(-20 * (u[1] - t * (2 - t)), 0.5 * u[1] * u[2])
lsrk_P1 <- ODEProblem(lsrk_f, c(1, 0.5), c(0, 3))
lsrk_P2 <- ODEProblem(lsrk_g, c(0, 1), c(0, 2))
lsrk_hex <- function(s) as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])

lsrk_solve <- function(alg, case) {
  switch(case,
    fixed = solve(lsrk_P1, alg, dt = 0.05, adaptive = FALSE),
    fixed_saveat = solve(lsrk_P1, alg, dt = 0.05, adaptive = FALSE, saveat = 0.37),
    adaptive = solve(lsrk_P1, alg, reltol = 1e-6, abstol = 1e-8),
    reject = solve(lsrk_P2, alg, reltol = 1e-4, abstol = 1e-6, dt = 0.5),
    saveat = solve(lsrk_P1, alg, reltol = 1e-6, abstol = 1e-8, saveat = 0.37))
}

lsrk_ref <- matrix(c(
  "AlJahdaliAdv3S113", "fixed", 61, 661, 60, 0, "-0x1.41c1feccf8761p-1,-0x1.c76217723f794p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S113", "fixed_saveat", 10, 661, 60, 0, "-0x1.41c1feccf8761p-1,-0x1.c76217723f794p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S154", "fixed", 61, 901, 60, 0, "-0x1.41c2049d22ec4p-1,-0x1.c762226215d9p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S154", "fixed_saveat", 10, 901, 60, 0, "-0x1.41c2049d22ec4p-1,-0x1.c762226215d9p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S165", "fixed", 61, 961, 60, 0, "-0x1.41c204aa64759p-1,-0x1.c762225947a4ep-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S165", "fixed_saveat", 10, 961, 60, 0, "-0x1.41c204aa64759p-1,-0x1.c762225947a4ep-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S42", "fixed", 61, 241, 60, 0, "-0x1.41dc200d40c3cp-1,-0x1.c75251e36df12p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S42", "fixed_saveat", 10, 241, 60, 0, "-0x1.41dc200d40c3cp-1,-0x1.c75251e36df12p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S53", "fixed", 61, 301, 60, 0, "-0x1.41c1de2718062p-1,-0x1.c761d8a10b8bcp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S53", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1de2718062p-1,-0x1.c761d8a10b8bcp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S64", "fixed", 61, 361, 60, 0, "-0x1.41c204235b20dp-1,-0x1.c76222ae852e1p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S64", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c204235b20dp-1,-0x1.c76222ae852e1p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S82", "fixed", 61, 481, 60, 0, "-0x1.41c887dcd3ca5p-1,-0x1.c75e25531542ep-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S82", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c887dcd3ca5p-1,-0x1.c75e25531542ep-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliAdv3S85", "fixed", 61, 481, 60, 0, "-0x1.41c204a9c0c2fp-1,-0x1.c76222592c669p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliAdv3S85", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204a9c0c2fp-1,-0x1.c76222592c669p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S113", "fixed", 61, 661, 60, 0, "-0x1.41c1fdb16d7bep-1,-0x1.c76215e28929p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S113", "fixed_saveat", 10, 661, 60, 0, "-0x1.41c1fdb16d7bep-1,-0x1.c76215e28929p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S154", "fixed", 61, 901, 60, 0, "-0x1.41c2049d55ea9p-1,-0x1.c76222620b715p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S154", "fixed_saveat", 10, 901, 60, 0, "-0x1.41c2049d55ea9p-1,-0x1.c76222620b715p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S165", "fixed", 61, 961, 60, 0, "-0x1.41c204aa48b04p-1,-0x1.c762225945927p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S165", "fixed_saveat", 10, 961, 60, 0, "-0x1.41c204aa48b04p-1,-0x1.c762225945927p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S42", "fixed", 61, 241, 60, 0, "-0x1.41dbc4b511726p-1,-0x1.c7526b24cee0dp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S42", "fixed_saveat", 10, 241, 60, 0, "-0x1.41dbc4b511726p-1,-0x1.c7526b24cee0dp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S53", "fixed", 61, 301, 60, 0, "-0x1.41c1de407bd25p-1,-0x1.c761d8dcad1ccp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S53", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1de407bd25p-1,-0x1.c761d8dcad1ccp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S64", "fixed", 61, 361, 60, 0, "-0x1.41c2042fd8bbap-1,-0x1.c76222a451962p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S64", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c2042fd8bbap-1,-0x1.c76222a451962p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S82", "fixed", 61, 481, 60, 0, "-0x1.41c892c17a444p-1,-0x1.c75e0e9561352p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S82", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c892c17a444p-1,-0x1.c75e0e9561352p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "AlJahdaliVor3S85", "fixed", 61, 481, 60, 0, "-0x1.41c204aa7a45ap-1,-0x1.c7622259e7b16p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "AlJahdaliVor3S85", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204aa7a45ap-1,-0x1.c7622259e7b16p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CFRLDDRK64", "fixed", 61, 361, 60, 0, "-0x1.41c204a4a82eap-1,-0x1.c762229f1ec71p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CFRLDDRK64", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c204a4a82eap-1,-0x1.c762229f1ec71p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK43_2", "fixed", 61, 241, 60, 0, "-0x1.41c1f717f32b7p-1,-0x1.c76232f206596p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK43_2", "fixed_saveat", 10, 241, 60, 0, "-0x1.41c1f717f32b7p-1,-0x1.c76232f206596p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK43_2", "adaptive", 175, 699, 174, 0, "-0x1.41c203f3348e5p-1,-0x1.c762239a7fdb6p-1", "87,174", "0x1.9bac6d4abe1b3p+0,0x1.7dfb85bc6c0acp+1",
  "CKLLSRK43_2", "saveat", 10, 699, 174, 0, "-0x1.41c203f3348e5p-1,-0x1.c762239a7fdb6p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3C", "fixed", 61, 301, 60, 0, "-0x1.41fdbbd2cbe86p-1,-0x1.c743f05244ab1p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3C", "fixed_saveat", 10, 301, 60, 0, "-0x1.41fdbbd2cbe86p-1,-0x1.c743f05244ab1p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3C", "adaptive", 250, 1248, 249, 0, "-0x1.41ce78c7909fcp-1,-0x1.c760bed19d701p-1", "125,249", "0x1.7f3a48cd93c87p-1,0x1.7ddbee236d81cp+1",
  "CKLLSRK54_3C", "saveat", 10, 1248, 249, 0, "-0x1.41ce78c7909fcp-1,-0x1.c760bed19d701p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3C_3R", "fixed", 61, 301, 60, 0, "-0x1.41c20369b2f05p-1,-0x1.c7622315a5a78p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3C_3R", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c20369b2f05p-1,-0x1.c7622315a5a78p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3C_3R", "adaptive", 40, 198, 39, 0, "-0x1.41c1fa641a936p-1,-0x1.c76227eb5b2f8p-1", "20,39", "0x1.56ee028a45a02p+0,0x1.79babf5d11438p+1",
  "CKLLSRK54_3C_3R", "saveat", 10, 198, 39, 0, "-0x1.41c1fa641a936p-1,-0x1.c76227eb5b2f8p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3M_3R", "fixed", 61, 301, 60, 0, "-0x1.41c204acabbeep-1,-0x1.c762225278e3ep-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3M_3R", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c204acabbeep-1,-0x1.c762225278e3ep-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3M_3R", "adaptive", 38, 188, 37, 0, "-0x1.41c204eec23afp-1,-0x1.c762225815ae2p-1", "19,37", "0x1.49210ed53b9bap+0,0x1.7f82d6d36384bp+1",
  "CKLLSRK54_3M_3R", "saveat", 10, 188, 37, 0, "-0x1.41c204eec23afp-1,-0x1.c762225815ae2p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3M_4R", "fixed", 61, 301, 60, 0, "-0x1.41c204ae9cbeep-1,-0x1.c762225fa4f8p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3M_4R", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c204ae9cbeep-1,-0x1.c762225fa4f8p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3M_4R", "adaptive", 35, 173, 34, 0, "-0x1.41c205559a9d7p-1,-0x1.c76223180b9b8p-1", "17,34", "0x1.2ba25f71b2f13p+0,0x1.794ec2beacaabp+1",
  "CKLLSRK54_3M_4R", "saveat", 10, 173, 34, 0, "-0x1.41c205559a9d7p-1,-0x1.c76223180b9b8p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3N_3R", "fixed", 61, 301, 60, 0, "-0x1.41c2036b68d08p-1,-0x1.c76223172c76dp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3N_3R", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c2036b68d08p-1,-0x1.c76223172c76dp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3N_3R", "adaptive", 40, 198, 39, 0, "-0x1.41c1f8f8c52ffp-1,-0x1.c762290c4a25cp-1", "20,39", "0x1.49ce8449c2c75p+0,0x1.7d87314752c6ep+1",
  "CKLLSRK54_3N_3R", "saveat", 10, 198, 39, 0, "-0x1.41c1f8f8c52ffp-1,-0x1.c762290c4a25cp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3N_4R", "fixed", 61, 301, 60, 0, "-0x1.41c2033f82973p-1,-0x1.c76223251ce63p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK54_3N_4R", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c2033f82973p-1,-0x1.c76223251ce63p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK54_3N_4R", "adaptive", 42, 208, 41, 0, "-0x1.41c1fa7c102e6p-1,-0x1.c76227dce5845p-1", "21,41", "0x1.4aec32f3a5592p+0,0x1.7d15dac430967p+1",
  "CKLLSRK54_3N_4R", "saveat", 10, 208, 41, 0, "-0x1.41c1fa7c102e6p-1,-0x1.c76227dce5845p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK65_4M_4R", "fixed", 61, 361, 60, 0, "-0x1.41c204aab3ce7p-1,-0x1.c762225a07a21p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK65_4M_4R", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c204aab3ce7p-1,-0x1.c762225a07a21p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK65_4M_4R", "adaptive", 22, 129, 21, 0, "-0x1.41c2065e490a3p-1,-0x1.c76224959516bp-1", "11,21", "0x1.024b1d318e785p+0,0x1.7277223173ab8p+1",
  "CKLLSRK65_4M_4R", "saveat", 10, 129, 21, 0, "-0x1.41c2065e490a3p-1,-0x1.c76224959516bp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK75_4M_5R", "fixed", 61, 421, 60, 0, "-0x1.41c204aa44dabp-1,-0x1.c762225945893p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK75_4M_5R", "fixed_saveat", 10, 421, 60, 0, "-0x1.41c204aa44dabp-1,-0x1.c762225945893p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK75_4M_5R", "adaptive", 19, 129, 18, 0, "-0x1.41c204be2efe7p-1,-0x1.c7622231f36ebp-1", "9,18", "0x1.e892bc4b27e92p-1,0x1.6ed26f864870bp+1",
  "CKLLSRK75_4M_5R", "saveat", 10, 129, 18, 0, "-0x1.41c204be2efe7p-1,-0x1.c7622231f36ebp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4C_3R", "fixed", 61, 481, 60, 0, "-0x1.41c204a985c13p-1,-0x1.c76222589bb9cp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK85_4C_3R", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204a985c13p-1,-0x1.c76222589bb9cp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4C_3R", "adaptive", 20, 155, 19, 0, "-0x1.41c202f19830ep-1,-0x1.c7622040564b1p-1", "10,19", "0x1.45ef85eb5ac8cp+0,0x1.71692863c4916p+1",
  "CKLLSRK85_4C_3R", "saveat", 10, 155, 19, 0, "-0x1.41c202f19830ep-1,-0x1.c7622040564b1p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4FM_4R", "fixed", 61, 481, 60, 0, "-0x1.41c204aa43543p-1,-0x1.c762225947336p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK85_4FM_4R", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204aa43543p-1,-0x1.c762225947336p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4M_3R", "fixed", 61, 481, 60, 0, "-0x1.41c204aadbf0fp-1,-0x1.c7622259fa052p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK85_4M_3R", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204aadbf0fp-1,-0x1.c7622259fa052p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4M_3R", "adaptive", 22, 171, 21, 0, "-0x1.41c20693cd6efp-1,-0x1.c76225a0149d3p-1", "11,21", "0x1.f8d4fef975c94p-1,0x1.7a3200df345bdp+1",
  "CKLLSRK85_4M_3R", "saveat", 10, 171, 21, 0, "-0x1.41c20693cd6efp-1,-0x1.c76225a0149d3p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4P_3R", "fixed", 61, 481, 60, 0, "-0x1.41c204a9d98d5p-1,-0x1.c76222597f87p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK85_4P_3R", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c204a9d98d5p-1,-0x1.c76222597f87p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK85_4P_3R", "adaptive", 22, 171, 21, 0, "-0x1.41c203bef287cp-1,-0x1.c7622251a11fp-1", "11,21", "0x1.55c31f289e766p+0,0x1.7fdf67adb6cd9p+1",
  "CKLLSRK85_4P_3R", "saveat", 10, 171, 21, 0, "-0x1.41c203bef287cp-1,-0x1.c7622251a11fp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK95_4C", "fixed", 61, 541, 60, 0, "-0x1.4248c7d885807p-1,-0x1.c86da81eb28a9p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK95_4C", "fixed_saveat", 10, 541, 60, 0, "-0x1.4248c7d885807p-1,-0x1.c86da81eb28a9p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK95_4M", "fixed", 61, 541, 60, 0, "-0x1.41c204aa43be6p-1,-0x1.c762225a4f257p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK95_4M", "fixed_saveat", 10, 541, 60, 0, "-0x1.41c204aa43be6p-1,-0x1.c762225a4f257p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK95_4M", "adaptive", 22, 192, 21, 0, "-0x1.41c2044bfb177p-1,-0x1.c76223ed1a3cp-1", "11,21", "0x1.473687f6ab821p+0,0x1.790c67af4fcc8p+1",
  "CKLLSRK95_4M", "saveat", 10, 192, 21, 0, "-0x1.41c2044bfb177p-1,-0x1.c76223ed1a3cp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK95_4S", "fixed", 61, 541, 60, 0, "-0x1.41c204a962129p-1,-0x1.c76222591883p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CKLLSRK95_4S", "fixed_saveat", 10, 541, 60, 0, "-0x1.41c204a962129p-1,-0x1.c76222591883p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CKLLSRK95_4S", "adaptive", 20, 174, 19, 0, "-0x1.41c202310d5e3p-1,-0x1.c762214e0bdc2p-1", "10,19", "0x1.345639ad5642ap+0,0x1.6abc7f3286a1ep+1",
  "CKLLSRK95_4S", "saveat", 10, 174, 19, 0, "-0x1.41c202310d5e3p-1,-0x1.c762214e0bdc2p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "CarpenterKennedy2N54", "fixed", 61, 301, 60, 0, "-0x1.41c2034cac10ep-1,-0x1.c76222d2df6eap-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "CarpenterKennedy2N54", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c2034cac10ep-1,-0x1.c76222d2df6eap-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "DGLDDRK73_C", "fixed", 61, 421, 60, 0, "-0x1.41c1edd97c444p-1,-0x1.c761fb10de447p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "DGLDDRK73_C", "fixed_saveat", 10, 421, 60, 0, "-0x1.41c1edd97c444p-1,-0x1.c761fb10de447p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "DGLDDRK84_C", "fixed", 61, 481, 60, 0, "-0x1.41c2044a217e2p-1,-0x1.c76222658cd49p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "DGLDDRK84_C", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c2044a217e2p-1,-0x1.c76222658cd49p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "DGLDDRK84_F", "fixed", 61, 481, 60, 0, "-0x1.41c20488739e3p-1,-0x1.c762225b5ca7dp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "DGLDDRK84_F", "fixed_saveat", 10, 481, 60, 0, "-0x1.41c20488739e3p-1,-0x1.c762225b5ca7dp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "NDBLSRK124", "fixed", 61, 721, 60, 0, "-0x1.41c2045c05057p-1,-0x1.c762225c517cp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "NDBLSRK124", "fixed_saveat", 10, 721, 60, 0, "-0x1.41c2045c05057p-1,-0x1.c762225c517cp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "NDBLSRK134", "fixed", 61, 781, 60, 0, "-0x1.41c2049412573p-1,-0x1.c76222581805ap-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "NDBLSRK134", "fixed_saveat", 10, 781, 60, 0, "-0x1.41c2049412573p-1,-0x1.c76222581805ap-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "NDBLSRK144", "fixed", 61, 841, 60, 0, "-0x1.41c20490b30efp-1,-0x1.c762225e0fbfbp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "NDBLSRK144", "fixed_saveat", 10, 841, 60, 0, "-0x1.41c20490b30efp-1,-0x1.c762225e0fbfbp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ORK256", "fixed", 61, 301, 60, 0, "-0x1.41c2f3f4a6be4p-1,-0x1.c7623b4004e78p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ORK256", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c2f3f4a6be4p-1,-0x1.c7623b4004e78p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S105", "fixed", 61, 601, 60, 0, "-0x1.41c204aa47235p-1,-0x1.c7622259504c1p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S105", "fixed_saveat", 10, 601, 60, 0, "-0x1.41c204aa47235p-1,-0x1.c7622259504c1p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S173", "fixed", 61, 1021, 60, 0, "-0x1.41c20184b807p-1,-0x1.c7621c8745a22p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S173", "fixed_saveat", 10, 1021, 60, 0, "-0x1.41c20184b807p-1,-0x1.c7621c8745a22p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S184", "fixed", 61, 1081, 60, 0, "-0x1.41c204a13b06fp-1,-0x1.c762225f0d546p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S184", "fixed_saveat", 10, 1081, 60, 0, "-0x1.41c204a13b06fp-1,-0x1.c762225f0d546p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S205", "fixed", 61, 1201, 60, 0, "-0x1.41c204aa46eb9p-1,-0x1.c76222594abcap-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S205", "fixed_saveat", 10, 1201, 60, 0, "-0x1.41c204aa46eb9p-1,-0x1.c76222594abcap-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S32", "fixed", 61, 181, 60, 0, "-0x1.41f1c2c72ca22p-1,-0x1.c74560f8281c3p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S32", "fixed_saveat", 10, 181, 60, 0, "-0x1.41f1c2c72ca22p-1,-0x1.c74560f8281c3p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S53", "fixed", 61, 301, 60, 0, "-0x1.41c1dd082280bp-1,-0x1.c761d71e8b807p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S53", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1dd082280bp-1,-0x1.c761d71e8b807p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S82", "fixed", 61, 481, 60, 0, "-0x1.41ca24f5188c9p-1,-0x1.c75e82e75f672p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S82", "fixed_saveat", 10, 481, 60, 0, "-0x1.41ca24f5188c9p-1,-0x1.c75e82e75f672p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "ParsaniKetchesonDeconinck3S94", "fixed", 61, 541, 60, 0, "-0x1.41c20480d16d9p-1,-0x1.c76222756986ap-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "ParsaniKetchesonDeconinck3S94", "fixed_saveat", 10, 541, 60, 0, "-0x1.41c20480d16d9p-1,-0x1.c76222756986ap-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp35", "fixed", 61, 301, 60, 0, "-0x1.41c1dd0823f4fp-1,-0x1.c761d71e8c021p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3Sp35", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1dd0823f4fp-1,-0x1.c761d71e8c021p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp35", "adaptive", 56, 278, 55, 0, "-0x1.41c1be9197f69p-1,-0x1.c76194719dbf1p-1", "28,55", "0x1.685a624085188p+0,0x1.78743ae42cc67p+1",
  "RDPK3Sp35", "reject", 28, 166, 27, 6, "0x1.874f12344070dp-4,0x1.f16fc982adcf7p+0", "14,27", "0x1.ba7627fee2d5p-2,0x1.f827a317c971ep+0",
  "RDPK3Sp35", "saveat", 10, 278, 55, 0, "-0x1.41c1be9197f69p-1,-0x1.c76194719dbf1p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp49", "fixed", 61, 541, 60, 0, "-0x1.41c20480d16fp-1,-0x1.c762227569a29p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3Sp49", "fixed_saveat", 10, 541, 60, 0, "-0x1.41c20480d16fp-1,-0x1.c762227569a29p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp49", "adaptive", 32, 282, 31, 0, "-0x1.41c1fb44964a3p-1,-0x1.c7622918ddbadp-1", "16,31", "0x1.e3f5cd4dde28ep-1,0x1.7aa86bc6d6d4cp+1",
  "RDPK3Sp49", "reject", 21, 235, 20, 6, "0x1.84f7774d1e6aap-4,0x1.f16fc81757afep+0", "10,20", "0x1.74febb2b99e99p-1,0x1.ead278f40ffe3p+0",
  "RDPK3Sp49", "saveat", 10, 282, 31, 0, "-0x1.41c1fb44964a3p-1,-0x1.c7622918ddbadp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp510", "fixed", 61, 601, 60, 0, "-0x1.41c204aa4719ep-1,-0x1.c76222595039bp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3Sp510", "fixed_saveat", 10, 601, 60, 0, "-0x1.41c204aa4719ep-1,-0x1.c76222595039bp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3Sp510", "adaptive", 13, 123, 12, 0, "-0x1.41c207460a77fp-1,-0x1.c762265caa58ep-1", "6,12", "0x1.f599eaf1dee04p-2,0x1.53f3efd8ef46bp+1",
  "RDPK3Sp510", "reject", 14, 161, 13, 3, "0x1.84fe4a5e52ab2p-4,0x1.f17046915ee1ep+0", "7,13", "0x1.d85ee6757eb33p-1,0x1.ff98a81474f66p+0",
  "RDPK3Sp510", "saveat", 10, 123, 12, 0, "-0x1.41c207460a77fp-1,-0x1.c762265caa58ep-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL35", "fixed", 61, 301, 60, 0, "-0x1.41c1dd081be45p-1,-0x1.c761d71e954ffp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3SpFSAL35", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1dd081be45p-1,-0x1.c761d71e954ffp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL35", "adaptive", 50, 248, 49, 0, "-0x1.41c171b7f50efp-1,-0x1.c760cebb064d5p-1", "25,49", "0x1.71cd28c550b1cp+0,0x1.79681c0a7b21dp+1",
  "RDPK3SpFSAL35", "reject", 42, 231, 41, 5, "0x1.855488e20285fp-4,0x1.f16fefaab9738p+0", "21,41", "0x1.82d75888ddb31p-1,0x1.fd5e98e5f87dp+0",
  "RDPK3SpFSAL35", "saveat", 10, 248, 49, 0, "-0x1.41c171b7f50efp-1,-0x1.c760cebb064d5p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL49", "fixed", 61, 541, 60, 0, "-0x1.41c20480d16c4p-1,-0x1.c7622275699fcp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3SpFSAL49", "fixed_saveat", 10, 541, 60, 0, "-0x1.41c20480d16c4p-1,-0x1.c7622275699fcp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL49", "adaptive", 24, 210, 23, 0, "-0x1.41c1e0a369b7cp-1,-0x1.c7623c07d740bp-1", "12,23", "0x1.b5fabe161dca4p-1,0x1.7d98d39d726fdp+1",
  "RDPK3SpFSAL49", "reject", 19, 208, 18, 5, "0x1.8434a3cf4202cp-4,0x1.f16fa59cf244p+0", "9,18", "0x1.271e73f62e13dp-1,0x1.f1eec0e549a0bp+0",
  "RDPK3SpFSAL49", "saveat", 10, 210, 23, 0, "-0x1.41c1e0a369b7cp-1,-0x1.c7623c07d740bp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL510", "fixed", 61, 601, 60, 0, "-0x1.41c204aa471c6p-1,-0x1.c762225950404p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RDPK3SpFSAL510", "fixed_saveat", 10, 601, 60, 0, "-0x1.41c204aa471c6p-1,-0x1.c762225950404p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RDPK3SpFSAL510", "adaptive", 14, 133, 13, 0, "-0x1.41c20684ea3c8p-1,-0x1.c762250f84b56p-1", "7,13", "0x1.6a5732528b2c7p-1,0x1.67ff554d9f30ap+1",
  "RDPK3SpFSAL510", "reject", 15, 181, 14, 4, "0x1.850a425d75ee8p-4,0x1.f1702eb53f5c1p+0", "7,14", "0x1.4fdf24c8a42aep-1,0x1.d6906521c8d2ep+0",
  "RDPK3SpFSAL510", "saveat", 10, 133, 13, 0, "-0x1.41c20684ea3c8p-1,-0x1.c762250f84b56p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "RK46NL", "fixed", 61, 361, 60, 0, "-0x1.41c204953a555p-1,-0x1.c762227a01453p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "RK46NL", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c204953a555p-1,-0x1.c762227a01453p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "SHLDDRK52", "fixed", 61, 301, 60, 0, "-0x1.41c1b122c33cbp-1,-0x1.c761fe0ee2cd5p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "SHLDDRK52", "fixed_saveat", 10, 301, 60, 0, "-0x1.41c1b122c33cbp-1,-0x1.c761fe0ee2cd5p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "SHLDDRK64", "fixed", 61, 361, 60, 0, "-0x1.41c1f72a90587p-1,-0x1.c7622cb1b44e5p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "SHLDDRK64", "fixed_saveat", 10, 361, 60, 0, "-0x1.41c1f72a90587p-1,-0x1.c7622cb1b44e5p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "SHLDDRK_2N", "fixed", 61, 331, 60, 0, "-0x1.41c1f9dca29c8p-1,-0x1.c7622ab180edcp-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "SHLDDRK_2N", "fixed_saveat", 10, 331, 60, 0, "-0x1.41c1f9dca29c8p-1,-0x1.c7622ab180edcp-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1",
  "TSLDDRK74", "fixed", 61, 421, 60, 0, "-0x1.41c204a3a4b3cp-1,-0x1.c762227f510c4p-1", "30,60", "0x1.7333333333336p+0,0x1.7999999999994p+1",
  "TSLDDRK74", "fixed_saveat", 10, 421, 60, 0, "-0x1.41c204a3a4b3cp-1,-0x1.c762227f510c4p-1", "5,9", "0x1.7ae147ae147aep+0,0x1.7ae147ae147aep+1"
), ncol = 9, byrow = TRUE)

test_that("the low storage methods reproduce OrdinaryDiffEqLowStorageRK", {
  for (i in seq_len(nrow(lsrk_ref))) {
    r <- lsrk_ref[i, ]
    lab <- paste(r[1], r[2])
    s <- lsrk_solve(get(r[1])(), r[2])
    if (!ref_platform()) {
      expect_equal(s$u[[length(s$u)]], lsrk_hex(r[7]), tolerance = 1e-6, label = paste(lab, "u"))
      next
    }
    expect_identical(length(s$t), as.integer(r[3]), label = paste(lab, "steps"))
    expect_identical(s$stats$nf, as.integer(r[4]), label = paste(lab, "nf"))
    expect_identical(s$stats$naccept, as.integer(r[5]), label = paste(lab, "naccept"))
    expect_identical(s$stats$nreject, as.integer(r[6]), label = paste(lab, "nreject"))
    expect_identical(s$u[[length(s$u)]], lsrk_hex(r[7]), label = paste(lab, "u"))
    expect_identical(s$t[as.integer(strsplit(r[8], ",")[[1]])], lsrk_hex(r[9]), label = paste(lab, "t"))
  }
})

test_that("the R+ methods start a retried step from the accepted state", {
  ref <- solve(lsrk_P2, Vern9(), reltol = 1e-12, abstol = 1e-14)
  uref <- ref$u[[length(ref$u)]]
  for (a in c("CKLLSRK43_2", "CKLLSRK54_3C", "CKLLSRK95_4C", "CKLLSRK95_4M", "CKLLSRK95_4S", "CKLLSRK54_3C_3R", "CKLLSRK54_3M_3R", "CKLLSRK54_3N_3R", "CKLLSRK85_4C_3R", "CKLLSRK85_4M_3R", "CKLLSRK85_4P_3R", "CKLLSRK54_3M_4R", "CKLLSRK54_3N_4R", "CKLLSRK65_4M_4R", "CKLLSRK85_4FM_4R", "CKLLSRK75_4M_5R")) {
    s <- lsrk_solve(get(a)(), "reject")
    expect_true(s$stats$nreject > 0L, label = a)
    expect_equal(s$u[[length(s$u)]], uref, tolerance = 1e-3, label = a)
  }
})

test_that("the low storage methods interpolate and restart after events", {
  s <- solve(lsrk_P1, RDPK3Sp35(), reltol = 1e-6, abstol = 1e-8)
  expect_equal(s(1.234), solve(lsrk_P1, Vern9(), reltol = 1e-12, abstol = 1e-14)(1.234), tolerance = 1e-5)
  cb <- DiscreteCallback(function(u, t, integ) t == 1, function(integ) integ$u <- integ$u * 0.5)
  for (a in c("SHLDDRK_2N", "CarpenterKennedy2N54", "RDPK3SpFSAL49")) {
    s <- solve(lsrk_P1, get(a)(), dt = 0.05, adaptive = FALSE, callback = cb, tstops = 1)
    expect_equal(s(3), solve(lsrk_P1, Vern9(), callback = cb, tstops = 1, reltol = 1e-12, abstol = 1e-14)(3),
                 tolerance = 1e-5, label = a)
  }
})
