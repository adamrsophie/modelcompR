# Reference values from OrdinaryDiffEqAdamsBashforthMoulton 2.2.1 (Julia),
# printed by tools/abm-probe.jl: saved points, nf, naccept, nreject and the
# final state on a Lotka-Volterra problem with a time dependent term. The
# fixed step methods take dt = 0.05; the adaptive ones reltol = 1e-6 and
# abstol = 1e-8 ("plain"), 1e-4 ("mid") or 1e-2 ("loose"). "tstop" adds a
# tstop at 0.73, "cb" a tstop at 1 where a DiscreteCallback halves u. The
# VCAB*, VCABM* values are those of the Julia methods with the corrections
# of R/solvers-abm.R (the "fix" runs of the probe).

abm_setup <- function() {
  lv <- function(u, p, t) c(p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2] + 0.1 * t)
  list(prob = ODEProblem(lv, c(1, 1), c(0, 2), c(1.5, 3)),
       cb = DiscreteCallback(function(u, t, integ) t == 1, function(integ) integ$u <- integ$u * 0.5))
}

abm_kwargs <- function(case, adaptive, cb) {
  base <- if (!adaptive) list(dt = 0.05, adaptive = FALSE)
    else switch(case, mid = list(reltol = 1e-4, abstol = 1e-4), loose = list(reltol = 1e-2, abstol = 1e-2),
                list(reltol = 1e-6, abstol = 1e-8))
  switch(case, tstop = c(base, list(tstops = 0.73)), cb = c(base, list(tstops = 1, callback = cb)), base)
}

abm_ref <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
alg   case  n   nf  na  nr u1                 u2
AB3   plain 41  43  40  0  5.9145026985122895 2.184419615386514
AB3   tstop 42  44  41  0  5.905715860887023  2.1861520150252933
AB3   cb    42  46  40  0  5.180690627635668  0.33251048777926145
AB4   plain 41  50  40  0  5.9057395752601405 2.181112638483607
AB4   tstop 42  51  41  0  5.897225474366834  2.1827281244709833
AB4   cb    42  60  40  0  5.179243623713821  0.3333358689630135
AB5   plain 41  53  40  0  5.906946716222575  2.179234284249052
AB5   tstop 42  54  41  0  5.898465120895818  2.1808575777526658
AB5   cb    42  66  40  0  5.178792492739055  0.3334507429867335
ABM32 plain 41  120 40  0  5.907607668487664  2.1782545643195492
ABM32 tstop 42  123 41  0  5.906168056337566  2.1788229845628453
ABM32 cb    42  120 40  0  5.178032343425095  0.3335752831335491
ABM43 plain 41  199 40  0  5.908558169283381  2.178343311730774
ABM43 tstop 42  204 41  0  5.906347881281552  2.179063879379895
ABM43 cb    42  198 40  0  5.1786057479927585 0.3335107253541619
ABM54 plain 41  198 40  0  5.9083773555293435 2.1785356616094274
ABM54 tstop 42  203 41  0  5.905722091476689  2.1793208406720583
ABM54 cb    42  196 40  0  5.178649557891228  0.3334987586341751
VCAB3 plain 443 448 442 0  5.908295845281072  2.178575599015274
VCAB3 mid   73  78  72  0  5.910657943751833  2.1795604200093974
VCAB3 loose 18  23  17  1  6.261933699139897  2.180003294251165
VCAB3 tstop 450 455 449 0  5.908295703008155  2.1785755325192686
VCAB3 cb    392 401 390 2  5.178663227838037  0.3334950025445969
VCAB4 plain 178 195 177 0  5.908286740435603  2.178574764665589
VCAB4 mid   44  61  43  0  5.907993765371094  2.179687138686123
VCAB4 loose 15  32  14  0  6.0075565336810675 2.218377798039913
VCAB4 tstop 181 198 180 0  5.908286786270012  2.178574677892119
VCAB4 cb    165 197 163 1  5.1786613775881385 0.33349561411444717
VCAB5 plain 95  117 94  0  5.908288910645755  2.178574015520456
VCAB5 mid   30  52  29  0  5.9078399961927275 2.180054174564244
VCAB5 loose 12  34  11  0  6.0602407854349805 2.2219712336031523
VCAB5 tstop 99  121 98  0  5.908288137642947  2.1785728685413344
VCAB5 cb    86  128 84  0  5.178663135713879  0.33349531378077324
VCABM3 plain 94 191 93  0  5.908411753640398  2.1785849702718507
VCABM3 mid  27  58  26  1  5.913364528666721  2.1803585935644874
VCABM3 loose 12 27  11  0  6.04183885363799   2.209354819954833
VCABM3 tstop 95 193 94  0  5.908410786460011  2.1785867397367564
VCABM3 cb   91  189 89  1  5.178634376522305  0.3334886879255225
VCABM4 plain 55 123 54  0  5.908224086698352  2.1786236495816067
VCABM4 mid  20  54  19  1  5.908813390376819  2.1842783396524053
VCABM4 loose 10 33  9   0  6.0458153804559425 2.2595219473460073
VCABM4 tstop 57 127 56  0  5.908249464984379  2.178609785460262
VCABM4 cb   52  134 50  1  5.1786668868883    0.33350477877001883
VCABM5 plain 44 105 43  0  5.908272809905584  2.1785859745087026
VCABM5 mid  19  57  18  2  5.908510718581784  2.180299595023148
VCABM5 loose 10 37  9   0  5.943273793410327  2.2109147217428244
VCABM5 tstop 44 105 43  0  5.908273021385784  2.1785858845036907
VCABM5 cb   49  130 47  0  5.178664460807733  0.33349705697192367
VCABM plain 32  68  31  3  5.908289984935434  2.1785781069863717
VCABM mid   17  38  16  3  5.906993480301624  2.179360989685579
VCABM loose 8   18  7   1  6.938860624778852  1.7477470109534596
VCABM tstop 32  69  31  4  5.908296919496497  2.1785796524950873
VCABM cb    41  88  39  6  5.178665925297634  0.3334985968040936
")

test_that("Adams-Bashforth-Moulton methods match OrdinaryDiffEqAdamsBashforthMoulton", {
  s <- abm_setup()
  for (i in seq_len(nrow(abm_ref))) {
    r <- abm_ref[i, ]
    adaptive <- startsWith(r$alg, "VC")
    sol <- do.call(solve, c(list(s$prob, get(r$alg)()), abm_kwargs(r$case, adaptive, s$cb)))
    lab <- paste(r$alg, r$case)
    expect_identical(c(length(sol$t), sol$stats$nf, sol$stats$naccept, sol$stats$nreject),
                     c(r$n, r$nf, r$na, r$nr), label = lab)
    # the final states agree to about 1e-14 (VCABM5 cb to 3e-12: the PI
    # controller of the R driver rounds q / gamma once more than Julia)
    tol <- if (!adaptive) 1e-13 else if (lab == "VCABM5 cb") 1e-11 else 1e-12
    expect_equal(sol$u[[length(sol$u)]], c(r$u1, r$u2), tolerance = tol, label = lab)
  }
})

test_that("VCABM with a fixed step and its first order step", {
  s <- abm_setup()
  sol <- solve(s$prob, VCABM(), dt = 0.05, adaptive = FALSE)
  expect_identical(c(length(sol$t), sol$stats$nf, sol$stats$naccept), c(41L, 81L, 40L))
  expect_equal(sol$u[[41]], c(5.1680079066478495, 2.4219465130857265), tolerance = 1e-13)
  # u' = 1: the order 1 corrector starts from u(t), so u(0.1) = 0.1
  one <- ODEProblem(function(u, p, t) 1, 0, c(0, 1))
  sol <- solve(one, VCABM(), dt = 0.1, adaptive = FALSE)
  expect_equal(sol(0.1), 0.1, tolerance = 1e-14)
  expect_equal(sol$u[[length(sol$u)]], 1, tolerance = 1e-14)
})
