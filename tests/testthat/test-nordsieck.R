# Reference values from OrdinaryDiffEqNordsieck 2.3.2 (Julia), printed by
# tools/nordsieck-probe.jl (its "...fix" runs, the Julia methods with the
# corrections of R/solvers-nordsieck.R): saved points, nf, naccept, nreject,
# the final state and the dense output at t = 0.005 (a), 0.3 (b) and 1.5 (c)
# on a Lotka-Volterra problem with a time dependent term on (0, 2). The
# adaptive cases take reltol = 1e-6 and abstol = 1e-8 ("plain"), 1e-4
# ("mid") or 1e-2 ("loose"); "tstop" adds a tstop at 0.73, "cb" a tstop at
# 1 where a DiscreteCallback halves u. "fixed" takes dt = 0.05 without
# adaptivity, "fixedcb" adds tstops at 0.73 and 1 and the callback. "long"
# integrates to t = 10 with tolerance 1e-4, where the Adams order of JVODE
# also decreases; "stiff" takes p[2] = 20 and tolerance 0.1, where
# functional iterations fail (AN5 once, JVODE_Adams 6 times, JVODE_BDF
# twice; as in Julia these failures do not count as rejections).

ns_setup <- function() {
  lv <- function(u, p, t) c(p[1] * u[1] - u[1] * u[2], -p[2] * u[2] + u[1] * u[2] + 0.1 * t)
  list(prob = ODEProblem(lv, c(1, 1), c(0, 2), c(1.5, 3)),
       prob10 = ODEProblem(lv, c(1, 1), c(0, 10), c(1.5, 3)),
       stiff = ODEProblem(lv, c(1, 1), c(0, 2), c(1.5, 20)),
       cb = DiscreteCallback(function(u, t, integ) t == 1, function(integ) integ$u <- integ$u * 0.5))
}

ns_kwargs <- function(case, cb) {
  tol <- switch(case, stiff = list(reltol = 1e-1, abstol = 1e-1), loose = list(reltol = 1e-2, abstol = 1e-2), mid = , long = list(reltol = 1e-4, abstol = 1e-4),
                list(reltol = 1e-6, abstol = 1e-8))
  switch(case, tstop = c(tol, list(tstops = 0.73)), cb = c(tol, list(tstops = 1, callback = cb)),
         fixed = list(dt = 0.05, adaptive = FALSE),
         fixedcb = list(dt = 0.05, adaptive = FALSE, tstops = c(0.73, 1), callback = cb), tol)
}

ns_ref <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
alg         case    n   nf  na  nr u1                 u2                  a1                 a2                 b1                 b2                 c1                 c2
AN5         plain   38  53  37  8  5.908224865337761  2.1785867441584004  1.0025280915374362 0.9900573115661637 1.248681892620222  0.5706520998351903 4.84576108142303   0.5054619094328928
AN5         mid     17  27  16  3  5.906181711449886  2.1835524418460306  1.0025280579101867 0.9900572843466353 1.2486899463730061 0.5706632860026919 4.847505641818171  0.5057190213051973
AN5         loose   8   15  7   0  6.220927231336909  2.259165821708317   1.0025278643481155 0.9900570863767096 1.248380474639066  0.5704832969069258 4.862522064731186  0.5040020816451254
AN5         tstop   40  51  39  4  5.908240046601443  2.1785782707423182  1.0025280915374362 0.9900573115661637 1.248681892620222  0.5706520998351903 4.845761746557803  0.5054572316956356
AN5         cb      36  52  34  4  5.17866890822397   0.3335000606845735  1.0025280915374362 0.9900573115661637 1.248681892620222  0.5706520998351903 2.711323304897063  0.14043738386505417
AN5         fixed   41  46  40  0  5.908181198013758  2.178534117610662   1.0025280658332298 0.9900572910969485 1.2486819799127116 0.5706527172338178 4.845753129864826  0.5054536361594312
AN5         fixedcb 43  53  41  0  5.178664974387549  0.3334977616671006  1.0025280658332298 0.9900572910969485 1.2486819799127116 0.5706527172338178 2.7113213808923917 0.14043707690689014
AN5         long    141 173 140 25 1.9667299997666752 1.4747973882030558  1.0025280579101867 0.9900572843466353 1.2486899463730061 0.5706632860026919 4.847505641818171  0.5057190213051973
AN5         stiff   60  124 59  9  13.845180333398627 0.28021801745669017 1.0027184046498878 0.9090427503468415 1.3749259524615454 -1.5040069429730192 7.079108945218507 -0.10181820899308021
JVODE_Adams plain   48  54  47  3  5.908241934581648  2.178577120646693   1.0025281696795023 0.9900574210012186 1.2486824052485423 0.5706531437760053 4.845759764742988  0.5054568805537818
JVODE_Adams mid     24  30  23  3  5.8967413306537475 2.185635313355928   1.0025486114023365 0.9900990886979915 1.2486496830616756 0.5705470152643933 4.847416145339986  0.5077021890266401
JVODE_Adams loose   14  19  13  2  6.036046391566264  2.235537250091809   1.0026062760947636 0.9902112446555702 1.252165364214628  0.5759772979347989 4.854538956236078  0.4948691980075403
JVODE_Adams tstop   49  54  48  2  5.908273406081316  2.1785666201485467  1.0025281696795023 0.9900574210012186 1.2486824052485423 0.5706531437760053 4.845754863512853  0.5054516153721361
JVODE_Adams cb      64  73  62  5  5.178658203184873  0.33350107889725383 1.0025281696795023 0.9900574210012186 1.2486824052485423 0.5706531437760053 2.711323218452681  0.14043815235335955
JVODE_Adams fixed   41  42  40  0  5.164090344984952  2.430253912666677   1.00260925         0.990216125        1.26373596813048   0.593006065693303  4.85192670332602   0.6043603609669055
JVODE_Adams fixedcb 43  45  41  0  5.301566713054654  0.4023501918955112  1.00260925         0.990216125        1.26373596813048   0.593006065693303  2.779027451148905  0.1564383809437627
JVODE_Adams long    121 132 120 8  1.966239883801778  1.4730070001945197  1.0025486114023365 0.9900990886979915 1.2486496830616756 0.5705470152643933 4.847416145339986  0.5077021890266401
JVODE_Adams stiff   34  66  33  0  17.978636136676275 0.02466908918728293 1.0033178474083975 0.9202350950172751 1.4959612145444874 0.18480492012050298 8.725368924063034 0.391958082085442
JVODE_BDF   plain   80  87  79  4  5.908274185645012  2.1785096242281967  1.002528665884793  0.9900584704166979 1.2486848019034267 0.5706523578144187 4.845720586022865  0.5054603542184504
JVODE_BDF   mid     50  56  49  3  5.900922858699407  2.1863887210818826  1.002549207986954  0.9901002910531259 1.2488292059205683 0.5698201268645967 4.851573184638444  0.5069142970297664
JVODE_BDF   loose   16  21  15  2  5.8079938929364605 2.3202295730518596  1.0026062760947636 0.9902112446555702 1.257508326725253  0.5827690616180276 4.906749011742782  0.5058608454125688
JVODE_BDF   tstop   82  89  81  4  5.908283639422584  2.178518015160357   1.002528665884793  0.9900584704166979 1.2486848019034267 0.5706523578144187 4.845727628861139  0.5054586444277397
JVODE_BDF   cb      91  100 89  5  5.178659806795965  0.3335007289779083  1.002528665884793  0.9900584704166979 1.2486848019034267 0.5706523578144187 2.711325979191178  0.14043765353776475
JVODE_BDF   fixed   41  42  40  0  5.661141285188454  2.252043117652268   1.00260925         0.990216125        1.26373596813048   0.593006065693303  4.829925494897085  0.5606991136612632
JVODE_BDF   fixedcb 43  45  41  0  5.315346725164629  0.40934474661317444 1.00260925         0.990216125        1.26373596813048   0.593006065693303  2.7921539732737255 0.15892733421646815
JVODE_BDF   long    134 143 133 6  1.9650699061511707 1.4760479881026487  1.002549207986954  0.9901002910531259 1.2488292059205683 0.5698201268645967 4.851573184638444  0.5069142970297664
JVODE_BDF   stiff   48  61  47  0  19.262623216562965 0.043973968525846334 1.0033178474083975 0.9202350950172751 1.4925887510621536 0.08272883039551107 8.957727812653758 0.016293112074465516
")

test_that("Nordsieck methods match OrdinaryDiffEqNordsieck", {
  s <- ns_setup()
  for (i in seq_len(nrow(ns_ref))) {
    r <- ns_ref[i, ]
    prob <- switch(r$case, long = s$prob10, stiff = s$stiff, s$prob)
    sol <- do.call(solve, c(list(prob, get(r$alg)()), ns_kwargs(r$case, s$cb)))
    lab <- paste(r$alg, r$case)
    expect_identical(c(length(sol$t), sol$stats$nf, sol$stats$naccept, sol$stats$nreject),
                     c(r$n, r$nf, r$na, r$nr), label = lab)
    # the final states and dense outputs agree to about 5e-15
    expect_equal(sol$u[[length(sol$u)]], c(r$u1, r$u2), tolerance = 1e-12, label = lab)
    expect_equal(c(sol(0.005), sol(0.3), sol(1.5)), c(r$a1, r$a2, r$b1, r$b2, r$c1, r$c2),
                 tolerance = 1e-12, label = lab)
  }
})

test_that("a failed functional iteration without adaptivity ends the solve", {
  # lines "F" of tools/nordsieck-probe.jl: dt = 0.2 on the "stiff" problem
  s <- ns_setup()
  for (a in list(AN5(), JVODE_Adams(), JVODE_BDF())) {
    sol <- suppressWarnings(solve(s$stiff, a, dt = 0.2, adaptive = FALSE))
    expect_identical(sol$retcode, "ConvergenceFailure")
    expect_identical(sol$stats$nreject, 0L)
    if (inherits(a, "AN5")) {
      expect_identical(c(length(sol$t), sol$stats$nf, sol$stats$naccept), c(2L, 9L, 1L))
      expect_equal(sol$u[[2]], c(1.780662006836855, 1.8602941944599025), tolerance = 1e-14)
    } else {
      expect_identical(c(length(sol$t), sol$stats$nf, sol$stats$naccept), c(1L, 4L, 0L))
    }
  }
})

test_that("JVODE selects the Adams or BDF formulas", {
  s <- ns_setup()
  ref <- solve(s$prob, JVODE_BDF(), reltol = 1e-4, abstol = 1e-4)
  sol <- solve(s$prob, JVODE(algorithm = "BDF"), reltol = 1e-4, abstol = 1e-4)
  expect_identical(sol$t, ref$t)
  expect_identical(sol$u, ref$u)
  sol <- solve(s$prob, JVODE(), reltol = 1e-4, abstol = 1e-4)
  expect_identical(c(length(sol$t), sol$stats$nf), c(24L, 30L))
})
