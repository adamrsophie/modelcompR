# Reference values from tools/exponentialrk-probe.jl (Julia,
# OrdinaryDiffEqExponentialRK): final states and statistics (nf, nf2, njacs,
# naccept, nreject). With dt = 0.1 the classical methods cannot change dt
# upstream and step past tf, then append tf; the values are those of the
# step before, and the counts leave that extra step out.
.EXPRK_REF <- list(
  `sp1_LawsonEuler` = list(u = c(0.16918173581494322, 0.1479958823573061, -0.1335572580803068, -0.002904450331952464),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `sp3_LawsonEuler` = list(u = c(0.28737102905594175, 0.08882148947786711, -0.31303938503062545, 0.0037126228134732983),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `sp1_NorsettEuler` = list(u = c(0.16910672659327258, 0.1499532849922745, -0.13340316764688598, -0.0039648868470729603),
    stats = c(11, 0, 0, 10, 0), nt = 11, tend = 1),
  `sp3_NorsettEuler` = list(u = c(0.28525740882018391, 0.091194667475786181, -0.31331079215316965, 0.0026675612714361112),
    stats = c(11, 0, 0, 10, 0), nt = 11, tend = 1),
  `sp1_ETDRK2` = list(u = c(0.17788576000959777, 0.14883662297706671, -0.13434080085526018, -0.0049710836070923289),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `sp3_ETDRK2` = list(u = c(0.30208497744389068, 0.088491717935511069, -0.31456109905135804, 0.0013589268894775314),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `sp1_ETDRK3` = list(u = c(0.17838074933982523, 0.14896209489444495, -0.13434200404897684, -0.0049156878894026307),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `sp3_ETDRK3` = list(u = c(0.30235511599405362, 0.088375529717632489, -0.31468786636089086, 0.0013771341744005152),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `sp1_ETDRK4` = list(u = c(0.17835418128044517, 0.14893835947003498, -0.13434861482831939, -0.0049190162681039453),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp3_ETDRK4` = list(u = c(0.3023373143482217, 0.088371636494989914, -0.31468965963829948, 0.0013770507238378433),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp1_HochOst4` = list(u = c(0.17835754697618647, 0.14894059059555728, -0.13434813826810416, -0.0049186961900271877),
    stats = c(11, 40, 0, 10, 0), nt = 11, tend = 1),
  `sp3_HochOst4` = list(u = c(0.30234195765833677, 0.088372049584900608, -0.31468967692556327, 0.0013770579662934559),
    stats = c(11, 40, 0, 10, 0), nt = 11, tend = 1),
  `sp1_Friedli` = list(u = c(0.17835482605400929, 0.14893885490151673, -0.13434853428125973, -0.0049189679259256554),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp3_Friedli` = list(u = c(0.3023385645358338, 0.088371847895754224, -0.31468967915758994, 0.0013770558893801896),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp1k_LawsonEuler` = list(u = c(0.16920783680956, 0.14802670943020604, -0.1335156917004271, -0.002874126707875531),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_LawsonEuler` = list(u = c(0.15321673743303385, 0.14527159687948987, -0.13269182526248996, -0.00086244217607597839),
    stats = c(5, 4, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_LawsonEuler` = list(u = c(0.16961910040877637, 0.18039949436402083, -0.10233567379733814, -0.0021270431787732678),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `np1k_LawsonEuler` = list(u = c(0.16588913815578599, 0.14962503273835934, -0.13467044029250569, -0.0028553684354119106),
    stats = c(21, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1k_NorsettEuler` = list(u = c(0.16910605961171904, 0.14995227444296241, -0.13340477531358921, -0.0039659005387507459),
    stats = c(11, 0, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_NorsettEuler` = list(u = c(0.14868387882898315, 0.14989116952521661, -0.13221915420417391, -0.0028644476995944856),
    stats = c(5, 0, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_NorsettEuler` = list(u = c(0.16955180974571604, 0.18237880705316009, -0.10223236718040728, -0.0031272291205412493),
    stats = c(11, 0, 0, 10, 0), nt = 11, tend = 1),
  `np1k_NorsettEuler` = list(u = c(0.17750643919801171, 0.14911290381853021, -0.1339593918868297, -0.003590388272497567),
    stats = c(11, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1k_ETDRK2` = list(u = c(0.17788509518888002, 0.14883563008497222, -0.13434248986610722, -0.0049721167940934624),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_ETDRK2` = list(u = c(0.17411388857505725, 0.14764594951173468, -0.13443513293447409, -0.0053469836199764829),
    stats = c(5, 4, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_ETDRK2` = list(u = c(0.17858692581041707, 0.18241071142849816, -0.10263980219983286, -0.0043947306905330061),
    stats = c(11, 10, 0, 10, 0), nt = 11, tend = 1),
  `np1k_ETDRK2` = list(u = c(0.17799123714366277, 0.14883981393998044, -0.13440749573876579, -0.0049756187020772655),
    stats = c(21, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1k_ETDRK3` = list(u = c(0.17838007825613378, 0.14896110031497961, -0.13434368365808955, -0.0049167222829933381),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_ETDRK3` = list(u = c(0.17849707721460981, 0.1490367413900098, -0.13431446370123834, -0.0048995929196347758),
    stats = c(5, 8, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_ETDRK3` = list(u = c(0.17910163808611695, 0.18264255365090737, -0.10258279209686938, -0.0043411876375227037),
    stats = c(11, 20, 0, 10, 0), nt = 11, tend = 1),
  `np1k_ETDRK3` = list(u = c(0.17835291151643326, 0.14893611703708523, -0.13435165816159242, -0.0049195141199288289),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1k_ETDRK4` = list(u = c(0.17835350908937819, 0.14893737151947142, -0.13435028131219276, -0.0049200479741814138),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_ETDRK4` = list(u = c(0.1783211398140406, 0.14886183619986557, -0.13436848154177877, -0.0049276707589027344),
    stats = c(5, 12, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_ETDRK4` = list(u = c(0.17907463298866633, 0.18261153603568311, -0.10259497010185736, -0.0043453677867454745),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `np1k_ETDRK4` = list(u = c(0.1783550722846271, 0.14893744959774496, -0.13435150038935695, -0.0049193523314040657),
    stats = c(41, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1k_Friedli` = list(u = c(0.1783541562785366, 0.14893786399750092, -0.13435021073418879, -0.0049200008610611381),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `sp1kd_Friedli` = list(u = c(0.17836339530549383, 0.14889167866002706, -0.13436328316587523, -0.0049243472518143566),
    stats = c(5, 12, 0, 4, 0), nt = 5, tend = 1),
  `sp2k_Friedli` = list(u = c(0.17907526803848273, 0.182612067144242, -0.10259480498012516, -0.0043453031853960301),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `np1k_Friedli` = list(u = c(0.17835524718115461, 0.14893755020343727, -0.13435148139841199, -0.0049193382951272149),
    stats = c(41, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1_ETD2` = list(u = c(0.17979279953505134, 0.14863569777132246, -0.13430360924460005, -0.0047549508423255721),
    stats = c(11, 11, 0, 10, 0), nt = 11, tend = 1),
  `sp3_ETD2` = list(u = c(0.31196658971806662, 0.08559920036578475, -0.31426524283525548, 0.0020683170449039209),
    stats = c(5, 5, 0, 4, 0), nt = 5, tend = 1),
  `np1_Exp4` = list(u = c(0.17977979792880397, 0.14909437519671517, -0.13453297894877284, -0.0058236989840665821),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_Exp4` = list(u = c(0.17977979802724273, 0.14909437501260653, -0.13453298125063717, -0.0058236967352631348),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_Exp4` = list(u = c(0.1819581639513761, 0.14922712104357727, -0.13484707200989815, -0.0075412577939819279),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EPIRK4s3A` = list(u = c(0.17802558633689339, 0.14891300642326674, -0.13430738027480277, -0.0047187548063211031),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EPIRK4s3A` = list(u = c(0.17802558633842608, 0.14891300642257732, -0.13430738027408823, -0.0047187548060742416),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EPIRK4s3A` = list(u = c(0.17860997332356676, 0.14937747203357191, -0.13419348675928777, -0.0043749843943010659),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EPIRK4s3B` = list(u = c(0.17725058272755964, 0.14882647993383102, -0.13427214288353623, -0.0043561457299688567),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EPIRK4s3B` = list(u = c(0.17725058273645278, 0.1488264799268364, -0.13427214288092493, -0.0043561457283002852),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EPIRK4s3B` = list(u = c(0.17691812419941155, 0.14943044486888071, -0.13425686717592752, -0.0036977102067036263),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EPIRK5s3` = list(u = c(0.18050095465828925, 0.14873489728360581, -0.13498323570153384, -0.0073706203080104874),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EPIRK5s3` = list(u = c(0.1805009546584492, 0.14873489728404232, -0.13498323570090973, -0.0073706203080785302),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EPIRK5s3` = list(u = c(0.16994359007936621, 0.14315182545898486, -0.13691075425450555, -0.012126320117717849),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EXPRB53s3` = list(u = c(0.17828820966953968, 0.14893059264111777, -0.13433798304490632, -0.0048823024334900805),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EXPRB53s3` = list(u = c(0.17828820964627637, 0.14893059262236319, -0.13433798305163064, -0.0048823024356115952),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EXPRB53s3` = list(u = c(0.17813593142592909, 0.14889732109204998, -0.13432430231371376, -0.0048605010704933851),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EPIRK5P1` = list(u = c(0.17835842218441872, 0.14894080818479313, -0.13434854853921502, -0.0049194719761747017),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EPIRK5P1` = list(u = c(0.17835842392266654, 0.14894081472936249, -0.13434853622455356, -0.0049194761602231541),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EPIRK5P1` = list(u = c(0.17830593180725793, 0.14888390925217374, -0.13440918475289176, -0.0049249939579378782),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1_EPIRK5P2` = list(u = c(0.17830151209644221, 0.14893266798234156, -0.13433979463849408, -0.0048891696726780717),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1t_EPIRK5P2` = list(u = c(0.17830151373618902, 0.14893267450089129, -0.13433978268378555, -0.0048891735088461815),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1n_EPIRK5P2` = list(u = c(0.17820927351859953, 0.14888832124299275, -0.13438014602398951, -0.0048772337905000959),
    stats = c(13, 0, 4, 4, 0), nt = 5, tend = 1),
  `np1a_Exprb32` = list(u = c(0.17813968942133057, 0.14890602270782832, -0.13431102944109308, -0.0048144269240111455),
    stats = c(83, 0, 40, 36, 4), nt = 37, tend = 1),
  `np1af_Exprb32` = list(u = c(0.17781999413615895, 0.14892583410258323, -0.13425201843014989, -0.0045327190427737344),
    stats = c(21, 0, 10, 10, 0), nt = 11, tend = 1),
  `np1a_Exprb43` = list(u = c(0.17835600872067517, 0.14894040855670723, -0.13434797770817034, -0.0049193332433221208),
    stats = c(141, 0, 46, 44, 2), nt = 45, tend = 1),
  `np1af_Exprb43` = list(u = c(0.17837248087497293, 0.14894225014240373, -0.13434552503630162, -0.0049282173236268731),
    stats = c(31, 0, 10, 10, 0), nt = 11, tend = 1),
  `npfdk_ETDRK4` = list(u = c(0.17835507228451972, 0.14893744959769717, -0.13435150038938357, -0.0049193523314008634),
    stats = c(91, 0, 10, 10, 0), nt = 11, tend = 1),
  `npfdk_LawsonEuler` = list(u = c(0.14454826562241385, 0.15024769790738016, -0.13501594387712143, -0.00057934144253967291),
    stats = c(29, 0, 4, 4, 0), nt = 5, tend = 1),
  `npfd_EPIRK4s3A` = list(u = c(0.17802558638539021, 0.14891300644548611, -0.13430738027073505, -0.0047187548008920319),
    stats = c(81, 0, 10, 10, 0), nt = 11, tend = 1),
  `npfd_EPIRK5P2` = list(u = c(0.17830151210234085, 0.14893266798439997, -0.13433979463921308, -0.0048891696722923378),
    stats = c(81, 0, 10, 10, 0), nt = 11, tend = 1),
  `sp1fd_EPIRK4s3A` = list(u = c(0.17802558638539021, 0.14891300644548611, -0.13430738027073505, -0.0047187548008920319),
    stats = c(81, 0, 10, 10, 0), nt = 11, tend = 1),
  `npfda_Exprb32` = list(u = c(0.17813968946415948, 0.14890602272377232, -0.13431102944661627, -0.0048144269273919429),
    stats = c(283, 0, 40, 36, 4), nt = 37, tend = 1),
  `npfda_Exprb43` = list(u = c(0.17835600872076093, 0.14894040855670226, -0.13434797770815024, -0.0049193332434125598),
    stats = c(371, 0, 46, 44, 2), nt = 45, tend = 1),
  `spu_ETDRK4` = list(u = c(0.17887415552918404, 0.14633043055810044, -0.13810630393310816, -0.0052600694683922504),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `spuk_ETDRK4` = list(u = c(0.17887426011592197, 0.14632997884482804, -0.13810665864336724, -0.0052603416363846212),
    stats = c(11, 30, 0, 10, 0), nt = 11, tend = 1),
  `spuk_LawsonEuler` = list(u = c(0.15321673743303385, 0.14527159687948987, -0.13269182526248996, -0.00086244217607597839),
    stats = c(5, 4, 0, 4, 0), nt = 5, tend = 1),
  `spu_EPIRK4s3A` = list(u = c(0.182063733922949, 0.13263384558146693, -0.1601226866710006, -0.0067519707158721631),
    stats = c(81, 0, 10, 10, 0), nt = 11, tend = 1)
)

exprk_A4 <- matrix(c(-2, 0.5, 0.1, 0, 0.3, -1.5, 0.4, 0.2, 0, 0.2, -1, 0.5, 0.1, 0, 0.3, -2.5), 4, byrow = TRUE)
exprk_f2 <- function(u, p, t) c(-u[1] * u[1] * u[1] + 0.1 * t, 0.5 * u[1] * u[2], -0.2 * u[3] * u[3] + 0.3 * u[4],
                                0.1 * u[4] * u[1] - 0.05 * t * t)
exprk_jf2 <- function(u, p, t) matrix(c(-3 * u[1] * u[1], 0, 0, 0,
                                        0.5 * u[2], 0.5 * u[1], 0, 0,
                                        0, 0, -0.4 * u[3], 0.3,
                                        0.1 * u[4], 0, 0, 0.1 * u[1]), 4, byrow = TRUE)
exprk_u0 <- c(1, 0.5, -0.5, 0.2)
exprk_sp1 <- SplitODEProblem(MatrixOperator(exprk_A4), exprk_f2, exprk_u0, c(0, 1))
exprk_sp2 <- SplitODEProblem(MatrixOperator((exprk_A4 + t(exprk_A4)) / 2), exprk_f2, exprk_u0, c(0, 1))
exprk_sp3 <- SplitODEProblem(DiagonalOperator(c(-1, -2, -0.5, -3)), exprk_f2, exprk_u0, c(0, 1))
exprk_fn <- function(u, p, t) as.numeric(exprk_A4 %*% u) + exprk_f2(u, p, t)
exprk_np1 <- ODEProblem(ODEFunction(exprk_fn, jac_fn = function(u, p, t) exprk_A4 + exprk_jf2(u, p, t)),
                        exprk_u0, c(0, 1))

# the final state, the statistics (nf, nf2, njacs, naccept, nreject), the
# number of saved points and the final time against the reference
exprk_check <- function(sol, nm, tol = 1e-12) {
  ref <- .EXPRK_REF[[nm]]
  expect_length(sol$u, ref$nt)
  expect_equal(sol$t[[ref$nt]], ref$tend, tolerance = 1e-14, label = nm)
  expect_equal(as.numeric(sol$u[[ref$nt]]), ref$u, tolerance = tol, label = nm)
  st <- sol$stats
  expect_equal(c(st$nf, st$nf2 %||% 0, st$njacs, st$naccept, st$nreject), ref$stats, label = nm)
}

exprk_classical <- c("LawsonEuler", "NorsettEuler", "ETDRK2", "ETDRK3", "ETDRK4", "HochOst4", "Friedli")

test_that("classical exponential Runge-Kutta methods match OrdinaryDiffEqExponentialRK", {
  for (nm in exprk_classical) {
    alg <- get(nm)
    exprk_check(solve(exprk_sp1, alg(), dt = 0.1), paste0("sp1_", nm))
    exprk_check(solve(exprk_sp3, alg(), dt = 0.1), paste0("sp3_", nm))
    if (nm == "HochOst4") next
    exprk_check(solve(exprk_sp1, alg(krylov = TRUE, m = 3, iop = 1), dt = 0.1), paste0("sp1k_", nm))
    exprk_check(solve(exprk_sp1, alg(krylov = TRUE), dt = 0.25), paste0("sp1kd_", nm))
    exprk_check(solve(exprk_sp2, alg(krylov = TRUE, m = 3), dt = 0.1), paste0("sp2k_", nm))
    exprk_check(solve(exprk_np1, alg(krylov = TRUE, m = 3), dt = 0.1), paste0("np1k_", nm))
  }
  exprk_check(solve(exprk_sp1, ETD2(), dt = 0.1), "sp1_ETD2")
  exprk_check(solve(exprk_sp3, ETD2(), dt = 0.25), "sp3_ETD2")
  expect_identical(ETD1()$name, NorsettEuler()$name)
})

test_that("Krylov forms agree with the precomputed operators", {
  # with a full Krylov space the two forms are the same method; upstream's
  # Krylov HochOst4 repeats the second stage in its fifth stage
  last <- function(s) s$u[[length(s$u)]]
  for (nm in exprk_classical)
    expect_equal(last(solve(exprk_sp1, get(nm)(krylov = TRUE), dt = 0.25)),
                 last(solve(exprk_sp1, get(nm)(), dt = 0.25)), tolerance = 1e-13, label = nm)
})

test_that("a last step shortened to the final time gets its own operators", {
  g <- function(u, p, t) c(0.1, 0, -0.2, 0.05)
  prob <- SplitODEProblem(MatrixOperator(exprk_A4), g, exprk_u0, c(1, 0))
  E <- .eu_exp_base(-exprk_A4)
  exact <- as.numeric(E %*% exprk_u0 + solve(exprk_A4, (E - diag(4)) %*% g()))
  for (nm in c("NorsettEuler", "ETDRK4", "HochOst4", "ETD2")) {
    sol <- solve(prob, get(nm)(), dt = 0.3)
    expect_equal(sol$t, c(1, 0.7, 0.4, 0.1, 0), tolerance = 1e-14)
    expect_equal(sol$u[[5]], exact, tolerance = 1e-13, label = nm)
  }
})

test_that("EPIRK methods match OrdinaryDiffEqExponentialRK", {
  for (nm in c("Exp4", "EPIRK4s3A", "EPIRK4s3B", "EPIRK5s3", "EXPRB53s3", "EPIRK5P1", "EPIRK5P2")) {
    alg <- get(nm)
    exprk_check(solve(exprk_np1, alg(), dt = 0.1), paste0("np1_", nm))
    exprk_check(solve(exprk_np1, alg(m = 3, iop = 1), dt = 0.1, reltol = 1e-8), paste0("np1t_", nm))
    exprk_check(solve(exprk_np1, alg(adaptive_krylov = FALSE, m = 3), dt = 0.25), paste0("np1n_", nm))
  }
})

test_that("exponential Rosenbrock methods match OrdinaryDiffEqExponentialRK", {
  for (nm in c("Exprb32", "Exprb43")) {
    alg <- get(nm)
    exprk_check(solve(exprk_np1, alg()), paste0("np1a_", nm))
    exprk_check(solve(exprk_np1, alg(), dt = 0.1, adaptive = FALSE), paste0("np1af_", nm))
    # a split problem is linearized as a whole, like the same ODEProblem
    # (with a finite difference Jacobian here)
    s1 <- solve(exprk_sp1, alg())
    s2 <- solve(exprk_np1, alg())
    expect_length(s1$t, length(s2$t))
    expect_equal(s1$u[[length(s1$u)]], s2$u[[length(s2$u)]], tolerance = 1e-8, label = nm)
  }
})

test_that("the finite difference Jacobian matches AutoFiniteDiff upstream", {
  # the forward differences are those of FiniteDiff.jl; the right-hand side
  # A u differs from Julia's in the last bit at some states (the order of
  # the sums), which the differences magnify by 1 / sqrt(eps)
  npfd <- ODEProblem(exprk_fn, exprk_u0, c(0, 1))
  exprk_check(solve(npfd, ETDRK4(krylov = TRUE, m = 3), dt = 0.1), "npfdk_ETDRK4")
  exprk_check(solve(npfd, LawsonEuler(krylov = TRUE), dt = 0.25), "npfdk_LawsonEuler", tol = 1e-9)
  exprk_check(solve(npfd, EPIRK4s3A(), dt = 0.1), "npfd_EPIRK4s3A", tol = 1e-9)
  exprk_check(solve(npfd, EPIRK5P2(autodiff = AutoFiniteDiff()), dt = 0.1), "npfd_EPIRK5P2", tol = 1e-9)
  exprk_check(solve(exprk_sp1, EPIRK4s3A(), dt = 0.1), "sp1fd_EPIRK4s3A", tol = 1e-9)
  exprk_check(solve(npfd, Exprb32()), "npfda_Exprb32", tol = 1e-6)
  exprk_check(solve(npfd, Exprb43()), "npfda_Exprb43", tol = 1e-6)
})

test_that("an operator with an update function is used as constructed, as upstream out of place", {
  upd <- function(A, u, p, t) exprk_A4 + t * diag(c(0.1, -0.2, 0.3, -0.1))
  spu <- SplitODEProblem(MatrixOperator(exprk_A4, upd), exprk_f2, exprk_u0, c(0, 1))
  exprk_check(solve(spu, ETDRK4(), dt = 0.1), "spu_ETDRK4")
  exprk_check(solve(spu, ETDRK4(krylov = TRUE, m = 3), dt = 0.1), "spuk_ETDRK4")
  exprk_check(solve(spu, LawsonEuler(krylov = TRUE), dt = 0.25), "spuk_LawsonEuler")
  exprk_check(solve(spu, EPIRK4s3A(), dt = 0.1), "spu_EPIRK4s3A", tol = 1e-9)
})

test_that("exponential methods check their problems", {
  expect_error(solve(exprk_np1, ETDRK4(), dt = 0.1), "SplitFunction")
  expect_error(solve(exprk_np1, ETD2(), dt = 0.1), "SplitODEProblem")
  prob <- SplitODEProblem(function(u, p, t) -u, exprk_f2, exprk_u0, c(0, 1))
  expect_error(solve(prob, ETDRK4(), dt = 0.1), "SciML operator")
})

test_that("EPIRK methods stop backward in time instead of returning garbage", {
  # upstream (tools/review-tdir-probe.jl) returns u0 unchanged, NaN or other
  # uninitialized values, varying between runs, often with retcode Success
  prob <- ODEProblem(function(u, p, t) -u, exp(-1), c(1, 0))
  for (alg in list(Exp4(), EPIRK4s3A(), EPIRK4s3B(), EPIRK5s3(), EXPRB53s3(), EPIRK5P1(), EPIRK5P2()))
    expect_error(solve(prob, alg, dt = -0.1), "cannot solve backward in time")
  expect_error(.eu_phiv_timestep(c(-0.5, -1), diag(2), diag(2)), "nonnegative times")
  expect_equal(tail(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), Exp4(), dt = 0.1)$u, 1)[[1]],
               exp(-1), tolerance = 1e-6)
})
