# The symplectic methods (R/solvers-symplectic.R) against
# OrdinaryDiffEqSymplecticRK 2.3.0: the reference values come from
# tools/rkn-probe.jl (out of place functions, state c(v, u)).

symp_test_probs <- list(
  P1 = function() DynamicalODEProblem(function(v, u, p, t) -u - u * u * u, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)),
  P3 = function() DynamicalODEProblem(function(v, u, p, t) -u - u * u * u - 0.1 * v, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)),
  P4 = function() DynamicalODEProblem(function(v, u, p, t) -u - u * u * u, function(v, u, p, t) 0.5 * v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)),
  P5 = function() DynamicalODEProblem(function(v, u, p, t) -u + t, function(v, u, p, t) v,
                                      c(0.5, -0.3), c(1, 0.2), c(0, 1)))

symp_refs <- list(
  list(alg = "SymplecticEuler", case = "P1", dt = 0.1, nf = 12, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.261057753557228, -0.3180670500319493, 0.5543336153483969, -0.13706077366240033),
       tq = 0.55, uq = c(-0.7198522110305612, -0.3534014436008864, 0.9931636960796856, 0.01829376849519228)),
  list(alg = "SymplecticEuler", case = "P5", dt = 0.1, nf = 12, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.17593229300305843, 0.1000170325110609, 1.1189224531862498, -0.020948873821482916),
       tq = 0.55, uq = c(0.011038071547687354, -0.23118436246561244, 1.154609453564869, 0.032898140708978736)),
  list(alg = "VelocityVerlet", case = "P1", dt = 0.1, nf = 11, nf2 = 1, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2061897216060409, -0.33056773587409505, 0.5065205331625162, -0.14579339862450158),
       tq = 0.55, uq = c(-0.6640833853096322, -0.36131485430084004, 0.9493210687732231, 0.012882382651317588)),
  list(alg = "VelocityVerlet", case = "P5", dt = 0.1, nf = 11, nf2 = 1, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.11061913790160416, 0.12993467129946906, 1.1185760567305767, 0.012414745259078241),
       tq = 0.55, uq = c(0.05156817311805606, -0.21257266512153875, 1.14063437818641, 0.039883184835039415)),
  list(alg = "VerletLeapfrog", case = "P1", dt = 0.1, nf = 11, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2061897216060407, -0.33056773587409505, 0.5065205331625162, -0.1457933986245016),
       tq = 0.55, uq = c(-0.6610044766075452, -0.3608622724698778, 0.9470714395476934, 0.012866281350848652)),
  list(alg = "VerletLeapfrog", case = "P5", dt = 0.1, nf = 11, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.11061913790160377, 0.129934671299469, 1.1185760567305765, 0.012414745259078234),
       tq = 0.55, uq = c(0.052755196681112476, -0.21105505227407753, 1.1398970069277499, 0.04052003479154998)),
  list(alg = "LeapfrogDriftKickDrift", case = "P1", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2132356861585696, -0.3310095575811115, 0.503333562600105, -0.14511917332648913),
       tq = 0.55, uq = c(-0.6698086088115676, -0.36109964475569556, 0.944770213004922, 0.013288239670992259)),
  list(alg = "LeapfrogDriftKickDrift", case = "P5", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.11272601387261805, 0.12951329610526618, 1.1196294947160834, 0.015153684021397345),
       tq = 0.55, uq = c(0.05144798966737513, -0.21131649367682495, 1.1405506104346181, 0.042219403909408745)),
  list(alg = "PseudoVerletLeapfrog", case = "P1", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2132356861585696, -0.3310095575811115, 0.503333562600105, -0.14511917332648913),
       tq = 0.55, uq = c(-0.6728888022971045, -0.3615514055593133, 0.9470261061060938, 0.013304894562924816)),
  list(alg = "PseudoVerletLeapfrog", case = "P5", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.1548635332929117, 0.08737577668497291, 1.0966270572627592, -0.007848753431927282),
       tq = 0.55, uq = c(0.02408547920436574, -0.2390074597440145, 1.1339109112829406, 0.034203349862335586)),
  list(alg = "McAte2", case = "P1", dt = 0.1, nf = 31, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.208971168251283, -0.3308262404374129, 0.5058554039896755, -0.14538351148764128),
       tq = 0.55, uq = c(-0.6674806540468347, -0.3614311209057323, 0.9484426002976938, 0.013137896165004716)),
  list(alg = "McAte2", case = "P5", dt = 0.1, nf = 31, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.13630466736265728, 0.10483952351849124, 1.1058540157093786, 0.0004986261054878323),
       tq = 0.55, uq = c(0.03566764327097157, -0.2280885493838829, 1.1367422694572138, 0.036553059227863546)),
  list(alg = "Ruth3", case = "P1", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083013649868397, -0.3308589454819584, 0.5067321502027955, -0.14531325767294573),
       tq = 0.55, uq = c(-0.666593542473854, -0.36142741615873053, 0.9487345951490415, 0.013180167812616318)),
  list(alg = "Ruth3", case = "P5", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.1925502496276994, 0.04836637908520236, 1.0788794386789842, -0.026520130167225053),
       tq = 0.55, uq = c(0.0012291948819042695, -0.26265191388987136, 1.128939706908539, 0.028782423925865337)),
  list(alg = "McAte3", case = "P1", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083386165156191, -0.33085382558287113, 0.5068258251460502, -0.14530976464559153),
       tq = 0.55, uq = c(-0.6666347079412013, -0.3614240817796939, 0.9487988294772237, 0.013182571035100683)),
  list(alg = "McAte3", case = "P5", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.15134890760541067, 0.08957806756171276, 1.096900758721698, -0.008509947987343145),
       tq = 0.55, uq = c(0.02625336987853174, -0.23762132656059856, 1.1336851130152203, 0.0335209349496838)),
  list(alg = "CandyRoz4", case = "P1", dt = 0.1, nf = 51, nf2 = 41, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083018928348075, -0.3308531975563514, 0.5069193151014074, -0.14530446319384271),
       tq = 0.55, uq = c(-0.6665926359127505, -0.3614226688971947, 0.948861955381959, 0.013185437372629156)),
  list(alg = "CandyRoz4", case = "P5", dt = 0.1, nf = 51, nf2 = 41, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.15258303299593295, 0.08834716594802074, 1.0971902289332252, -0.008227420295180941),
       tq = 0.55, uq = c(0.025607961023677548, -0.23826421232161066, 1.1340013715598303, 0.033834806688973616)),
  list(alg = "McAte4", case = "P1", dt = 0.1, nf = 51, nf2 = 41, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083610633025181, -0.33085110179528404, 0.5068432969237597, -0.14530799856358406),
       tq = 0.55, uq = c(-0.6666709320590303, -0.3614222885992872, 0.9488346165544769, 0.013183398856756972)),
  list(alg = "McAte4", case = "P5", dt = 0.1, nf = 51, nf2 = 41, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.13186985760159015, 0.10906510577940998, 1.1081734080581922, 0.002754803135550507),
       tq = 0.55, uq = c(0.0384317468945451, -0.22543800541769624, 1.137448826726787, 0.03727973275857657)),
  list(alg = "CalvoSanz4", case = "P1", dt = 0.1, nf = 51, nf2 = 51, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083555244174788, -0.3308514999841825, 0.5068249193980983, -0.14530839290673456),
       tq = 0.55, uq = c(-0.666668482141321, -0.36142253900189986, 0.9488313467817493, 0.013182969273639556)),
  list(alg = "CalvoSanz4", case = "P5", dt = 0.1, nf = 51, nf2 = 51, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.13206722267784768, 0.1088677346748198, 1.1087885428321071, 0.003369931572624844),
       tq = 0.55, uq = c(0.03840051100206961, -0.2254692453646554, 1.1378133989736714, 0.03764429997266046)),
  list(alg = "McAte42", case = "P1", dt = 0.1, nf = 61, nf2 = 51, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.218219915279062, -0.3282967565720282, 0.5141834442716314, -0.14393242487094166),
       tq = 0.55, uq = c(-0.6766943477125934, -0.3598732246060503, 0.9566803905626984, 0.014037018137566056)),
  list(alg = "McAte42", case = "P5", dt = 0.1, nf = 61, nf2 = 51, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.16126834045459854, 0.0863861015213679, 1.099208425581031, -0.012659877649011858),
       tq = 0.55, uq = c(0.020256124877264463, -0.23949247567298035, 1.1373437452665525, 0.033020192592383)),
  list(alg = "McAte5", case = "P1", dt = 0.1, nf = 71, nf2 = 61, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083551132349524, -0.3308514913514574, 0.5068250686899355, -0.14530841337205638),
       tq = 0.55, uq = c(-0.6666679417342838, -0.3614225355961255, 0.9488315380617033, 0.013182955410036886)),
  list(alg = "McAte5", case = "P5", dt = 0.1, nf = 71, nf2 = 61, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.16001098909167086, 0.0809239553080463, 1.0951247800880897, -0.010293853725326792),
       tq = 0.55, uq = c(0.021245436478110967, -0.2426243285477173, 1.1338004899142762, 0.03363137478586851)),
  list(alg = "Yoshida6", case = "P1", dt = 0.1, nf = 81, nf2 = 81, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083578485808384, -0.3308513227819585, 0.5068317718223128, -0.14530826370658378),
       tq = 0.55, uq = c(-0.6666696878063856, -0.3614224333469482, 0.9488325396633164, 0.013183124593727313)),
  list(alg = "Yoshida6", case = "P5", dt = 0.1, nf = 81, nf2 = 81, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.3580552420315543, -0.11712029337739195, 0.9849397924795978, -0.12047883781688323),
       tq = 0.55, uq = c(-0.1020233311106954, -0.365893093796717, 1.0979910344630528, -0.002178079588800728)),
  list(alg = "KahanLi6", case = "P1", dt = 0.1, nf = 101, nf2 = 101, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083574314996635, -0.33085132287503144, 0.5068322799491403, -0.14530826299368435),
       tq = 0.55, uq = c(-0.6666689975546312, -0.3614224329414246, 0.9488327685841333, 0.013183125017795928)),
  list(alg = "KahanLi6", case = "P5", dt = 0.1, nf = 101, nf2 = 101, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.20332018993665535, 0.03761475382482601, 1.0694719390055099, -0.03594669327803788),
       tq = 0.55, uq = c(-0.005908293489668461, -0.2697780589424184, 1.1251098079451134, 0.024940693647564682)),
  list(alg = "McAte8", case = "P1", dt = 0.1, nf = 161, nf2 = 161, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083573836078572, -0.3308513228850246, 0.5068323376403419, -0.14530826292524265),
       tq = 0.55, uq = c(-0.6666689134281947, -0.36142243290896886, 0.948832795717532, 0.013183125056486797)),
  list(alg = "McAte8", case = "P5", dt = 0.1, nf = 161, nf2 = 161, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.2549047925648707, -0.013969849413045868, 1.041291142382908, -0.06412749015794125),
       tq = 0.55, uq = c(-0.037950522371391376, -0.30182028816980944, 1.1160691215517127, 0.015900007214255458)),
  list(alg = "KahanLi8", case = "P1", dt = 0.1, nf = 181, nf2 = 181, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083573835521446, -0.3308513228849953, 0.5068323377044913, -0.14530826292528315),
       tq = 0.55, uq = c(-0.6666689133223765, -0.36142243290894704, 0.9488327957519066, 0.013183125056464285)),
  list(alg = "KahanLi8", case = "P5", dt = 0.1, nf = 181, nf2 = 181, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.1989070453149361, 0.04202789783687951, 1.0718828511312744, -0.033535781409556445),
       tq = 0.55, uq = c(-0.0031670298399488163, -0.26703679563836996, 1.125883253136144, 0.025714138798707615)),
  list(alg = "SofSpa10", case = "P1", dt = 0.1, nf = 361, nf2 = 361, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2083573835402257, -0.33085132288498936, 0.5068323377178426, -0.14530826292528828),
       tq = 0.55, uq = c(-0.6666689133037587, -0.3614224329089403, 0.9488327957583587, 0.013183125056461522)),
  list(alg = "SofSpa10", case = "P5", dt = 0.1, nf = 361, nf2 = 361, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-0.18753163127669986, 0.05340331187510565, 1.0780972681433865, -0.027321364397443374),
       tq = 0.55, uq = c(0.0038989084999475777, -0.25997085729847896, 1.1278769013752488, 0.027707787037817434)),
  list(alg = "SymplecticEuler", case = "P3", dt = 0.1, nf = 12, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.2148747686554624, -0.29085229164632714, 0.5605553803330047, -0.12357872437250424),
       tq = 0.55, uq = c(-0.7158280338137859, -0.33645098835648535, 0.9906731532874972, 0.022111493916745807)),
  list(alg = "SymplecticEuler", case = "P4", dt = 0.1, nf = 12, nf2 = 11, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.3813310897029345, -0.4067635634550078, 0.7862839518284969, 0.018002933415157386),
       tq = 0.55, uq = c(-0.66574093663896, -0.38275232056057934, 1.0050648498710995, 0.1069442934817494)),
  list(alg = "LeapfrogDriftKickDrift", case = "P3", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.1654921944366305, -0.30319217540495574, 0.5143180475600834, -0.12991373375723014),
       tq = 0.55, uq = c(-0.6621338035910229, -0.3437434964735416, 0.943810285254923, 0.018114533949779573)),
  list(alg = "LeapfrogDriftKickDrift", case = "P4", dt = 0.1, nf = 21, nf2 = 21, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.3471161526589976, -0.41393339258288603, 0.7499768224395109, 0.013424339259257309),
       tq = 0.55, uq = c(-0.6386940610789924, -0.38680499634168397, 0.9784964754848218, 0.10420872591649803)),
  list(alg = "McAte3", case = "P3", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.164822752812821, -0.3030911003434694, 0.514120977048355, -0.1302737868546605),
       tq = 0.55, uq = c(-0.6635430745843806, -0.3443036258939032, 0.9461912518263154, 0.017869921709111357)),
  list(alg = "McAte3", case = "P4", dt = 0.1, nf = 41, nf2 = 31, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.3430593630557375, -0.41381481480509535, 0.7510065658290961, 0.013378001352031006),
       tq = 0.55, uq = c(-0.6368187257138789, -0.38697379696956613, 0.98015514001294, 0.10424373710369604)),
  list(alg = "KahanLi6", case = "P3", dt = 0.1, nf = 101, nf2 = 101, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.1696978836834073, -0.3029868952323758, 0.5099832081453511, -0.13040800602328428),
       tq = 0.55, uq = c(-0.6691811329168299, -0.3445233496313977, 0.9445299609491019, 0.017780627059118505)),
  list(alg = "KahanLi6", case = "P4", dt = 0.1, nf = 101, nf2 = 101, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.343064350589938, -0.4138139735539972, 0.7510168212599124, 0.013378401687896644),
       tq = 0.55, uq = c(-0.6368289512336333, -0.3869732424784287, 0.9801658947584834, 0.1042439557654839)),
  list(alg = "SofSpa10", case = "P3", dt = 0.1, nf = 361, nf2 = 361, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.1682165396848296, -0.30301161771458096, 0.5112979873065412, -0.13036234961975707),
       tq = 0.55, uq = c(-0.6674620498603556, -0.3444508347903992, 0.9450814402337916, 0.017810749917166444)),
  list(alg = "SofSpa10", case = "P4", dt = 0.1, nf = 361, nf2 = 361, naccept = 10, nreject = 0, nt = 11,
       t2 = 0.1, u = c(-1.3430643397726603, -0.4138139735458216, 0.7510168253188669, 0.013378401694405096),
       tq = 0.55, uq = c(-0.6368289402911871, -0.3869732424700585, 0.980165896358728, 0.10424395576832367))
)
symp_vl_055 <- c(-0.6610044766075452, -0.3608622724698778, 0.9470714395476934, 0.012866281350848652)
symp_mcate5_saveat <- list(
  c(0.18164179125267615, -0.3273583626181145, 1.0515241858305258, 0.1528512956895723),
  c(-0.6666679417342837, -0.3614225355961255, 0.9488315380617033, 0.01318295541003689))

symp_check <- function(r) {
  prob <- symp_test_probs[[r$case]]()
  sol <- solve(prob, get(r$alg)(), dt = r$dt, adaptive = FALSE)
  lab <- paste(r$alg, r$case)
  st <- sol$stats
  expect_identical(c(st$nf, st$nf2, st$naccept, st$nreject, length(sol$t)),
                   as.integer(c(r$nf, r$nf2, r$naccept, r$nreject, r$nt)), label = lab)
  expect_equal(unname(sol$u[[length(sol$t)]]), r$u, tolerance = 1e-12, label = lab)
  if (!is.null(r$tq)) expect_equal(unname(sol(r$tq)), r$uq, tolerance = 1e-12, label = lab)
}

test_that("the symplectic methods match OrdinaryDiffEqSymplecticRK", {
  for (r in symp_refs) symp_check(r)
})

test_that("VerletLeapfrog interpolates linearly by default and Hermite at saveat", {
  # tools/rkn-probe.jl: VerletLeapfrog dense = false; McAte5 with saveat
  sol <- solve(symp_test_probs$P1(), VerletLeapfrog(), dt = 0.1)
  expect_false(sol$dense)
  expect_equal(unname(sol(0.55)), symp_vl_055, tolerance = 1e-12)
  expect_true(solve(symp_test_probs$P1(), VerletLeapfrog(), dt = 0.1, dense = TRUE)$dense)
  sol <- solve(symp_test_probs$P1(), McAte5(), dt = 0.1, saveat = c(0.15, 0.55))
  expect_equal(sol$t, c(0.15, 0.55))
  expect_equal(unname(sol$u[[1]]), symp_mcate5_saveat[[1]], tolerance = 1e-12)
  expect_equal(unname(sol$u[[2]]), symp_mcate5_saveat[[2]], tolerance = 1e-12)
})

test_that("VelocityVerlet needs du/dt = v and partitioned methods a dynamical problem", {
  expect_error(solve(symp_test_probs$P4(), VelocityVerlet(), dt = 0.1), "not applicable if f2")
  expect_error(solve(ODEProblem(function(u, p, t) -u, 1, c(0, 1)), McAte3(), dt = 0.1),
               "designed for partitioned ODE")
  expect_error(solve(symp_test_probs$P1(), McAte3()), "fixed step method")
})

test_that("SecondOrderODEProblem works with the symplectic methods", {
  prob <- SecondOrderODEProblem(function(v, u, p, t) -u - u * u * u, c(0.5, -0.3), c(1, 0.2), c(0, 1))
  sol <- solve(prob, KahanLi8(), dt = 0.1)
  ref <- solve(symp_test_probs$P1(), KahanLi8(), dt = 0.1)
  expect_equal(sol$u, ref$u)
})

test_that("after an event the partitioned methods restart from the new state", {
  # Upstream refreshes only fsalfirst after an event, so VelocityVerlet,
  # VerletLeapfrog and the compositions step from the derivatives before it;
  # here a run with an event equals a restart from the modified state.
  # A kick of the position (index 3) changes f1, one of the velocity
  # (index 1) changes f2.
  for (ki in c(1L, 3L)) {
    kick <- function(integ) integ$u[ki] <- integ$u[ki] + 0.1
    for (alg in list(SymplecticEuler(), VelocityVerlet(), VerletLeapfrog(), McAte3(), DPRKN6(), Nystrom4(),
                     IRKN4(), NewmarkBeta())) {
      cb <- DiscreteCallback(function(u, t, integ) t == 0.5, kick)
      sol <- solve(symp_test_probs$P1(), alg, dt = 0.1, adaptive = FALSE, callback = cb, tstops = 0.5)
      first <- solve(remake(symp_test_probs$P1(), tspan = c(0, 0.5)), alg, dt = 0.1, adaptive = FALSE)
      z <- first$u[[length(first$t)]]
      z[ki] <- z[ki] + 0.1
      second <- solve(DynamicalODEProblem(function(v, u, p, t) -u - u * u * u, function(v, u, p, t) v,
                                          z[1:2], z[3:4], c(0.5, 1)), alg, dt = 0.1, adaptive = FALSE)
      expect_equal(sol$u[[length(sol$t)]], second$u[[length(second$t)]], tolerance = 1e-14,
                   label = paste(alg$name, ki))
    }
  }
})
