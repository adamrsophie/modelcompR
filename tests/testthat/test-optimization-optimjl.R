# Optim.jl methods driven by OptimizationOptimJL. The reference rows come from
# tools/optimjl-probe.jl (Optim 2.3.2, LineSearches 7.8.2, OptimizationOptimJL
# 0.4.23), with the analytic gradient and Hessian given to both sides:
# problem, method, maxiters (0: the default), retcode, u, objective,
# iterations, fevals, gevals and hevals.
# long reference comparisons: run on CI (NOT_CRAN = true), not on CRAN
skip_on_cran()

ojl_ref <- c(
  "quad\tGradientDescent\t3\tMaxIters\t0.569274095941856\t-0.7505857983942367\t-0.697660410927028\t3\t5\t5\t0",
  "quad\tGradientDescent\t0\tSuccess\t0.5999999964948312\t-0.7999999943285176\t-0.6999999999999998\t33\t35\t35\t0",
  "quad\tConjugateGradient\t3\tSuccess\t0.6000000000000205\t-0.8000000000000245\t-0.7\t2\t5\t3\t0",
  "quad\tConjugateGradient\t0\tSuccess\t0.6000000000000205\t-0.8000000000000245\t-0.7\t2\t5\t3\t0",
  "quad\tAcceleratedGradientDescent\t3\tMaxIters\t0.5840715228328497\t-0.7775362450892404\t-0.6994726185451973\t3\t7\t7\t0",
  "quad\tAcceleratedGradientDescent\t0\tSuccess\t0.5999999965363522\t-0.7999999943957\t-0.7\t40\t81\t81\t0",
  "quad\tMomentumGradientDescent\t3\tMaxIters\t0.5704279683160138\t-0.7515264219904845\t-0.697772016832811\t3\t8\t8\t0",
  "quad\tMomentumGradientDescent\t0\tSuccess\t0.5999999959360753\t-0.7999999934244318\t-0.7000000000000001\t32\t66\t66\t0",
  "quad\tBFGS\t3\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "quad\tBFGS\t0\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "quad\tBFGSstep\t3\tMaxIters\t0.5264078326004735\t-0.7390708000522995\t-0.6926478338219637\t3\t7\t7\t0",
  "quad\tBFGSstep\t0\tSuccess\t0.6000000006224367\t-0.7999999998551609\t-0.7\t8\t12\t12\t0",
  "quad\tLBFGS\t3\tSuccess\t0.5999999999999999\t-0.7999999999999999\t-0.7\t3\t5\t5\t0",
  "quad\tLBFGS\t0\tSuccess\t0.5999999999999999\t-0.7999999999999999\t-0.7\t3\t5\t5\t0",
  "quad\tLBFGS3\t3\tSuccess\t0.5999999999999999\t-0.7999999999999999\t-0.7\t3\t5\t5\t0",
  "quad\tLBFGS3\t0\tSuccess\t0.5999999999999999\t-0.7999999999999999\t-0.7\t3\t5\t5\t0",
  "quad\tNewton\t3\tSuccess\t0.6\t-0.8\t-0.7\t1\t2\t2\t2",
  "quad\tNewton\t0\tSuccess\t0.6\t-0.8\t-0.7\t1\t2\t2\t2",
  "quad\tNewtonTrustRegion\t3\tSuccess\t0.6000000000000001\t-0.8\t-0.7\t1\t2\t2\t2",
  "quad\tNewtonTrustRegion\t0\tSuccess\t0.6000000000000001\t-0.8\t-0.7\t1\t2\t2\t2",
  "quad\tNGMRES\t3\tSuccess\t0.6000000000000001\t-0.8000000000000002\t-0.7000000000000001\t3\t7\t7\t0",
  "quad\tNGMRES\t0\tSuccess\t0.6000000000000001\t-0.8000000000000002\t-0.7000000000000001\t3\t7\t7\t0",
  "quad\tOACCEL\t3\tSuccess\t0.6\t-0.8\t-0.7\t3\t7\t7\t0",
  "quad\tOACCEL\t0\tSuccess\t0.6\t-0.8\t-0.7\t3\t7\t7\t0",
  "quad\tGDback\t3\tMaxIters\t0.569274095941856\t-0.7505857983942367\t-0.697660410927028\t3\t5\t4\t0",
  "quad\tGDback\t0\tSuccess\t0.59999999497157\t-0.7999999918638294\t-0.7\t34\t37\t35\t0",
  "quad\tBFGSback2\t3\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t4\t0",
  "quad\tBFGSback2\t0\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t4\t0",
  "quad\tCGquad\t3\tMaxIters\t0.5568238869071211\t-0.7892059717267802\t-0.6975532680263755\t3\t5\t5\t0",
  "quad\tCGquad\t0\tSuccess\t0.6\t-0.7999999999999999\t-0.7\t6\t10\t10\t0",
  "quad\tMGDcc\t3\tMaxIters\t0.6167381703817816\t-0.8371173203348146\t-0.6988233310414477\t3\t9\t9\t0",
  "quad\tMGDcc\t0\tSuccess\t0.5999999999639405\t-0.7999999961753544\t-0.7\t24\t63\t63\t0",
  "quad\tGDstatic\t3\tMaxIters\t0.6979089884999999\t-0.402692508\t-0.48886752709183856\t3\t4\t4\t0",
  "quad\tGDstatic\t0\tMaxIters\t0.5687589289155818\t-0.7344732655285944\t-0.6962893656553123\t1000\t1001\t1001\t0",
  "quad\tNGMRESlbfgs\t3\tSuccess\t0.599999999999937\t-0.7999999999999248\t-0.7\t2\t5\t5\t0",
  "quad\tNGMRESlbfgs\t0\tSuccess\t0.599999999999937\t-0.7999999999999248\t-0.7\t2\t5\t5\t0",
  "rosen\tGradientDescent\t3\tMaxIters\t-1.032037865062163\t1.0676086750652847\t4.129806149369229\t3\t6\t6\t0",
  "rosen\tGradientDescent\t0\tMaxIters\t0.9804139031961087\t0.9610281637199409\t0.00038697353234713\t1000\t1431\t1431\t0",
  "rosen\tConjugateGradient\t3\tMaxIters\t-1.0227784515001652\t1.0623903039062883\t4.118249095356918\t3\t8\t6\t0",
  "rosen\tConjugateGradient\t0\tSuccess\t1.0000000011561798\t1.0000000023119946\t1.3367650983421405e-18\t27\t68\t47\t0",
  "rosen\tAcceleratedGradientDescent\t3\tMaxIters\t-1.027690963082564\t1.06867355521915\t4.127217802511311\t3\t8\t8\t0",
  "rosen\tAcceleratedGradientDescent\t0\tMaxIters\t3.9891881608603135\t15.917316083191674\t8.936610351076173\t1000\t2125\t2125\t0",
  "rosen\tMomentumGradientDescent\t3\tMaxIters\t-1.031385663507855\t1.0678038672533108\t4.128165923634761\t3\t8\t8\t0",
  "rosen\tMomentumGradientDescent\t0\tMaxIters\t0.9846337496792824\t0.9694617717335943\t0.0002362967850914613\t1000\t2293\t2293\t0",
  "rosen\tBFGS\t3\tMaxIters\t-0.941840255659282\t0.8593425997244014\t3.84758601009649\t3\t9\t9\t0",
  "rosen\tBFGS\t0\tSuccess\t1.0000000000000007\t1.0000000000000016\t5.374114916818143e-30\t41\t57\t57\t0",
  "rosen\tBFGSstep\t3\tMaxIters\t-1.0327206498143497\t1.068550689941779\t4.132368890088708\t3\t5\t5\t0",
  "rosen\tBFGSstep\t0\tSuccess\t1.0000000004389917\t1.000000000900043\t2.413757418579879e-19\t44\t64\t64\t0",
  "rosen\tLBFGS\t3\tMaxIters\t-1.029534506181162\t1.0677753099802143\t4.125147483929388\t3\t6\t6\t0",
  "rosen\tLBFGS\t0\tSuccess\t0.9999999998934164\t0.9999999997940162\t1.652030315296824e-20\t44\t64\t64\t0",
  "rosen\tLBFGS3\t3\tMaxIters\t-1.029534506181162\t1.0677753099802143\t4.125147483929388\t3\t6\t6\t0",
  "rosen\tLBFGS3\t0\tSuccess\t1.0000000000225073\t1.0000000000482614\t1.5607094775514643e-21\t40\t64\t64\t0",
  "rosen\tNewton\t3\tMaxIters\t-0.7832924977869852\t0.5997008876960553\t3.1993039948895428\t3\t11\t11\t4",
  "rosen\tNewton\t0\tSuccess\t0.9999999999999981\t0.999999999999995\t1.5270621491848618e-28\t24\t49\t49\t25",
  "rosen\tNewtonTrustRegion\t3\tMaxIters\t-1.0745376631668053\t1.151871184291669\t4.304468278807907\t3\t4\t4\t3",
  "rosen\tNewtonTrustRegion\t0\tSuccess\t0.9999999999999304\t0.9999999999998552\t7.927176954904439e-27\t26\t27\t27\t24",
  "rosen\tNGMRES\t3\tMaxIters\t-1.0348770294809473\t1.067373658490525\t4.142018227641194\t3\t6\t6\t0",
  "rosen\tNGMRES\t0\tMaxIters\t-0.985314660026175\t0.978855802875389\t3.9478916288105332\t1000\t1004\t1004\t0",
  "rosen\tOACCEL\t3\tMaxIters\t-1.0348845627481584\t1.0673718748661059\t4.142061415838306\t3\t6\t6\t0",
  "rosen\tOACCEL\t0\tSuccess\t0.9999999998855726\t0.9999999997808762\t2.2562627851801628e-20\t86\t231\t231\t0",
  "rosen\tGDback\t3\tMaxIters\t1.4533385266242784\t2.113604180152446\t0.20571499851807307\t3\t10\t4\t0",
  "rosen\tGDback\t0\tMaxIters\t1.4057083161565604\t1.9772969981891078\t0.16476336671366063\t1000\t1007\t1001\t0",
  "rosen\tBFGSback2\t3\tMaxIters\t-0.5458601897188354\t0.29483563228089027\t2.390661985918743\t3\t10\t4\t0",
  "rosen\tBFGSback2\t0\tSuccess\t1.0000000000310387\t1.0000000000669311\t3.319216760314283e-21\t36\t50\t37\t0",
  "rosen\tCGquad\t3\tMaxIters\t-1.0449931560145909\t1.0631588529261562\t4.265239893699355\t3\t9\t9\t0",
  "rosen\tCGquad\t0\tSuccess\t0.9999999998849469\t0.9999999997543485\t3.7402979046060516e-20\t65\t123\t123\t0",
  "rosen\tMGDcc\t3\tMaxIters\t-1.0281702876549232\t1.0674587805909648\t4.124134535199662\t3\t11\t11\t0",
  "rosen\tMGDcc\t0\tMaxIters\t0.908364050535666\t0.8246293233134251\t0.008421741394012576\t1000\t2027\t2027\t0",
  "rosen\tGDstatic\t3\tMaxIters\t-1.0268830682337082\t1.0624243118951573\t4.114552150338377\t3\t4\t4\t0",
  "rosen\tGDstatic\t0\tMaxIters\t0.3272627747528461\t0.10401280036956503\t0.45352902482768304\t1000\t1001\t1001\t0",
  "rosen\tNGMRESlbfgs\t3\tMaxIters\t-1.0281815447183373\t1.0664540281806858\t4.122163314462505\t3\t7\t7\t0",
  "rosen\tNGMRESlbfgs\t0\tSuccess\t1.0000000002217213\t1.000000000447674\t5.09507131415956e-20\t57\t185\t185\t0",
  "rosenbox\tFminboxLBFGS\t3\tSuccess\t0.49999999347278695\t0.2499999934569644\t0.25000000652721305\t3\t50\t50\t0",
  "rosenbox\tFminboxLBFGS\t0\tSuccess\t0.49999999347278695\t0.2499999934569644\t0.25000000652721305\t3\t50\t50\t0",
  "rosenbox\tFminboxGD\t3\tSuccess\t0.4999999934684133\t0.2499999935008139\t0.25000000653158677\t3\t2477\t2477\t0",
  "rosenbox\tFminboxGD\t0\tSuccess\t0.4999999934684133\t0.2499999935008139\t0.25000000653158677\t3\t2477\t2477\t0",
  "rosenbox\tFminboxCG\t3\tFailure\t0.4999999934736008\t0.2499999934155579\t0.25000000652639925\t3\t113\t81\t0",
  "rosenbox\tFminboxCG\t0\tSuccess\t0.49999999999999994\t0.2500000000645566\t0.25\t7\t473\t403\t0",
  "rosenbox\tFminboxBFGS\t3\tSuccess\t0.499999993472787\t0.24999999345696358\t0.25000000652721305\t3\t59\t59\t0",
  "rosenbox\tFminboxBFGS\t0\tSuccess\t0.499999993472787\t0.24999999345696358\t0.25000000652721305\t3\t59\t59\t0",
  "rosenbox\tGradientDescent\t3\tSuccess\t0.4999999934684133\t0.2499999935008139\t0.25000000653158677\t3\t2477\t2477\t0",
  "rosenbox\tGradientDescent\t0\tSuccess\t0.4999999934684133\t0.2499999935008139\t0.25000000653158677\t3\t2477\t2477\t0",
  "quad\tNelderMead\t3\tMaxIters\t0.5476562500000001\t-0.6515625000000002\t-0.6816262817382813\t3\t9\t0\t0",
  "quad\tNelderMead\t0\tSuccess\t0.5999995883209037\t-0.8000563074636131\t-0.6999999968060348\t30\t62\t0\t0",
  "rosen\tNelderMead\t3\tMaxIters\t-1.18203125\t1.369140625\t4.839981309212745\t3\t10\t0\t0",
  "rosen\tNelderMead\t0\tSuccess\t1.0000117499532974\t1.0000167773369513\t4.657541291131408e-9\t78\t149\t0\t0",
  "rosenbox\tFminboxNM\t3\tSuccess\t0.4935782166380482\t0.2436023109623764\t0.25646305205852393\t2\t84\t2\t0",
  "rosenbox\tFminboxNM\t0\tSuccess\t0.4935782166380482\t0.2436023109623764\t0.25646305205852393\t2\t84\t2\t0",
  "rosenbox\tIPNewton\t3\tMaxIters\t0.3322609023800811\t0.08201985982066987\t0.5264034547531197\t3\t6\t6\t4",
  "rosenbox\tIPNewton\t0\tSuccess\t0.4999999999999917\t0.24999999999999165\t0.2500000000000083\t33\t50\t50\t33",
  "rosenbox\tIPNewtonmu\t3\tMaxIters\t0.33062373433986314\t0.08035481966897612\t0.5319167253530344\t3\t6\t6\t4",
  "rosenbox\tIPNewtonmu\t0\tSuccess\t0.49999999999999994\t0.24999999999999994\t0.25\t40\t56\t56\t40",
  "rosencons\tIPNewton\t3\tMaxIters\t0.24612346217209752\t0.048828846532581\t0.5821311781559653\t3\t6\t6\t4",
  "rosencons\tIPNewton\t0\tSuccess\t0.724701840447071\t0.5242206047577945\t0.07588358475785828\t37\t56\t56\t37",
  "rosencons\tIPNewtonmu\t3\tMaxIters\t0.2409443338031399\t0.047248900174901354\t0.5878408942893497\t3\t6\t6\t4",
  "rosencons\tIPNewtonmu\t0\tSuccess\t0.7247018404470709\t0.5242206047577944\t0.07588358475785831\t36\t64\t64\t36",
  "quadeq\tIPNewton\t3\tMaxIters\t0.8166666666666667\t-0.3416666666666666\t-0.32020833333333315\t3\t7\t7\t4",
  "quadeq\tIPNewton\t0\tSuccess\t0.8333333333333334\t-0.33333333333333337\t-0.29166666666666674\t25\t46\t46\t25",
  "quadeq\tIPNewtonmu\t3\tMaxIters\t0.8166666666666667\t-0.3416666666666667\t-0.3202083333333332\t3\t7\t7\t4",
  "quadeq\tIPNewtonmu\t0\tSuccess\t0.8333333333333334\t-0.33333333333333337\t-0.29166666666666674\t27\t50\t50\t26",
  "rosenboth\tIPNewton\t3\tMaxIters\t0.15261346180426866\t0.0290688973830816\t0.7214025066339976\t3\t5\t5\t4",
  "rosenboth\tIPNewton\t0\tSuccess\t0.49999999999999994\t0.24999999999999997\t0.25\t45\t63\t63\t45",
  "rosenboth\tIPNewtonmu\t3\tMaxIters\t0.2071049984270005\t0.042731212625433945\t0.6286850842481019\t3\t5\t5\t4",
  "rosenboth\tIPNewtonmu\t0\tSuccess\t0.5\t0.249999999981331\t0.25\t44\t57\t57\t45",
  "quad\tBFGSmt\t3\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "quad\tBFGSmt\t0\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "quad\tGDmt\t3\tMaxIters\t0.569274095941856\t-0.7505857983942367\t-0.697660410927028\t3\t5\t5\t0",
  "quad\tGDmt\t0\tSuccess\t0.5999999957351806\t-0.7999999930993772\t-0.7\t33\t37\t37\t0",
  "quad\tCGsw\t3\tSuccess\t0.6000000000000205\t-0.8000000000000245\t-0.7\t2\t5\t3\t0",
  "quad\tCGsw\t0\tSuccess\t0.6000000000000205\t-0.8000000000000245\t-0.7\t2\t5\t3\t0",
  "quad\tLBFGSsw\t3\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "quad\tLBFGSsw\t0\tSuccess\t0.6\t-0.8\t-0.7\t3\t5\t5\t0",
  "rosen\tBFGSmt\t3\tMaxIters\t-0.7124017331587208\t0.5083557691030471\t2.932390178415022\t3\t12\t12\t0",
  "rosen\tBFGSmt\t0\tSuccess\t0.9999999999995596\t0.9999999999990832\t3.2336750801143314e-25\t34\t51\t51\t0",
  "rosen\tGDmt\t3\tMaxIters\t-1.025575942676001\t1.0606058482324539\t4.110701607455092\t3\t8\t8\t0",
  "rosen\tGDmt\t0\tMaxIters\t0.9776617966473211\t0.9556773672902635\t0.0005011042525932075\t1000\t1457\t1457\t0",
  "rosen\tCGsw\t3\tMaxIters\t-1.025291101941779\t1.0607415264247253\t4.110866483482529\t3\t11\t9\t0",
  "rosen\tCGsw\t0\tSuccess\t0.9999999994824655\t0.9999999989615072\t2.690142451733964e-19\t29\t78\t56\t0",
  "rosen\tLBFGSsw\t3\tMaxIters\t-1.0082497284807113\t1.006648605514182\t4.0429054484447935\t3\t13\t13\t0",
  "rosen\tLBFGSsw\t0\tSuccess\t1.0000000000629952\t1.0000000001314304\t6.927851757796163e-21\t34\t51\t51\t0"
)

ojl_quad <- function(x, p) 1.5 * x[1] * x[1] + x[1] * x[2] + x[2] * x[2] - x[1] + x[2]
ojl_quadg <- function(x, p) c(3 * x[1] + x[2] - 1, x[1] + 2 * x[2] + 1)
ojl_quadh <- function(x, p) matrix(c(3, 1, 1, 2), 2, 2)
ojl_rosen <- function(x, p) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
ojl_rosg <- function(x, p) c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2), 200 * (x[2] - x[1]^2))
ojl_rosh <- function(x, p) matrix(c(2 - 400 * (x[2] - x[1]^2) + 800 * x[1]^2, -400 * x[1], -400 * x[1], 200), 2, 2)
ojl_c1 <- function(x, p) x[1]^2 + x[2]^2
ojl_c1j <- function(x, p) matrix(c(2 * x[1], 2 * x[2]), 1, 2)
ojl_c1h <- function(x, p) list(matrix(c(2, 0, 0, 2), 2, 2))
ojl_probs <- list(
  quad = OptimizationProblem(ojl_quad, c(0.7, -0.4), 1, grad = ojl_quadg, hess = ojl_quadh),
  rosen = OptimizationProblem(ojl_rosen, c(-1.2, 1), 1, grad = ojl_rosg, hess = ojl_rosh),
  rosenbox = OptimizationProblem(ojl_rosen, c(0.1, 0.2), 1, grad = ojl_rosg, hess = ojl_rosh, lb = c(-0.5, -0.5),
                                 ub = c(0.5, 0.8)),
  rosencons = OptimizationProblem(ojl_rosen, c(0.1, 0.2), 1, grad = ojl_rosg, hess = ojl_rosh, cons = ojl_c1,
                                  cons_j = ojl_c1j, cons_h = ojl_c1h, lcons = -Inf, ucons = 0.8),
  quadeq = OptimizationProblem(ojl_quad, c(0.7, -0.4), 1, grad = ojl_quadg, hess = ojl_quadh,
                               cons = function(x, p) x[1] + x[2], cons_j = function(x, p) matrix(c(1, 1), 1, 2),
                               cons_h = function(x, p) list(matrix(0, 2, 2)), lcons = 0.5, ucons = 0.5),
  rosenboth = OptimizationProblem(ojl_rosen, c(0.1, 0.2), 1, grad = ojl_rosg, hess = ojl_rosh, lb = c(-0.5, -0.5),
                                  ub = c(0.5, 0.8), cons = ojl_c1, cons_j = ojl_c1j, cons_h = ojl_c1h, lcons = 0.01,
                                  ucons = 0.5))
ojl_LS <- Optim$LineSearches
ojl_methods <- list(
  GradientDescent = function() GradientDescent(), ConjugateGradient = function() ConjugateGradient(),
  AcceleratedGradientDescent = function() AcceleratedGradientDescent(),
  MomentumGradientDescent = function() MomentumGradientDescent(), BFGS = function() Optim$BFGS(),
  BFGSstep = function() Optim$BFGS(initial_stepnorm = 0.01), LBFGS = function() Optim$LBFGS(),
  LBFGS3 = function() Optim$LBFGS(m = 3), Newton = function() Newton(),
  NewtonTrustRegion = function() NewtonTrustRegion(), NGMRES = function() NGMRES(), OACCEL = function() OACCEL(),
  GDback = function() GradientDescent(linesearch = ojl_LS$BackTracking()),
  BFGSback2 = function() Optim$BFGS(linesearch = ojl_LS$BackTracking(order = 2)),
  CGquad = function() ConjugateGradient(alphaguess = ojl_LS$InitialQuadratic()),
  MGDcc = function() MomentumGradientDescent(alphaguess = ojl_LS$InitialConstantChange()),
  GDstatic = function() GradientDescent(alphaguess = ojl_LS$InitialStatic(alpha = 0.001), linesearch = ojl_LS$Static()),
  NGMRESlbfgs = function() NGMRES(nlprecon = Optim$LBFGS()),
  BFGSmt = function() Optim$BFGS(linesearch = ojl_LS$MoreThuente()),
  GDmt = function() GradientDescent(linesearch = ojl_LS$MoreThuente()),
  CGsw = function() ConjugateGradient(linesearch = ojl_LS$StrongWolfe()),
  LBFGSsw = function() Optim$LBFGS(linesearch = ojl_LS$StrongWolfe()),
  FminboxLBFGS = function() Fminbox(), FminboxGD = function() Fminbox(GradientDescent()),
  FminboxCG = function() Fminbox(ConjugateGradient()), FminboxBFGS = function() Fminbox(Optim$BFGS()),
  NelderMead = function() Optim$NelderMead(), IPNewton = function() IPNewton(),
  IPNewtonmu = function() IPNewton(mu0 = 0.1), FminboxNM = function() Fminbox(Optim$NelderMead()))

test_that("the Optim methods follow the iterates of Optim.jl", {
  for (row in strsplit(ojl_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[1:3], collapse = " ")
    prob <- ojl_probs[[row[1]]]
    s <- suppressWarnings(if (n == 0L) solve(prob, ojl_methods[[row[2]]]())
                          else solve(prob, ojl_methods[[row[2]]](), maxiters = n))
    ju <- as.numeric(row[5:6])
    jo <- as.numeric(row[7])
    rel <- max(abs(unname(s$u) - ju) / pmax(1, abs(ju)), abs(s$objective - jo) / max(1, abs(jo)))
    # Newton and NewtonTrustRegion factor the Hessian with LAPACK, whose
    # rounding differs from Julia's in the last bits
    expect_lt(rel, ref_tol(1e-10, 1e-6), label = info)
    expect_identical(s$retcode, row[4], info = info)
    st <- s$stats
    if (ref_platform()) {
      expect_identical(as.integer(c(st$iterations, st$fevals, st$gevals, st$hevals)), as.integer(row[8:11]),
                       info = info)
    }
  }
})

test_that("the randomized Optim methods converge", {
  set.seed(42)
  box <- ojl_probs$rosenbox
  # the minimizer on the box is (0.5, 0.25), with the value 0.25
  s <- solve(box, ParticleSwarm(lower = box$lb, upper = box$ub, n_particles = 20), maxiters = 100)
  expect_lt(abs(s$objective - 0.25), 1e-4)
  expect_true(all(s$u >= box$lb & s$u <= box$ub))
  s <- suppressWarnings(solve(box, SAMIN(), maxiters = 4000))
  expect_lt(abs(s$objective - 0.25), 1e-3)
  s <- solve(ojl_probs$quad, SimulatedAnnealing(), maxiters = 3000)
  expect_lt(s$objective, ojl_quad(c(0.7, -0.4), 1))
  expect_lt(s$objective, -0.69)
})

test_that("the Optim methods check their problems", {
  expect_error(solve(ojl_probs$rosencons, GradientDescent()), "IPNewton")
  expect_error(solve(ojl_probs$rosen, SAMIN()), "bounds")
  expect_error(solve(ojl_probs$rosenbox, SimulatedAnnealing()), "bounds")
  expect_error(Fminbox(Newton()), "Newton is not supported")
  expect_error(NewtonTrustRegion(delta_hat = -1), "positive")
  expect_error(IPNewton(mu0 = "large"), "auto")
  expect_error(GradientDescent(manifold = "Sphere"), "must be NULL or a manifold")
  expect_error(Optim$Stiefel("QR"), "SVD")
  expect_error(NGMRES(manifold = Optim$Sphere(), nlprecon = GradientDescent()), "manifold == nlprecon.manifold")
  expect_error(solve(ojl_probs$quad, Optim$KrylovTrustRegion(eta = 0.3)), "AssertionError")
  expect_error(Optim$LineSearches$StrongWolfe(c_1 = 0.95), "0 < c_1 < c_2 < 1")
  expect_error(Optim$LineSearches$StrongWolfe(rho = 1), "bracket growth factor")
})

test_that("BFGS, LBFGS and NelderMead of Optim are in the Optim list", {
  expect_true(all(c("BFGS", "LBFGS", "NelderMead", "LineSearches") %in% names(Optim)))
  expect_s3_class(Optim$BFGS(), "OptimizationAlgorithm")
  expect_identical(Optim$LBFGS(m = 3)$name, "LBFGS")
})

ojl_ext_ref <- c(
  "quad\tAdam\t3\tMaxIters\t0.6981022068658238\t-0.40189676104055777\t-0.4880229403397301\t3\t4\t4\t0",
  "quad\tAdam\t0\tMaxIters\t0.6116074653371643\t-0.4951181626766023\t-0.6033060600343167\t1000\t1001\t1001\t0",
  "quad\tAdam2\t3\tMaxIters\t0.48531619862956865\t-0.6053623175155592\t-0.6647094004291738\t3\t4\t4\t0",
  "quad\tAdam2\t0\tSuccess\t0.5999999995053055\t-0.7999999985483311\t-0.7000000000000001\t288\t289\t289\t0",
  "quad\tAdamsched\t3\tMaxIters\t0.5253653690571295\t-0.568312892817224\t-0.6552574739001169\t3\t4\t4\t0",
  "quad\tAdamsched\t0\tMaxIters\t0.5999764420213131\t-0.7999697750388424\t-0.6999999989660234\t1000\t1001\t1001\t0",
  "quad\tAdaMax\t3\tMaxIters\t0.6623933554440168\t-0.43774203644101545\t-0.5403272817499909\t3\t4\t4\t0",
  "quad\tAdaMax\t0\tMaxIters\t0.5996267743996121\t-0.7993050064479421\t-0.6999995674273253\t1000\t1001\t1001\t0",
  "quad\tKrylovTrustRegion\t3\tSuccess\t0.6000000000000001\t-0.7999999999999999\t-0.7000000000000001\t1\t2\t2\t0",
  "quad\tKrylovTrustRegion\t0\tSuccess\t0.6000000000000001\t-0.7999999999999999\t-0.7000000000000001\t1\t2\t2\t0",
  "rosen\tAdam\t3\tMaxIters\t-1.198101045702814\t1.001898667147572\t23.62798722298106\t3\t4\t4\t0",
  "rosen\tAdam\t0\tMaxIters\t-1.1172623976413614\t1.0835440475362603\t7.196437467214458\t1000\t1001\t1001\t0",
  "rosen\tAdam2\t3\tMaxIters\t-0.9832164401465344\t1.2158829354240501\t10.141634972251964\t3\t4\t4\t0",
  "rosen\tAdam2\t0\tMaxIters\t-0.5591667267355083\t0.3182964392059543\t2.4341694581504485\t1000\t1001\t1001\t0",
  "rosen\tAdamsched\t3\tMaxIters\t-1.023135199045846\t1.1766311117398558\t5.778541461018588\t3\t4\t4\t0",
  "rosen\tAdamsched\t0\tMaxIters\t-1.0406276061993276\t1.089735044603734\t4.1688248651756785\t1000\t1001\t1001\t0",
  "rosen\tAdaMax\t3\tMaxIters\t-1.1625338497884354\t1.037495600493354\t14.535483931403615\t3\t4\t4\t0",
  "rosen\tAdaMax\t0\tMaxIters\t-1.042397293428487\t1.0928579565020409\t4.17531277823565\t1000\t1001\t1001\t0",
  "rosen\tKrylovTrustRegion\t3\tMaxIters\t-1.0740775360522827\t1.1520743336965065\t4.302043556876904\t3\t4\t3\t0",
  "rosen\tKrylovTrustRegion\t0\tMaxIters\t0.999999968268468\t0.9999999364101126\t1.0084985642282296e-15\t1000\t26\t23\t0",
  "rosen\tGDsphere\t3\tMaxIters\t-0.783930185935964\t0.6208489861298037\t3.1863789955259727\t3\t10\t10\t0",
  "rosen\tGDsphere\t0\tSuccess\t-0.7839301861678557\t0.6208489858369999\t3.1863789955259723\t4\t11\t11\t0",
  "rosen\tCGsphere\t3\tSuccess\t-0.7839301861617111\t0.6208489858447587\t3.186378995525973\t3\t9\t7\t0",
  "rosen\tCGsphere\t0\tSuccess\t-0.7839301861617111\t0.6208489858447587\t3.186378995525973\t3\t9\t7\t0",
  "rosen\tAGDsphere\t3\tMaxIters\t-0.7839301193522076\t0.6208490702034059\t3.186378995529536\t3\t12\t12\t0",
  "rosen\tAGDsphere\t0\tSuccess\t-0.783930186166909\t0.6208489858381951\t3.1863789955259727\t7\t20\t20\t0",
  "rosen\tMGDsphere\t3\tMaxIters\t-0.7839316596785372\t0.6208471252680923\t3.1863789972592724\t3\t12\t12\t0",
  "rosen\tMGDsphere\t0\tSuccess\t-0.7839301861624297\t0.6208489858438512\t3.1863789955259727\t9\t24\t24\t0",
  "rosen\tBFGSsphere\t3\tMaxIters\t-0.7839301860432096\t0.620848985994387\t3.186378995525972\t3\t10\t10\t0",
  "rosen\tBFGSsphere\t0\tSuccess\t-0.7839301861672502\t0.6208489858377644\t3.186378995525973\t4\t11\t11\t0",
  "rosen\tLBFGSsphere\t3\tSuccess\t-0.7839301861659207\t0.6208489858394431\t3.1863789955259727\t3\t10\t10\t0",
  "rosen\tLBFGSsphere\t0\tSuccess\t-0.7839301861659207\t0.6208489858394431\t3.1863789955259727\t3\t10\t10\t0",
  "rosen\tNGMRESsphere\t3\tSuccess\t-0.7839301861668613\t0.6208489858382554\t3.1863789955259727\t3\t7\t7\t0",
  "rosen\tNGMRESsphere\t0\tSuccess\t-0.7839301861668613\t0.6208489858382554\t3.1863789955259727\t3\t7\t7\t0",
  "rosen\tOACCELsphere\t3\tSuccess\t-0.7839301861663348\t0.6208489858389202\t3.1863789955259727\t3\t7\t7\t0",
  "rosen\tOACCELsphere\t0\tSuccess\t-0.7839301861663348\t0.6208489858389202\t3.1863789955259727\t3\t7\t7\t0",
  "rosen\tAdamsphere\t3\tMaxIters\t-0.8774225390127907\t0.5384102636562206\t8.882090786426254\t3\t4\t4\t0",
  "rosen\tAdamsphere\t0\tMaxIters\t-13.77061077521612\t14.659617009631313\t3.0616719042694247e6\t1000\t1001\t1001\t0",
  "stiefel\tGDstiefel\t1\tMaxIters\t0.9944653987679074 -0.1035557329457828 0.017742063787739907 0.1000593599894983 0.8819873499214662 -0.46052843458032117\t5.8015961939028315\t1\t2\t2\t0",
  "stiefel\tGDstiefel\t3\tMaxIters\t0.8795706101585434 -0.1748455952846103 -0.44247549034369993 -0.2356927286043036 0.6477487384728259 -0.7244794748575845\t4.855074952045019\t3\t4\t4\t0",
  "stiefel\tGDstiefel\t10\tMaxIters\t0.5752046029664382 0.30203485031226035 -0.760206954666341 -0.7861166005926645 0.46107275417034627 -0.411621920740871\t4.446712543465365\t10\t13\t13\t0",
  "stiefel\tGDstiefel\t0\tSuccess\t-0.9691223456333198 0.1787520674900226 -0.16985163396981118 -0.2145356382963838 -0.9508158505261824 0.223436072040607\t2.8030825825534516\t112\t119\t119\t0",
  "stiefel\tLBFGSstiefel\t1\tMaxIters\t0.9944653987679078 -0.10355573294578302 0.017742063787739897 0.1000593599894985 0.8819873499214649 -0.46052843458032394\t5.801596193902828\t1\t2\t2\t0",
  "stiefel\tLBFGSstiefel\t3\tMaxIters\t0.8353264744665272 -0.3528830430203234 -0.42154861997643744 -0.21838867764512732 0.49069601340963026 -0.8435187063132085\t4.870540458331618\t3\t4\t4\t0",
  "stiefel\tLBFGSstiefel\t10\tMaxIters\t-0.9485312650526819 0.2488352943422567 -0.19588117701086383 -0.28941820986362915 -0.9322278440835644 0.21722878841130147\t2.807725387015526\t10\t14\t14\t0",
  "stiefel\tLBFGSstiefel\t0\tSuccess\t-0.9691223460526057 0.1787520666967267 -0.16985163241235574 -0.2145356369499178 -0.9508158511369067 0.2234360707345458\t2.803082582553451\t17\t21\t21\t0",
  "stiefel\tCGstiefelchol\t1\tMaxIters\t0.79405677299211 0.46918204785114925 -0.38644798775457284 -0.6000806300067658 0.5038043686543564 -0.6213568987429594\t4.986614099874749\t1\t3\t2\t0",
  "stiefel\tCGstiefelchol\t3\tMaxIters\t-0.006489411652118044 0.5445832514766309 -0.8386816855920665 -0.9942155775655279 0.08640233754759226 0.06379671929147743\t3.9426654467081437\t3\t8\t5\t0",
  "stiefel\tCGstiefelchol\t10\tMaxIters\t-0.9690302788782447 0.17905121518816317 -0.17006169750062283 -0.21489759512312315 -0.9507079112494995 0.22354751418414845\t2.803082742708522\t10\t22\t12\t0",
  "stiefel\tCGstiefelchol\t0\tSuccess\t-0.9691223462632422 0.1787520661960588 -0.1698516317374323 -0.21453563619328106 -0.9508158514762868 0.2234360700168354\t2.8030825825534516\t24\t51\t27\t0"
)

# Optim's Adam, AdaMax and KrylovTrustRegion and the first order methods on
# manifolds (the last section of tools/optimjl-probe.jl): problem, method,
# maxiters, retcode, u, objective, iterations, fevals, gevals, hevals; the
# Stiefel rows give u as one space separated field
ojl_SA <- matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), 3, 3)
ojl_SB <- matrix(c(1, 0, 1, 0, 1, 1), 3, 2)
ojl_stief <- function(X, p) {
  s <- 0
  for (j in 1:2) for (i in 1:3) {
    ax <- 0
    for (k in 1:3) ax <- ax + ojl_SA[i, k] * X[k, j]
    s <- s + X[i, j] * ax + ojl_SB[i, j] * X[i, j]
  }
  s
}
ojl_stiefg <- function(X, p) {
  G <- matrix(0, 3, 2)
  for (j in 1:2) for (i in 1:3) {
    ax <- 0
    for (k in 1:3) ax <- ax + ojl_SA[i, k] * X[k, j]
    G[i, j] <- 2 * ax + ojl_SB[i, j]
  }
  G
}
ojl_ext_methods <- local({
  S <- Optim$Sphere()
  list(
    Adam = function() Optim$Adam(), Adam2 = function() Optim$Adam(alpha = 0.01),
    Adamsched = function() Optim$Adam(alpha = function(t) 0.01 / sqrt(t)), AdaMax = function() Optim$AdaMax(),
    KrylovTrustRegion = function() Optim$KrylovTrustRegion(),
    GDsphere = function() GradientDescent(manifold = S), CGsphere = function() ConjugateGradient(manifold = S),
    AGDsphere = function() AcceleratedGradientDescent(manifold = S),
    MGDsphere = function() MomentumGradientDescent(manifold = S), BFGSsphere = function() Optim$BFGS(manifold = S),
    LBFGSsphere = function() Optim$LBFGS(manifold = S), NGMRESsphere = function() NGMRES(manifold = S),
    OACCELsphere = function() OACCEL(manifold = S), Adamsphere = function() Optim$Adam(alpha = 0.01, manifold = S),
    GDstiefel = function() GradientDescent(manifold = Optim$Stiefel()),
    LBFGSstiefel = function() Optim$LBFGS(manifold = Optim$Stiefel()),
    CGstiefelchol = function() ConjugateGradient(manifold = Optim$Stiefel("CholQR")))
})

test_that("Optim's Adam, AdaMax, KrylovTrustRegion and the manifolds follow Optim.jl", {
  probs <- c(ojl_probs[c("quad", "rosen")],
             list(stiefel = OptimizationProblem(ojl_stief, matrix(c(1, 0, 0.5, 0, 1, 0.5), 3, 2), 1,
                                                grad = ojl_stiefg)))
  for (row in strsplit(ojl_ext_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[1:3], collapse = " ")
    s <- suppressWarnings(if (n == 0L) solve(probs[[row[1]]], ojl_ext_methods[[row[2]]]())
                          else solve(probs[[row[1]]], ojl_ext_methods[[row[2]]](), maxiters = n))
    st <- s$stats
    cnt <- as.integer(c(st$iterations, st$fevals, st$gevals, st$hevals))
    if (row[1] == "stiefel") {
      ju <- as.numeric(strsplit(row[5], " ")[[1]]); jo <- as.numeric(row[6]); jc <- as.integer(row[7:10])
    } else {
      ju <- as.numeric(row[5:6]); jo <- as.numeric(row[7]); jc <- as.integer(row[8:11])
    }
    rel <- max(abs(as.numeric(s$u) - ju) / pmax(1, abs(ju)), abs(s$objective - jo) / max(1, abs(jo)))
    if (row[1] == "stiefel") {
      # the Stiefel retractions factor with LAPACK (svd, chol), whose rounding
      # differs from Julia's in the last bits: the first iterates agree to
      # 1e-9, and the converged runs stop a few iterations apart
      if (n > 0L) {
        expect_lt(rel, 1e-9, label = info)
        expect_identical(cnt, jc, info = info)
      } else {
        expect_lt(rel, 1e-6, label = info)
      }
    } else {
      expect_lt(rel, 1e-12, label = info)
      expect_identical(cnt, jc, info = info)
    }
    expect_identical(s$retcode, row[4], info = info)
  }
  # the iterates stay on the manifold
  s <- solve(ojl_probs$rosen, Optim$LBFGS(manifold = Optim$Sphere()))
  expect_equal(sum(s$u * s$u), 1, tolerance = 1e-14)
  s <- solve(probs$stiefel, GradientDescent(manifold = Optim$Stiefel()), maxiters = 5)
  X <- matrix(s$u, 3, 2)
  expect_equal(crossprod(X), diag(2), tolerance = 1e-12)
})

# PowerManifold (two unit circles) and ProductManifold (a unit circle and the
# plane) on f(x) = x'Px + c'x (tools/optimjl-probe.jl, the manifolds rows)
ojl_pm_ref <- c(
  "manifolds\tGDpower\t1\tMaxIters\t-0.6461615499719696 0.7632006625637989 -0.35631009252831203 0.9343677637646034\t0.3213893230356144\t1\t2\t2\t0",
  "manifolds\tGDpower\t2\tMaxIters\t-0.544925386885822 0.8384845393490193 -0.5135949181900541 0.8580327849268649\t0.2032493243155571\t2\t5\t5\t0",
  "manifolds\tGDpower\t3\tMaxIters\t-0.5579750156080174 0.8298577480250654 -0.50423156580341 0.8635684848623422\t0.20238409771900123\t3\t6\t6\t0",
  "manifolds\tGDpower\t5\tMaxIters\t-0.5560215308394474 0.8311678875190964 -0.5039585285547761 0.8637278515232125\t0.20235180093066063\t5\t8\t8\t0",
  "manifolds\tGDpower\t10\tMaxIters\t-0.5555874324955307 0.8314581197228302 -0.5037288432320106 0.8638618248864459\t0.20235047181419552\t10\t13\t13\t0",
  "manifolds\tGDpower\t25\tSuccess\t-0.5555954800311215 0.8314527422343304 -0.5037331401231601 0.8638593193001166\t0.20235047133806794\t22\t25\t25\t0",
  "manifolds\tGDpower\t0\tSuccess\t-0.5555954800311215 0.8314527422343304 -0.5037331401231601 0.8638593193001166\t0.20235047133806794\t22\t25\t25\t0",
  "manifolds\tLBFGSpower\t1\tMaxIters\t-0.6461615499719697 0.7632006625637989 -0.35631009252831275 0.9343677637646031\t0.3213893230356135\t1\t2\t2\t0",
  "manifolds\tLBFGSpower\t2\tMaxIters\t-0.5185482665246521 0.8550483584478005 -0.4964765081825395 0.8680501580109716\t0.2102130348407314\t2\t4\t4\t0",
  "manifolds\tLBFGSpower\t3\tMaxIters\t-0.5631166164323214 0.8263774417890494 -0.5109413339887674 0.8596155845619474\t0.2029504571583034\t3\t5\t5\t0",
  "manifolds\tLBFGSpower\t5\tMaxIters\t-0.5555923615421844 0.8314548260705321 -0.5037369227474393 0.8638571135673657\t0.20235047143692997\t5\t7\t7\t0",
  "manifolds\tLBFGSpower\t10\tSuccess\t-0.5555954806153437 0.8314527418439399 -0.5037331404323292 0.863859319119834\t0.20235047133806794\t7\t9\t9\t0",
  "manifolds\tLBFGSpower\t25\tSuccess\t-0.5555954806153437 0.8314527418439399 -0.5037331404323292 0.863859319119834\t0.20235047133806794\t7\t9\t9\t0",
  "manifolds\tLBFGSpower\t0\tSuccess\t-0.5555954806153437 0.8314527418439399 -0.5037331404323292 0.863859319119834\t0.20235047133806794\t7\t9\t9\t0",
  "manifolds\tCGpower\t1\tMaxIters\t-0.5045333142728983 0.8633922253465136 -0.6115215076172106 0.7912278089915533\t0.2600830914455893\t1\t3\t2\t0",
  "manifolds\tCGpower\t2\tMaxIters\t-0.5629402513424894 0.8264975943210331 -0.5036557906066966 0.863904418664671\t0.20264716269991068\t2\t5\t3\t0",
  "manifolds\tCGpower\t3\tMaxIters\t-0.5558497075145229 0.8312828054615466 -0.5025971353214932 0.8645207455964423\t0.2023559018125367\t3\t7\t4\t0",
  "manifolds\tCGpower\t5\tMaxIters\t-0.555595547123695 0.8314526974015551 -0.5037331249132181 0.863859328169329\t0.2023504713380928\t5\t11\t6\t0",
  "manifolds\tCGpower\t10\tSuccess\t-0.555595479926728 0.8314527423040885 -0.5037331421122827 0.8638593181402205\t0.20235047133806816\t10\t20\t11\t0",
  "manifolds\tCGpower\t25\tSuccess\t-0.555595479926728 0.8314527423040885 -0.5037331421122827 0.8638593181402205\t0.20235047133806816\t10\t20\t11\t0",
  "manifolds\tCGpower\t0\tSuccess\t-0.555595479926728 0.8314527423040885 -0.5037331421122827 0.8638593181402205\t0.20235047133806816\t10\t20\t11\t0",
  "manifolds\tGDproduct\t1\tMaxIters\t-0.40906152401801577 0.9125068052173958 0.5432181854798125 -0.48035968618248637\t0.3191790489238835\t1\t3\t3\t0",
  "manifolds\tGDproduct\t2\tMaxIters\t-0.4569173422270205 0.8895091581159781 -0.008247513012149568 -0.6202115168551261\t-1.1032535772674108\t2\t5\t5\t0",
  "manifolds\tGDproduct\t3\tMaxIters\t-0.4334835032524266 0.9011615018452594 -0.017800364516072427 -0.5956316147777925\t-1.1121186254556916\t3\t6\t6\t0",
  "manifolds\tGDproduct\t5\tMaxIters\t-0.4291831481486077 0.9032174850750233 -0.030666162642395807 -0.5617370419262616\t-1.1174973801055512\t5\t8\t8\t0",
  "manifolds\tGDproduct\t10\tMaxIters\t-0.43513442765353505 0.9003654979343836 -0.044951420359871416 -0.5130177278456882\t-1.1219353782533124\t10\t13\t13\t0",
  "manifolds\tGDproduct\t25\tMaxIters\t-0.44000176347179554 0.8979969087595515 -0.05534365000815217 -0.47695964341689573\t-1.1230455535935266\t25\t28\t28\t0",
  "manifolds\tGDproduct\t0\tSuccess\t-0.44047988158045487 0.8977624819086998 -0.05637372415817482 -0.47338634224730125\t-1.1230546544632638\t104\t107\t107\t0",
  "manifolds\tLBFGSproduct\t1\tMaxIters\t-0.40906152401801565 0.9125068052173958 0.5432181854798125 -0.48035968618248637\t0.3191790489238835\t1\t3\t3\t0",
  "manifolds\tLBFGSproduct\t2\tMaxIters\t-0.39967308829543224 0.9166577455585 -0.03936370320515714 -0.5637557143438647\t-1.113655810003909\t2\t4\t4\t0",
  "manifolds\tLBFGSproduct\t3\tMaxIters\t-0.4330776479159719 0.9013566169256035 -0.03237031953676873 -0.5626643476742585\t-1.117421202611642\t3\t5\t5\t0",
  "manifolds\tLBFGSproduct\t5\tMaxIters\t-0.4404447780515742 0.8977797043189935 -0.056598258846431104 -0.47332421562278093\t-1.1230544735848138\t5\t7\t7\t0",
  "manifolds\tLBFGSproduct\t10\tSuccess\t-0.44047988307176705 0.8977624811769995 -0.05637372737464741 -0.47338633109105993\t-1.123054654463264\t9\t11\t11\t0",
  "manifolds\tLBFGSproduct\t25\tSuccess\t-0.44047988307176705 0.8977624811769995 -0.05637372737464741 -0.47338633109105993\t-1.123054654463264\t9\t11\t11\t0",
  "manifolds\tLBFGSproduct\t0\tSuccess\t-0.44047988307176705 0.8977624811769995 -0.05637372737464741 -0.47338633109105993\t-1.123054654463264\t9\t11\t11\t0",
  "manifolds\tCGproduct\t1\tMaxIters\t-0.25311549228111957 0.9674360689819698 0.22733602556298885 -0.2743495818889057\t-0.4648921711456644\t1\t3\t2\t0",
  "manifolds\tCGproduct\t2\tMaxIters\t-0.5165311618982991 0.8562683917954658 0.031822969452863316 -0.48029432953108897\t-1.0718714725430005\t2\t5\t3\t0",
  "manifolds\tCGproduct\t3\tMaxIters\t-0.43673197097843613 0.8995916771098377 -0.052419168792561914 -0.4819938950031592\t-1.1229549037406366\t3\t7\t4\t0",
  "manifolds\tCGproduct\t5\tMaxIters\t-0.4410037748306609 0.8975052482214841 -0.05533500191232158 -0.4742102746737143\t-1.1230501057850475\t5\t11\t6\t0",
  "manifolds\tCGproduct\t10\tMaxIters\t-0.4404797165218283 0.8977625628933577 -0.05637383695604433 -0.47338655749013697\t-1.1230546544630504\t10\t21\t11\t0",
  "manifolds\tCGproduct\t25\tSuccess\t-0.4404798833576851 0.8977624810367162 -0.056373727175418865 -0.4733863359979632\t-1.1230546544632638\t19\t40\t22\t0",
  "manifolds\tCGproduct\t0\tSuccess\t-0.4404798833576851 0.8977624810367162 -0.056373727175418865 -0.4733863359979632\t-1.1230546544632638\t19\t40\t22\t0"
)
ojl_PA <- matrix(c(3, 1, 0, 0.5, 1, 2, 0.5, 0, 0, 0.5, 4, 1, 0.5, 0, 1, 1), 4, byrow = TRUE)
ojl_PC <- c(1, -2, 0.5, 1.5)
ojl_pmf <- function(x, p) {
  s <- 0
  for (i in 1:4) {
    ax <- 0
    for (k in 1:4) ax <- ax + ojl_PA[i, k] * x[k]
    s <- s + x[i] * ax + ojl_PC[i] * x[i]
  }
  s
}
ojl_pmg <- function(x, p) {
  G <- numeric(4)
  for (i in 1:4) {
    ax <- 0
    for (k in 1:4) ax <- ax + ojl_PA[i, k] * x[k]
    G[i] <- 2 * ax + ojl_PC[i]
  }
  G
}

test_that("PowerManifold and ProductManifold follow Optim.jl", {
  prob <- OptimizationProblem(ojl_pmf, c(0.6, 0.8, -0.5, 0.2), 1, grad = ojl_pmg)
  PW <- Optim$PowerManifold(Optim$Sphere(), 2, 2)
  PR <- Optim$ProductManifold(Optim$Sphere(), Optim$Flat(), 2, 2)
  ms <- list(GDpower = function() GradientDescent(manifold = PW), LBFGSpower = function() Optim$LBFGS(manifold = PW),
             CGpower = function() ConjugateGradient(manifold = PW),
             GDproduct = function() GradientDescent(manifold = PR), LBFGSproduct = function() Optim$LBFGS(manifold = PR),
             CGproduct = function() ConjugateGradient(manifold = PR))
  for (row in strsplit(ojl_pm_ref, "\t")) {
    n <- as.integer(row[3])
    info <- paste(row[2:3], collapse = " ")
    s <- if (n == 0L) solve(prob, ms[[row[2]]]()) else solve(prob, ms[[row[2]]](), maxiters = n)
    ju <- as.numeric(strsplit(row[5], " ")[[1]])
    jo <- as.numeric(row[6])
    rel <- max(abs(s$u - ju) / pmax(1, abs(ju)), abs(s$objective - jo) / max(1, abs(jo)))
    expect_lt(rel, 1e-12, label = info)
    st <- s$stats
    expect_identical(as.integer(c(st$iterations, st$fevals, st$gevals, st$hevals)), as.integer(row[7:10]),
                     info = info)
    expect_identical(s$retcode, row[4], info = info)
  }
  s <- solve(prob, Optim$LBFGS(manifold = PW))
  expect_equal(c(sum(s$u[1:2] * s$u[1:2]), sum(s$u[3:4] * s$u[3:4])), c(1, 1), tolerance = 1e-14)
  expect_error(solve(OptimizationProblem(ojl_pmf, c(0.6, 0.8, -0.5), 1, grad = ojl_pmg),
                    GradientDescent(manifold = PR)), "AssertionError")
  expect_error(Optim$PowerManifold("Sphere", 2, 2), "inner manifolds")
})
