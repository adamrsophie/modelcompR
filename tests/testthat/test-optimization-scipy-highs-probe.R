# Reference results of SciPy 1.18.1 (HiGHS 1.12.0) on Linux x86-64 from
# tools/highs-probe.py: status, simplex iterations and the bits of x for the
# dual simplex method with and without presolve and the default linprog.
.hs_probe_lps <- list(
  t_lp1 = list(
    c = c(-0x1p+0, -0x1p+1, 0x1p-1),
    A_ub = matrix(c(0x1p+0, 0x1p+0, 0x1p+0, 0x1p+1, 0x1p-1, 0x0p+0, 0x0p+0, 0x1p+0, -0x1p+0), 3, 3, byrow = TRUE), b_ub = c(0x1p+2, 0x1.4p+2, 0x1.8p+0),
    A_eq = matrix(c(0x1p+0, -0x1p+0, 0x1p+1), 1, 3, byrow = TRUE), b_eq = c(0x1p-1),
    lb = c(0x0p+0, -0x1p+0, -Inf), ub = c(0x1.8p+1, Inf, 0x1p+1),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 4L, x = c(0x1.7ffffffffffffp+0, 0x1p+1, 0x1.0000000000001p-1)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 4L, x = c(0x1.7ffffffffffffp+0, 0x1p+1, 0x1.0000000000001p-1)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 4L, x = c(0x1.7ffffffffffffp+0, 0x1p+1, 0x1.0000000000001p-1)))),
  t_lp2 = list(
    c = c(0x1p+0, 0x1p+0),
    A_ub = matrix(c(-0x1p+0, -0x1p+1, -0x1.8p+1, -0x1p+0), 2, 2, byrow = TRUE), b_ub = c(-0x1p+2, -0x1.8p+2),
    A_eq = NULL, b_eq = NULL,
    lb = c(0x0p+0, 0x0p+0), ub = c(Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 2L, x = c(0x1.9999999999999p+0, 0x1.3333333333333p+0)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 2L, x = c(0x1.9999999999999p+0, 0x1.3333333333333p+0)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 2L, x = c(0x1.9999999999999p+0, 0x1.3333333333333p+0)))),
  t_lpmax = list(
    c = c(-0x1p+0, -0x1p+0),
    A_ub = matrix(c(0x1p+0, 0x1p+1, 0x1.8p+1, 0x1p+0), 2, 2, byrow = TRUE), b_ub = c(0x1p+2, 0x1.8p+2),
    A_eq = NULL, b_eq = NULL,
    lb = c(0x0p+0, 0x0p+0), ub = c(Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 2L, x = c(0x1.999999999999ap+0, 0x1.333333333333p+0)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 2L, x = c(0x1.999999999999ap+0, 0x1.333333333333p+0)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 2L, x = c(0x1.999999999999ap+0, 0x1.333333333333p+0)))),
  t_lpinf = list(
    c = c(0x1p+0, 0x1p+0),
    A_ub = matrix(c(0x1p+0, 0x1p+0), 1, 2, byrow = TRUE), b_ub = c(-0x1p+0),
    A_eq = NULL, b_eq = NULL,
    lb = c(0x0p+0, 0x0p+0), ub = c(Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Infeasible", nit = 0L, x = NULL),
      ds_default = list(presolve = TRUE, status = "Infeasible", nit = 0L, x = NULL),
      highs = list(presolve = TRUE, status = "Infeasible", nit = 0L, x = NULL))),
  t_lpunb = list(
    c = c(-0x1p+0, 0x0p+0),
    A_ub = matrix(c(0x1p+0, -0x1p+0), 1, 2, byrow = TRUE), b_ub = c(0x1p+0),
    A_eq = NULL, b_eq = NULL,
    lb = c(0x0p+0, 0x0p+0), ub = c(Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Unbounded", nit = 2L, x = NULL),
      ds_default = list(presolve = TRUE, status = "Unbounded", nit = 1L, x = NULL),
      highs = list(presolve = TRUE, status = "Unbounded", nit = 1L, x = NULL))),
  t_lp3 = list(
    c = c(0x1p+0, -0x1p+0, 0x1p+1, 0x1p-1),
    A_ub = matrix(c(0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, -0x1p+1, 0x0p+0, 0x1p+0), 2, 4, byrow = TRUE), b_ub = c(0x1.4p+2, 0x1p+1),
    A_eq = matrix(c(0x0p+0, 0x1p+0, 0x1p+0, -0x1p+0), 1, 4, byrow = TRUE), b_eq = c(0x1p+0),
    lb = c(0x1p+0, -0x1p+1, -Inf, 0x0p+0), ub = c(0x1p+2, 0x1.8p+1, Inf, 0x1p+1),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 1L, x = c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(0x1p+0, 0x1.8p+1, -0x1p+1, 0x0p+0)))),
  lp_dense = list(
    c = c(0x1.baba5af540692p-1, -0x1.ca2249a622daap-1, 0x1.43291d47563fp-4, -0x1.1f0040eed85e4p-1, 0x1.087519624621ep-1, 0x1.1b2ca5f114fdp-2, 0x1.6b3c4fcb2646cp-1, -0x1.24f2eca87ed08p-2, 0x1.251673fc5254ap-1, -0x1.2c347354035dp-1, -0x1.484c277b1c34p-2, 0x1.446945c6ba44p-3),
    A_ub = matrix(c(-0x1.4394faf443e08p-1, 0x1.59b8161f616dp-1, -0x1.28342e70f321cp-1, -0x1.4e4814af773c6p-1, 0x1.7fc5c4df55a7p-4, 0x1.ae1dd44200654p-2, -0x1.2540fb60ff28p-7, -0x1.ce938c0440f12p-1, -0x1.a97370f305f54p-2, -0x1.513a9164d6fb6p-1, 0x1.3752bc4aed7dap-1, -0x1.95f3a339064dp-2, -0x1.ef44e1fea3b92p-1, 0x1.29bf0a9ab5b34p-1, 0x1.0c35ae5c79764p-1, -0x1.0f953e7bbbd5p-4, 0x1.1ddd45779d6a6p-1, 0x1.f3f3bf9105bap-1, 0x1.7236bc8c40e98p-3, 0x1.a341068afdd6ep-1, -0x1.5bee22958a8f8p-1, -0x1.8c108c4054364p-2, 0x1.ac6037ea086b8p-1, -0x1.4884c41bca606p-1, -0x1.41cc4c8e6813ap-1, 0x1.3a19e7fc3d99p-4, 0x1.7843768aa0534p-2, 0x1.53f235d058dd4p-1, 0x1.264d01f7ab502p-1, 0x1.37f5b5f073fb8p-2, -0x1.643b55c14c468p-3, -0x1.ddea051aed27p-4, 0x1.43535f093c948p-1, -0x1.830b9d20808c2p-1, 0x1.60740137d8cdcp-1, -0x1.cdcd286972af8p-2, -0x1.a1c8f9ad17d7ap-1, 0x1.5ca8bf049484p-5, 0x1.8084656eba42p-2, -0x1.f7989ff9eaacp-2, -0x1.af5728a07c744p-2, 0x1.72e118830bb8cp-1, -0x1.3053fc44362cp-6, 0x1.f5021cbdf638p-5, -0x1.0f0523611495cp-2, 0x1.f1e9dc0625914p-1, 0x1.be502526cb8p-6, 0x1.472add0ee9f58p-2, 0x1.77304d580b7cp-6, -0x1.f500c0848dep-2, 0x1.d7e16d1cbc1eap-1, -0x1.bf2bfb1005abp-4, -0x1.d09fc7c3b1fd4p-1, -0x1.fc4dcd612388cp-1, -0x1.b68d35d547be4p-1, 0x1.6be5d86c30f6cp-1, -0x1.9c8a04fa81488p-3, 0x1.8e875d65a41bp-4, 0x1.0938b4109782p-4, 0x1.5308cce9686fp-4, -0x1.77454d947b59p-4, 0x1.cc12d81957408p-2, -0x1.c5dd779a75d04p-2, 0x1.ff751da09779p-2, -0x1.d91b7d89c145cp-1, 0x1.3509e84f58ce8p-3, 0x1.b9cc26042df6cp-2, -0x1.22d0841089edp-1, -0x1.f112b098fcd14p-1, 0x1.8335e2b6fd89p-4, -0x1.fb6a881303ebcp-2, -0x1.3c20e0f55a654p-1, -0x1.d324227e86e0ap-1, -0x1.8a57cf563f42cp-2, 0x1.3d5a570700a4p-5, 0x1.2dd751ffcc144p-1, 0x1.02f077419e3fep-1, -0x1.439065f7b1d84p-1, -0x1.3919d0910603p-4, 0x1.f67d0db5cbafep-1, -0x1.214136b29ae84p-2, -0x1.bafe0798d76ep-3, -0x1.b46a943ef4ae4p-1, -0x1.edee9c17b16acp-2, 0x1.4266833264d28p-1, -0x1.480552fab3524p-2, -0x1.15c9163f9f5fap-1, -0x1.88109b3d9efb6p-1, -0x1.b1d400cb1b518p-2, -0x1.db4067281d18ap-1, 0x1.8b299ff2e179p-2, -0x1.9da0d277a6436p-1, -0x1.01be9a0cca0e2p-1, 0x1.365af0d4fe0b4p-2, -0x1.253d1c460ba4p-3, 0x1.1147a8930f4f2p-1), 8, 12, byrow = TRUE), b_ub = c(0x1.060dcf0b02b0cp-2, 0x1.74b38abf5500ep+1, 0x1.a8d96fec9fcc6p+0, 0x1.08f89352e00eep+0, -0x1.34e471446fdep-1, 0x1.323ea377ce1ep-4, 0x1.31536547e55d4p-3, -0x1.a0f9bb055cc14p+0),
    A_eq = matrix(c(0x1.b146094ca6fp-2, 0x1.1cd897231c5ap-1, -0x1.bdbc1e8c3c038p-3, 0x1.4a12e3faa921p-2, 0x1.1b196cf8c279p-4, -0x1.d0c51af759f74p-1, -0x1.8d0ee6bb70492p-1, -0x1.44f8144fb3502p-1, 0x1.a2407adec0f5cp-1, -0x1.091d489e2f342p-1, -0x1.01431d0edba8p-4, 0x1.8d1eefd6f0af4p-1, 0x1.36e09d37244d4p-1, 0x1.0ab50f047e63cp-1, -0x1.6578bd0159014p-2, 0x1.9fd64db3c1404p-1, -0x1.53a7f8e37a268p-2, 0x1.7a069431ed334p-1, 0x1.228dfc15b162cp-2, -0x1.7a2ad852a92c4p-2, 0x1.edff77d63bd6p-4, 0x1.c414df766fe6p-2, -0x1.f8a919baa507p-4, -0x1.659b839f05ffp-3, 0x1.430ebda7eb46cp-1, -0x1.9f654911fdd2cp-1, 0x1.633eed18b0798p-3, 0x1.54216037a63e4p-2, -0x1.011138b7e69e4p-1, -0x1.ff31dd7372fe8p-1, 0x1.9e0930afb59dep-1, 0x1.953f09d523c78p-1, 0x1.888ea6672d5ep-5, -0x1.ca1fad8c5d4p-8, -0x1.d0c6e666aa57p-3, 0x1.55511ac43e33cp-2), 3, 12, byrow = TRUE), b_eq = c(0x1.a19efabee9231p-3, 0x1.9be6af31565ebp-1, -0x1.9f811f9edca46p-2),
    lb = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0), ub = c(0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1, 0x1p+1),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 10L, x = c(0x0p+0, 0x1.b20de3bae3bedp+0, 0x0p+0, 0x1.a50de51db8296p-1, 0x1.4f8dd4d589223p+0, 0x0p+0, 0x0p+0, 0x1.ab3baa1231348p+0, 0x0p+0, 0x1.7c04c5457dd5ep+0, 0x1.5a8abbaa9f56cp+0, 0x1.0df15c16f1b64p+0)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 10L, x = c(0x0p+0, 0x1.b20de3bae3bedp+0, 0x0p+0, 0x1.a50de51db8296p-1, 0x1.4f8dd4d589223p+0, 0x0p+0, 0x0p+0, 0x1.ab3baa1231348p+0, 0x0p+0, 0x1.7c04c5457dd5ep+0, 0x1.5a8abbaa9f56cp+0, 0x1.0df15c16f1b64p+0)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 10L, x = c(0x0p+0, 0x1.b20de3bae3bedp+0, 0x0p+0, 0x1.a50de51db8296p-1, 0x1.4f8dd4d589223p+0, 0x0p+0, 0x0p+0, 0x1.ab3baa1231348p+0, 0x0p+0, 0x1.7c04c5457dd5ep+0, 0x1.5a8abbaa9f56cp+0, 0x1.0df15c16f1b64p+0)))),
  lp_klee_minty = list(
    c = c(-0x1p+4, -0x1p+3, -0x1p+2, -0x1p+1, -0x1p+0),
    A_ub = matrix(c(0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+2, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+3, 0x1p+2, 0x1p+0, 0x0p+0, 0x0p+0, 0x1p+4, 0x1p+3, 0x1p+2, 0x1p+0, 0x0p+0, 0x1p+5, 0x1p+4, 0x1p+3, 0x1p+2, 0x1p+0), 5, 5, byrow = TRUE), b_ub = c(0x1.4p+2, 0x1.9p+4, 0x1.f4p+6, 0x1.388p+9, 0x1.86ap+11),
    A_eq = NULL, b_eq = NULL,
    lb = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0), ub = c(Inf, Inf, Inf, Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 1L, x = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.86ap+11)))),
  lp_free_fixed = list(
    c = c(0x1p+0, -0x1p+0, 0x1p+1, 0x1p-1),
    A_ub = matrix(c(0x1p+0, 0x1p+0, 0x0p+0, 0x1p+0, -0x1p+0, 0x1p+1, 0x1p+0, 0x0p+0), 2, 4, byrow = TRUE), b_ub = c(0x1p+2, 0x1.8p+1),
    A_eq = matrix(c(0x1p+0, -0x1p+0, 0x1p+0, 0x0p+0, 0x1p+1, -0x1p+1, 0x1p+1, 0x0p+0), 2, 4, byrow = TRUE), b_eq = c(0x1p+0, 0x1p+1),
    lb = c(-Inf, -0x1p+0, 0x0p+0, 0x0p+0), ub = c(Inf, Inf, 0x0p+0, 0x1.8p+1),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 1L, x = c(-0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(-0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 0L, x = c(-0x0p+0, -0x1p+0, 0x0p+0, 0x0p+0)))),
  lp_transport = list(
    c = c(0x1p+3, 0x1.8p+2, 0x1.4p+3, 0x1.2p+3, 0x1.2p+3, 0x1.8p+3, 0x1.ap+3, 0x1.cp+2, 0x1.cp+3, 0x1.2p+3, 0x1p+4, 0x1.4p+2),
    A_ub = matrix(c(0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0), 3, 12, byrow = TRUE), b_ub = c(0x1.4p+4, 0x1.ep+4, 0x1.9p+4),
    A_eq = matrix(c(0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1p+0), 4, 12, byrow = TRUE), b_eq = c(0x1.4p+3, 0x1.9p+4, 0x1.ep+3, 0x1.9p+4),
    lb = c(0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0), ub = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 8L, x = c(0x0p+0, 0x1.4p+4, 0x0p+0, 0x0p+0, 0x1.4p+3, 0x0p+0, 0x1.ep+3, 0x1.4p+2, 0x0p+0, 0x1.4p+2, 0x0p+0, 0x1.4p+4)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 8L, x = c(0x0p+0, 0x1.4p+4, 0x0p+0, 0x0p+0, 0x1.4p+3, 0x0p+0, 0x1.ep+3, 0x1.4p+2, 0x0p+0, 0x1.4p+2, 0x0p+0, 0x1.4p+4)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 8L, x = c(0x0p+0, 0x1.4p+4, 0x0p+0, 0x0p+0, 0x1.4p+3, 0x0p+0, 0x1.ep+3, 0x1.4p+2, 0x0p+0, 0x1.4p+2, 0x0p+0, 0x1.4p+4)))),
  lp_sparse = list(
    c = c(-0x1.91887a11f442p-4, 0x1.6766a566bc9ap-3, 0x1.930b7edbefa8p-3, 0x1.fb43f1e484f68p-3, 0x1.463e4683e7cfep-1, 0x1.652f1591b3938p-1, 0x1.b19c15e103f08p-3, 0x1.1151dbd1c242ap-1, -0x1.f5b6bf064628p-6, -0x1.63a7a594130b4p-2, -0x1.867889795042ap-1, -0x1.48f22b8e906bp-3, 0x1.adf0f3a4871eep-1, 0x1.2da64eaef0b72p-1, -0x1.ddcc1ce8bd948p-2, -0x1.469a491d4d2fp-4, 0x1.1b3ec4f9883fp-3, -0x1.a933b429cf78p-2, 0x1.33a91fd496df8p-1, 0x1.0eb8758d8163p-2),
    A_ub = matrix(c(0x1.a8d91c6d80d98p-3, 0x0p+0, -0x1.5c3ee506be52p-3, 0x1.655290723b1eep-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.2942ef17ebddcp-1, 0x1.038fc4a752888p-3, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.6a7230254e52p-4, 0x1.029f830fe4738p-3, 0x1.8e3cab359aab2p-1, 0x0p+0, 0x0p+0, 0x1.fe9056a0c6fcep-1, 0x1.469c098fa04d8p-2, -0x1.39511e28161cp-6, -0x1.d83163885a5ep-5, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.5650f1d67331cp-1, 0x1.d8fd96048cdbcp-1, 0x0p+0, 0x0p+0, 0x1.076a0faf07a14p-2, -0x1.67c097ea0109ap-1, -0x1.2854a2864141p-2, 0x0p+0, 0x0p+0, 0x1.bc47828526ad6p-1, -0x1.2ea567a692264p-2, 0x1.81de3017f4dep-1, -0x1.c8815d9109fp-8, 0x0p+0, 0x0p+0, -0x1.8dae52e021642p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.72b2382f66ed4p-1, 0x1.a641b9f41ecp-9, 0x1.4ec1cc480f2b2p-1, 0x0p+0, 0x1.2f46f861a1a1p-2, 0x1.589ef6804e13ep-1, 0x1.4397b626db2f4p-1, -0x1.7ede42e36893p-3, 0x1.18b7f8c1566c4p-2, 0x1.7fd8d60df9c78p-3, 0x1.23ce1d18522aap-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.738a396b5b402p-1, -0x1.3a947f815a9fcp-2, -0x1.2e9e70b6742bap-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.1fc2a0a9a50bp-1, 0x0p+0, -0x1.1de33c9fa44cp-4, 0x0p+0, 0x1.ce149ba9b86p-6, -0x1.4b199757a6634p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.ff911e430e8fp-3, 0x0p+0, 0x0p+0, 0x1.5dd84d289cacp-2, 0x0p+0, 0x1.3a49130ada1c8p-2, 0x1.ce74e1df4c13ap-1, 0x1.339c136a371f4p-2, 0x1.fe405e236a4bp-1, -0x1.1516e69c24f94p-2, 0x0p+0, -0x1.65efb1cf29156p-1, -0x1.c82395c11bb4ap-1, 0x0p+0, 0x1.6614be1cc0d5ep-1, -0x1.d4cf48e18c964p-1, 0x1.67b75ac753f86p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.9c8dcd5c21084p-2, 0x0p+0, 0x0p+0, 0x1.a2a89fdebb6e8p-3, -0x1.e368eb621e8cp-4, -0x1.a105e6249c716p-1, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.9bb5a06e5d188p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.5317426ed532p-5, 0x0p+0, 0x0p+0, -0x1.c74b6499eccd6p-1, 0x0p+0, 0x0p+0, 0x1.13246d94cf1cp-2, 0x1.4ab6a667de3ap-1, 0x0p+0, 0x1.2574486297728p-3, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.8e091146b4c82p-1, 0x0p+0, 0x1.be2fb5877584ep-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.b9834c8c55a98p-1, 0x1.b4dbfc8a149bcp-2, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.5dab6f4daa088p-1, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.b4f55ef85e9ap-4, -0x1.c26a3d6978d6cp-1, -0x1.8442eeac4acbp-4, 0x0p+0, 0x0p+0, 0x1.4b9300e60e964p-1, -0x1.d53f7b265631ap-1, 0x0p+0, 0x0p+0, -0x1.de15696cc238cp-1, 0x0p+0, -0x1.844fb756fafbp-4, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.a2a0ef2180d76p-1, 0x1.b4e7f92c17e5cp-1, -0x1.2a79e5a482b58p-3, -0x1.49df0bd7420dp-1, -0x1.7f50c7f077ap-6, -0x1.304ac826fe1fep-1, 0x0p+0, 0x1.e9352039dd93cp-2, -0x1.1a9c9d79d0d8p-1, 0x1.21dce004eacfp-3, 0x0p+0, 0x1.df11413ab208p-6, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.3afbe87b657ap-3, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.e6d611026fcaep-1, 0x1.a324b28917bep-4, -0x1.c24ea812cfa18p-3, -0x1.0fc0a9105cf04p-2, -0x1.58a8082d422fp-4, 0x0p+0, 0x1.4d2e92c39028ap-1, 0x0p+0, -0x1.6b19cbf129d98p-1, 0x1.507a2abfc4cbap-1, -0x1.ee2c19456f188p-1, 0x1.8a5ca0fc74a68p-1, -0x1.9c0a84c67ba4p-1, 0x0p+0, 0x0p+0, 0x1.574cee075e4bp-1, 0x1.36a2a01b3502p-3, 0x0p+0, -0x1.b475b114c3e68p-2, 0x0p+0, 0x1.8b85457ad9318p-2, 0x0p+0, 0x0p+0, -0x1.2f046dd6914ep-3, 0x0p+0, 0x1.b3adbc3e6291ep-1, 0x0p+0, -0x1.0684ee7875824p-2, -0x1.8728f9a0a8eep-4, -0x1.ed7d01c533f6ep-1, 0x0p+0, 0x0p+0, -0x1.33e90ae561028p-2, 0x0p+0, -0x1.bf8d593b247fp-2, 0x0p+0, 0x1.cd330756f2f8p-6, -0x1.7b082a5876d46p-1, 0x1.cfae734a0d954p-1, 0x0p+0, -0x1.efcfbe51cef8cp-2, 0x1.be62756595eap-2, -0x1.34b0920d8c2d8p-3, 0x1.52d2e82e226d4p-1, 0x0p+0, 0x0p+0, -0x1.f9438b78ba8bep-1, -0x1.c8b302f0c37e4p-1, -0x1.23350d0b2a0aap-1, -0x1.bcdcefeb5f52p-3, -0x1.25b40c1a78e5ap-1, 0x1.581154f855d92p-1, 0x1.4848f2d009212p-1, 0x0p+0, 0x1.9d57b5950e91cp-2, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.78e04f1e1e664p-1, 0x0p+0, 0x0p+0, -0x1.7722ab676c714p-1, 0x1.011dc998bd9cp-6, 0x1.ef635c43e03e4p-2), 12, 20, byrow = TRUE), b_ub = c(0x1.7b143ebd93418p-2, 0x1.2a1fbf6133a8ep+0, 0x1.08730fc1817b5p+2, 0x1.9dfad61a08d16p-1, -0x1.3f79bce78604ap-1, -0x1.4c4ca2a261d26p-2, 0x1.019fb9b771e39p+1, -0x1.bcc5bab1864c4p+1, 0x1.8b9ab1c84627ep+0, -0x1.04b07a4d4931cp+0, -0x1.10ded6734bb6ep+0, -0x1.00b88799f6a1fp-1),
    A_eq = matrix(c(0x1.0c240223bb08p-6, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.81babd06ca854p-1, 0x1.ed24fa1bfd62cp-1, 0x0p+0, -0x1.533673751e98p-5, 0x0p+0, 0x1.12be049aa0c6ap-1, 0x0p+0, 0x1.59fb669464416p-1, 0x0p+0, 0x1.81648e94fdc38p-2, 0x0p+0, 0x1.04f56efac3c9p-3, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.850a7b734ef48p-1, 0x0p+0, 0x1.0f9527d82240ap-1, 0x1.b973f701514cap-1, 0x1.14cfa25c4848p-7, 0x0p+0, 0x0p+0, -0x1.17a5ec3a982ep-5, 0x0p+0, 0x0p+0, -0x1.47ee82bb9cf48p-2, -0x1.82ee275150e1p-4, 0x1.39dfd1f1c5114p-1, 0x0p+0, 0x0p+0, 0x1.11f7cfb813ebp-3, 0x1.032293cffdd7p-4, 0x0p+0, -0x1.50c90edb94256p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -0x1.b58354eda0038p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.4c33992779616p-1, 0x0p+0, 0x0p+0, 0x1.9bcb5e3a0a9e8p-3, -0x1.8514c0ce6cf6p-3, 0x0p+0, 0x1.7af934971ff78p-1, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.5e87ce34e30bcp-1, 0x0p+0, 0x0p+0), 3, 20, byrow = TRUE), b_eq = c(0x1.27b724592490ep+1, 0x1.53d4ad2f77004p-5, 0x1.d55629fb06beap+0),
    lb = c(-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -Inf, -Inf, -Inf, 0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0, -Inf, -Inf, -Inf), ub = c(0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, 0x1.8p+1, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0x1p+1, 0x1p+1, 0x1p+1),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Optimal", nit = 10L, x = c(0x1.8p+1, 0x1.8p+1, -0x1p+0, 0x1.9fc93124a9389p-4, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.8p+1, 0x1.8p+1, -0x1.70b25a57bf037p+2, -0x1.57dca04ccc8e8p+1, -0x1.6cbd455c3404ap+4, 0x1.3f835ac95bcb2p+1, 0x0p+0, 0x1.a9d466dea278p+3, 0x0p+0, 0x1p+1, -0x1.382edcbb2830ep+4, -0x1.6f10519255de6p+3)),
      ds_default = list(presolve = TRUE, status = "Optimal", nit = 9L, x = c(0x1.8p+1, 0x1.8p+1, -0x1p+0, 0x1.9fc93124a93cp-4, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.8p+1, 0x1.8p+1, -0x1.70b25a57bf038p+2, -0x1.57dca04ccc8eap+1, -0x1.6cbd455c34049p+4, 0x1.3f835ac95bcb5p+1, 0x0p+0, 0x1.a9d466dea2781p+3, 0x0p+0, 0x1p+1, -0x1.382edcbb28311p+4, -0x1.6f10519255de3p+3)),
      highs = list(presolve = TRUE, status = "Optimal", nit = 9L, x = c(0x1.8p+1, 0x1.8p+1, -0x1p+0, 0x1.9fc93124a93cp-4, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.8p+1, 0x1.8p+1, -0x1.70b25a57bf038p+2, -0x1.57dca04ccc8eap+1, -0x1.6cbd455c34049p+4, 0x1.3f835ac95bcb5p+1, 0x0p+0, 0x1.a9d466dea2781p+3, 0x0p+0, 0x1p+1, -0x1.382edcbb28311p+4, -0x1.6f10519255de3p+3)))),
  lp_infeas_eq = list(
    c = c(0x1p+0, 0x1p+1, 0x1.8p+1),
    A_ub = NULL, b_ub = NULL,
    A_eq = matrix(c(0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0, 0x1p+0), 2, 3, byrow = TRUE), b_eq = c(0x1p+0, 0x1p+1),
    lb = c(0x0p+0, 0x0p+0, 0x0p+0), ub = c(Inf, Inf, Inf),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Infeasible", nit = 1L, x = NULL),
      ds_default = list(presolve = TRUE, status = "Infeasible", nit = 0L, x = NULL),
      highs = list(presolve = TRUE, status = "Infeasible", nit = 0L, x = NULL))),
  lp_unb_free = list(
    c = c(0x1p+0, -0x1p+0, 0x0p+0),
    A_ub = matrix(c(0x1p+0, 0x1p+0, 0x1p+0), 1, 3, byrow = TRUE), b_ub = c(0x1.8p+1),
    A_eq = NULL, b_eq = NULL,
    lb = c(-Inf, -Inf, 0x0p+0), ub = c(Inf, Inf, 0x1p+0),
    ref = list(
      ds_nopresolve = list(presolve = FALSE, status = "Unbounded", nit = 0L, x = NULL),
      ds_default = list(presolve = TRUE, status = "Unbounded", nit = 0L, x = NULL),
      highs = list(presolve = TRUE, status = "Unbounded", nit = 0L, x = NULL)))
)

test_that("the HiGHS LP port reproduces SciPy bit for bit", {
  for (nm in names(.hs_probe_lps)) {
    p <- .hs_probe_lps[[nm]]
    A <- rbind(p$A_ub, p$A_eq)
    rl <- c(rep(-Inf, length(p$b_ub)), p$b_eq); ru <- c(p$b_ub, p$b_eq)
    for (ly in names(p$ref)) {
      r <- p$ref[[ly]]
      s <- .hs_highs_run(p$c, A, rl, ru, p$lb, p$ub, list(presolve = r$presolve))
      expect_identical(s$status, r$status, info = paste(nm, ly))
      expect_identical(as.integer(s$nit), r$nit, info = paste(nm, ly))
      if (!is.null(r$x)) expect_identical(s$col_value, r$x, info = paste(nm, ly))
    }
  }
})
