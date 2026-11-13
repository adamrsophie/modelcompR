# The SSP methods of OrdinaryDiffEqSSPRK in the out of place form of
# ssprk_perform_step.jl (perform_step! of the ConstantCaches), with the
# coefficients of ssprk_caches.jl. Under @muladd a sum of products becomes a
# chain of muladd calls (MuladdMacro: the first non product summand, or the
# first product, starts the chain; a product a * b * c is split into
# muladd(a * b, c, .)); for Julia arrays muladd(x, y, z) is x * y + z, so the
# vector arithmetic below is unfused and evaluated in that order, while the
# stage times t + c dt are fused. The stage limiter is applied to every stage
# and to the result as the in place forms do (the out of place forms of
# OrdinaryDiffEq ignore it).
#
# Each form returns the new state u, f0 = f(uprev) (integrator.fsalfirst),
# for the FSAL methods f1 = f(u, t + dt), for the adaptive ones err (the
# utilde passed to calculate_residuals) and nf, the count Julia adds.

.SSPRK_NF0 <- c(SSPRK53_H = 0L, KYKSSPRK42 = 0L)

.SSPRK_FORMS <- list(
  SSPRK22 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    f0 <- f(uprev, p, t)
    u <- lim(uprev + dt * f0, t + dt)
    k <- f(u, p, t + dt)
    u <- lim((uprev + u + dt * k) / 2, t + dt)
    list(u = u, f0 = f0, nf = 2L)
  },
  SSPRK33 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    f0 <- f(uprev, p, t)
    u <- lim(uprev + dt * f0, t + dt)
    k <- f(u, p, t + dt)
    u <- lim((u + 3 * uprev + dt * k) / 4, t + dt / 2)
    k <- f(u, p, t + dt / 2)
    u <- lim((uprev + 2 * u + (2 * dt) * k) / 3, t + dt)
    list(u = u, f0 = f0, nf = 3L)
  },
  SSPRK53 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a30 <- 0.355909775063327; a32 <- 0.644090224936674; a40 <- 0.367933791638137; a43 <- 0.632066208361863
    a52 <- 0.237593836598569; a54 <- 0.762406163401431; b10 <- 0.377268915331368; b21 <- 0.377268915331368
    b32 <- 0.242995220537396; b43 <- 0.23845893284629; b54 <- 0.287632146308408
    c1 <- 0.377268915331368; c2 <- 0.754537830662736; c3 <- 0.728985661612188; c4 <- 0.69922613593167
    f0 <- f(uprev, p, t)
    tmp <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(tmp, p, .fma(c1, dt, t))
    u <- lim(tmp + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u, p, .fma(c2, dt, t))
    tmp <- lim(a30 * uprev + a32 * u + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(tmp, p, .fma(c3, dt, t))
    tmp <- lim(a40 * uprev + a43 * tmp + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(tmp, p, .fma(c4, dt, t))
    u <- lim(a52 * u + a54 * tmp + (b54 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 5L)
  },
  SSPRK53_2N1 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a40 <- 0.571403511494104; a43 <- 0.428596488505896; b10 <- 0.443568244942995; b21 <- 0.291111420073766
    b32 <- 0.270612601278217; b43 <- 0.110577759392786; b54 <- 0.458557505351052
    c1 <- 0.443568244942995; c2 <- 0.734679665016762; c3 <- 1.005292266294979; c4 <- 0.541442494648948
    f0 <- f(uprev, p, t)
    u <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u, p, .fma(c1, dt, t))
    u <- lim(u + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u, p, .fma(c2, dt, t))
    u <- lim(u + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(u, p, .fma(c3, dt, t))
    u <- lim(a40 * uprev + a43 * u + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(u, p, .fma(c4, dt, t))
    u <- lim(u + (b54 * dt) * k, t + dt)
    # OrdinaryDiffEq counts one evaluation per step for this method
    list(u = u, f0 = f0, nf = 1L)
  },
  SSPRK53_2N2 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a30 <- 0.682342861037239; a32 <- 0.317657138962761; a50 <- 0.0452309744824; a54 <- 0.9547690255176
    b10 <- 0.465388589249323; b21 <- 0.465388589249323; b32 <- 0.124745797313998; b43 <- 0.465388589249323
    b54 <- 0.154263303748666
    c1 <- 0.465388589249323; c2 <- 0.930777178498646; c3 <- 0.42041381284771; c4 <- 0.885802402097033
    f0 <- f(uprev, p, t)
    u <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u, p, .fma(c1, dt, t))
    u <- lim(u + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u, p, .fma(c2, dt, t))
    u <- lim(a30 * uprev + a32 * u + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(u, p, .fma(c3, dt, t))
    u <- lim(u + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(u, p, .fma(c4, dt, t))
    u <- lim(a50 * uprev + a54 * u + (b54 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 5L)
  },
  SSPRK53_H = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a30 <- 0.308684154602513; a32 <- 0.691315845397487; a40 <- 0.280514990468574; a41 <- 0.270513101776498
    a43 <- 0.448971907754928; b10 <- 0.377268915331368; b21 <- 0.377268915331368; b32 <- 0.260811979144498
    b43 <- 0.169383144652957; b54 <- 0.377268915331368
    c1 <- 0.377268915331368; c2 <- 0.754537830662737; c3 <- 0.782435937433493; c4 <- 0.622731084668631
    f0 <- f(uprev, p, t)
    tmp <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(tmp, p, .fma(c1, dt, t))
    u <- lim(tmp + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u, p, .fma(c2, dt, t))
    u <- lim(a30 * uprev + a32 * u + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(u, p, .fma(c3, dt, t))
    u <- lim(a40 * uprev + a41 * tmp + a43 * u + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(u, p, .fma(c4, dt, t))
    u <- lim(u + (b54 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 5L)
  },
  SSPRK63 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a40 <- 0.476769811285196; a41 <- 0.098511733286064; a43 <- 0.42471845542874; a62 <- 0.155221702560091
    a65 <- 0.844778297439909; b10 <- 0.284220721334261; b21 <- 0.284220721334261; b32 <- 0.284220721334261
    b43 <- 0.12071378576593; b54 <- 0.284220721334261; b65 <- 0.2401034970659
    c1 <- 0.284220721334261; c2 <- 0.568441442668522; c3 <- 0.852662164002783; c4 <- 0.510854218958172
    c5 <- 0.795074940292433
    f0 <- f(uprev, p, t)
    u <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u, p, .fma(c1, dt, t))
    u2 <- lim(u + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u2, p, .fma(c2, dt, t))
    tmp <- lim(u2 + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(tmp, p, .fma(c3, dt, t))
    tmp <- lim(a40 * uprev + a41 * u + a43 * tmp + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(tmp, p, .fma(c4, dt, t))
    tmp <- lim(tmp + (b54 * dt) * k, .fma(c5, dt, t))
    k <- f(tmp, p, .fma(c5, dt, t))
    u <- lim(a62 * u2 + a65 * tmp + (b65 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 6L)
  },
  SSPRK73 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a40 <- 0.184962588071072; a43 <- 0.815037411928928; a50 <- 0.18071865657038; a51 <- 0.314831034403793
    a54 <- 0.504450309025826; a73 <- 0.120199; a76 <- 0.879801; b10 <- 0.233213863663009
    b21 <- 0.233213863663009; b32 <- 0.233213863663009; b43 <- 0.190078023865845; b54 <- 0.117644805593912
    b65 <- 0.233213863663009; b76 <- 0.205181790464579
    c1 <- 0.233213863663009; c2 <- 0.466427727326018; c3 <- 0.699641590989027; c4 <- 0.760312095463379
    c5 <- 0.574607439040817; c6 <- 0.807821302703826
    f0 <- f(uprev, p, t)
    u1 <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u1, p, .fma(c1, dt, t))
    tmp <- lim(u1 + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(tmp, p, .fma(c2, dt, t))
    u <- lim(tmp + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(u, p, .fma(c3, dt, t))
    tmp <- lim(a40 * uprev + a43 * u + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(tmp, p, .fma(c4, dt, t))
    tmp <- lim(a50 * uprev + a51 * u1 + a54 * tmp + (b54 * dt) * k, .fma(c5, dt, t))
    k <- f(tmp, p, .fma(c5, dt, t))
    tmp <- lim(tmp + (b65 * dt) * k, .fma(c6, dt, t))
    k <- f(tmp, p, .fma(c6, dt, t))
    u <- lim(a73 * u + a76 * tmp + (b76 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 7L)
  },
  SSPRK83 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a50 <- 0.421366967085359; a51 <- 0.005949401107575; a54 <- 0.572683631807067; a61 <- 0.004254010666365
    a65 <- 0.995745989333635; a72 <- 0.104380143093325; a73 <- 0.243265240906726; a76 <- 0.65235461599995
    b10 <- 0.195804015330143; b21 <- 0.195804015330143; b32 <- 0.195804015330143; b43 <- 0.195804015330143
    b54 <- 0.112133754621673; b65 <- 0.194971062960412; b76 <- 0.127733653231944; b87 <- 0.195804015330143
    c1 <- 0.195804015330143; c2 <- 0.391608030660286; c3 <- 0.587412045990429; c4 <- 0.783216061320572
    c5 <- 0.561833689734037; c6 <- 0.755247658555329; c7 <- 0.804195984669857
    f0 <- f(uprev, p, t)
    u <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u, p, .fma(c1, dt, t))
    u2 <- lim(u + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u2, p, .fma(c2, dt, t))
    u3 <- lim(u2 + (b32 * dt) * k, .fma(c3, dt, t))
    k <- f(u3, p, .fma(c3, dt, t))
    tmp <- lim(u3 + (b43 * dt) * k, .fma(c4, dt, t))
    k <- f(tmp, p, .fma(c4, dt, t))
    tmp <- lim(a50 * uprev + a51 * u + a54 * tmp + (b54 * dt) * k, .fma(c5, dt, t))
    k <- f(tmp, p, .fma(c5, dt, t))
    tmp <- lim(a61 * u + a65 * tmp + (b65 * dt) * k, .fma(c6, dt, t))
    k <- f(tmp, p, .fma(c6, dt, t))
    tmp <- lim(a72 * u2 + a73 * u3 + a76 * tmp + (b76 * dt) * k, .fma(c7, dt, t))
    k <- f(tmp, p, .fma(c7, dt, t))
    u <- lim(tmp + (b87 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 8L)
  },
  SSPRK43 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    one_third <- 1 / 3; two_thirds <- 2 * one_third
    dt_2 <- 0.5 * dt
    f0 <- f(uprev, p, t)
    u <- lim(uprev + dt_2 * f0, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_2 * k, t + dt)
    k <- f(u, p, t + dt)
    u <- u + dt_2 * k
    if (adaptive) utilde <- one_third * uprev + two_thirds * u
    u <- lim(two_thirds * uprev + one_third * u, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_2 * k, t + dt)
    list(u = u, f0 = f0, nf = 4L, err = if (adaptive) 0.5 * (utilde - u))
  },
  SSPRK432 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    dt_2 <- dt / 2
    f0 <- f(uprev, p, t)
    u <- lim(uprev + dt_2 * f0, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_2 * k, t + dt)
    k <- f(u, p, t + dt)
    u <- u + dt_2 * k
    if (adaptive) utilde <- (uprev + 2 * u) / 3
    u <- lim((2 * uprev + u) / 3, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_2 * k, t + dt)
    list(u = u, f0 = f0, nf = 4L, err = if (adaptive) utilde - u)
  },
  SSPRK932 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    dt_6 <- dt / 6; dt_3 <- dt / 3; dt_2 <- dt / 2
    f0 <- f(uprev, p, t)
    u <- lim(uprev + dt_6 * f0, t + dt_6)
    k <- f(u, p, t + dt_6)
    u <- lim(u + dt_6 * k, t + dt_3)
    k <- f(u, p, t + dt_3)
    u <- lim(u + dt_6 * k, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_6 * k, .fma(2, dt_3, t))
    k <- f(u, p, .fma(2, dt_3, t))
    u <- lim(u + dt_6 * k, .fma(5, dt_6, t))
    k <- f(u, p, .fma(5, dt_6, t))
    u <- u + dt_6 * k
    nf <- 9L
    if (adaptive) {
      k <- f(u, p, t + dt)
      nf <- 10L
      utilde <- (uprev + 6 * u + (6 * dt) * k) / 7
    }
    u <- lim((3 * uprev + dt_2 * f0 + 2 * u) / 5, t + dt_2)
    k <- f(u, p, t + dt_2)
    u <- lim(u + dt_6 * k, .fma(2, dt_3, t))
    k <- f(u, p, .fma(2, dt_3, t))
    u <- lim(u + dt_6 * k, .fma(5, dt_6, t))
    k <- f(u, p, .fma(5, dt_6, t))
    u <- lim(u + dt_6 * k, t + dt)
    list(u = u, f0 = f0, nf = nf, err = if (adaptive) utilde - u)
  },
  SSPRK54 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    b10 <- 0.39175222657189; a20 <- 0.444370493651235; a21 <- 0.555629506348765; b21 <- 0.368410593050371
    a30 <- 0.620101851488403; a32 <- 0.379898148511597; b32 <- 0.251891774271694; a40 <- 0.178079954393132
    a43 <- 0.821920045606868; b43 <- 0.544974750228521; a52 <- 0.517231671970585; a53 <- 0.096059710526147
    b53 <- 0.06369246866629; a54 <- 0.386708617503269; b54 <- 0.226007483236906
    c1 <- 0.39175222657189; c2 <- 0.58607968931154; c3 <- 0.4745423631214; c4 <- 0.935010630967653
    f0 <- f(uprev, p, t)
    u2 <- lim(uprev + (b10 * dt) * f0, .fma(c1, dt, t))
    k <- f(u2, p, .fma(c1, dt, t))
    u2 <- lim(a20 * uprev + a21 * u2 + (b21 * dt) * k, .fma(c2, dt, t))
    k <- f(u2, p, .fma(c2, dt, t))
    u3 <- lim(a30 * uprev + a32 * u2 + (b32 * dt) * k, .fma(c3, dt, t))
    k3 <- f(u3, p, .fma(c3, dt, t))
    tmp <- lim(a40 * uprev + a43 * u3 + (b43 * dt) * k3, .fma(c4, dt, t))
    k <- f(tmp, p, .fma(c4, dt, t))
    u <- lim(a52 * u2 + a53 * u3 + (b53 * dt) * k3 + a54 * tmp + (b54 * dt) * k, t + dt)
    list(u = u, f0 = f0, nf = 5L)
  },
  SSPRK104 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    dt_6 <- dt / 6; dt_3 <- dt / 3; dt_2 <- dt / 2
    f0 <- f(uprev, p, t)
    tmp <- lim(uprev + dt_6 * f0, t + dt_6)
    k <- f(tmp, p, t + dt_6)
    tmp <- lim(tmp + dt_6 * k, t + dt_3)
    k <- f(tmp, p, t + dt_3)
    tmp <- lim(tmp + dt_6 * k, t + dt_2)
    k <- f(tmp, p, t + dt_2)
    u <- lim(tmp + dt_6 * k, .fma(2, dt_3, t))
    k4 <- f(u, p, .fma(2, dt_3, t))
    tmp <- lim((3 * uprev + 2 * u + (2 * dt_6) * k4) / 5, t + dt_3)
    k <- f(tmp, p, t + dt_3)
    tmp <- lim(tmp + dt_6 * k, t + dt_2)
    k <- f(tmp, p, t + dt_2)
    tmp <- lim(tmp + dt_6 * k, .fma(2, dt_3, t))
    k <- f(tmp, p, .fma(2, dt_3, t))
    tmp <- lim(tmp + dt_6 * k, .fma(5, dt_6, t))
    k <- f(tmp, p, .fma(5, dt_6, t))
    tmp <- lim(tmp + dt_6 * k, t + dt)
    k <- f(tmp, p, t + dt)
    u <- lim((uprev + 9 * (u + dt_6 * k4) + 15 * (tmp + dt_6 * k)) / 25, t + dt)
    list(u = u, f0 = f0, nf = 10L)
  },
  # the FSAL methods start from the derivative fs of the previous step
  KYKSSPRK42 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a20 <- 0.394806441339829; a21 <- 0.605193558660171; a30 <- 0.00279730708739; a32 <- 0.99720269291261
    a40 <- 0.252860909354373; a43 <- 0.747139090645627; b10 <- 0.406584463657504; b21 <- 0.246062298456822
    b30 <- 0.013637216641451; b32 <- 0.405447122055692; b40 <- 0.016453567333598; b43 <- 0.303775146447707
    c1 <- 0.406584463657504; c2 <- 0.4921245969136438; c3 <- 0x1.d1d58a71aa2b7p-1
    u <- lim(uprev + (dt * b10) * fs, .fma(c1, dt, t))
    k <- f(u, p, .fma(c1, dt, t))
    tmp <- lim(a20 * uprev + a21 * u + (dt * b21) * k, .fma(c2, dt, t))
    k <- f(tmp, p, .fma(c2, dt, t))
    tmp <- lim(a30 * uprev + a32 * tmp + (dt * b30) * fs + (dt * b32) * k, .fma(c3, dt, t))
    k <- f(tmp, p, .fma(c3, dt, t))
    u <- lim(a40 * uprev + a43 * tmp + (dt * b40) * fs + (dt * b43) * k, t + dt)
    list(u = u, f0 = fs, f1 = f(u, p, t + dt), nf = 4L)
  },
  KYK2014DGSSPRK_3S2 = function(f, p, t, uprev, dt, lim, fs, adaptive) {
    a10 <- 1; a20 <- 0.087353119859156; a21 <- 0.912646880140844; a30 <- 0.344956917166841
    a32 <- 0.655043082833159; b10 <- 0.528005024856522; b21 <- 0.481882138633993; b30 <- 0.022826837460491
    b32 <- 0.345866039233415
    c1 <- b10; c2 <- a21 * b10 + b21
    u1 <- lim(a10 * uprev + (dt * b10) * fs, .fma(c1, dt, t))
    u2 <- lim(a20 * uprev + a21 * u1 + (dt * b21) * f(u1, p, .fma(c1, dt, t)), .fma(c2, dt, t))
    u <- lim(a30 * uprev + (dt * b30) * fs + a32 * u2 + (dt * b32) * f(u2, p, .fma(c2, dt, t)), t + dt)
    list(u = u, f0 = fs, f1 = f(u, p, t + dt), nf = 3L)
  }
)

# step and init of the method `name` for .erk_alg, with the dense output
# `dense` of .ERK_DENSE ("fsal" for the FSAL methods)
.ssprk_parts <- function(name, stage_limiter, dense) {
  form <- .SSPRK_FORMS[[name]]
  fsal <- identical(dense, "fsal")
  nf0 <- if (is.na(.SSPRK_NF0[name])) 1L else .SSPRK_NF0[[name]]
  step <- function(ctx, t, u, dt, cache) {
    p <- ctx$p
    lim <- if (is.null(stage_limiter)) function(x, tt) x else function(x, tt) stage_limiter(x, NULL, p, tt)
    fs <- cache$fsal
    nf <- 0L
    if (fsal && is.null(fs)) {
      fs <- ctx$f(u, p, t)
      nf <- 1L
    }
    r <- form(ctx$f, p, t, u, dt, lim, fs, isTRUE(ctx$adaptive))
    k <- switch(dense, ssp = list(r$f0), ends = .lsrk_k(ctx, t, dt), fsal = list(r$f0, r$f1))
    list(u = r$u, err = r$err, k = k, cache = if (fsal) list(fsal = r$f1) else list(), nf = r$nf + nf)
  }
  init <- if (fsal) function(ctx, t, u) list(fsal = ctx$f(u, ctx$p, t), nf0 = nf0) else
    function(ctx, t, u) list(nf0 = nf0)
  list(step = step, init = init)
}
