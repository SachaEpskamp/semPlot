# Regression tests for semPlot (plain stopifnot-style, no test framework).
# Structure: baseline pipeline tests, then per-fix regression blocks.
# Long-running blocks (psychonetrics model fitting) are skipped on CRAN
# (set NOT_CRAN=true to run them).
library(semPlot)
library(lavaan)
options(warn = 1)

ok <- 0L; fail <- 0L; fails <- character(0)
check <- function(label, expr){
  res <- tryCatch(isTRUE(expr), error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) { ok <<- ok + 1L; cat("PASS:", label, "\n") }
  else { fail <<- fail + 1L; fails <<- c(fails, label)
         cat("FAIL:", label, if (!is.null(attr(res,"msg"))) paste0(" [", attr(res,"msg"), "]") else "", "\n") }
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
not_cran <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

HS <- HolzingerSwineford1939

## ================= BASELINE =================
fit_cfa <- quiet(cfa(" visual =~ x1 + x2 + x3\n textual =~ x4 + x5 + x6 ", data = HS))
check("B1 lavaan CFA imports", { m <- quiet(semPlotModel(fit_cfa)); nrow(m@Pars) > 10 })
check("B2 lavaan CFA renders", { p <- quiet(semPaths(fit_cfa, DoNotPlot = TRUE)); inherits(p, "qgraph") })

fit_mg <- quiet(cfa(" visual =~ x1 + x2 + x3 ", data = HS, group = "school"))
check("B3 multigroup renders (list of qgraph)", {
  p <- quiet(semPaths(fit_mg, DoNotPlot = TRUE)); is.list(p) && all(sapply(p, inherits, "qgraph")) })

check("B4 semSyntax roundtrip", {
  syn <- quiet(capture.output(s <- semSyntax(fit_cfa, "lavaan")))
  fit2 <- quiet(cfa(s, data = HS))
  inherits(fit2, "lavaan") && length(coef(fit2)) == length(coef(fit_cfa)) })

d_lm <- local({ set.seed(1); d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- 0.5*d$x - 0.7*d$z + rnorm(200); d })
fit_lm <- lm(y ~ x + z, data = d_lm)
check("B5 lm imports", { m <- quiet(semPlotModel(fit_lm)); sum(m@Pars$edge == "->") == 2 })

fit_sam <- quiet(sam(" visual =~ x1+x2+x3\n textual =~ x4+x5+x6\n visual ~ textual ", data = HS))
check("B6 sam() renders", { p <- quiet(semPaths(fit_sam, DoNotPlot = TRUE)); inherits(p, "qgraph") || is.list(p) })

check("B7 lavaan model string imports", {
  m <- quiet(semPlotModel(" visual =~ x1 + x2 + x3 "))
  inherits(m, "semPlotModel") })



## ================= item01 =================

# Item 1: AST-based `+` handling (formerly deparse-split crash)
d1 <- local({ set.seed(2); d <- data.frame(x = rnorm(100), z = rnorm(100)); d$y <- d$x - d$z + rnorm(100); d })

check("T1a inline lm with + in formula renders", {
  p <- quiet(semPaths(lm(y ~ x + z, data = d1), DoNotPlot = TRUE)); inherits(p, "qgraph") })

check("T1b inline cfa call with + in model string renders", {
  p <- quiet(semPaths(cfa("visual =~ x1 + x2 + x3", data = HS), DoNotPlot = TRUE)); inherits(p, "qgraph") })

check("T1c index arithmetic in argument works", {
  mods <- list(fit_cfa, fit_cfa)
  m <- quiet(semPlotModel(mods[[1 + 1]])); inherits(m, "semPlotModel") })

check("T1d combining two lavaan fits with + still works", {
  fitA <- quiet(cfa(" visual =~ x1 + x2 + x3 ", data = HS))
  fitB <- quiet(cfa(" textual =~ x4 + x5 + x6 ", data = HS))
  m <- quiet(semPlotModel(fitA + fitB))
  inherits(m, "semPlotModel") && length(m@Original) == 2 &&
    all(c("visual","textual") %in% m@Vars$name) })

check("T1e combining via semPaths renders", {
  fitA <- quiet(cfa(" visual =~ x1 + x2 + x3 ", data = HS))
  fitB <- quiet(cfa(" textual =~ x4 + x5 + x6 ", data = HS))
  p <- quiet(semPaths(fitA + fitB, DoNotPlot = TRUE))
  inherits(p, "qgraph") || (is.list(p) && all(sapply(p, inherits, "qgraph"))) })

check("T1f combining two semPlotModel objects still works", {
  mA <- quiet(semPlotModel(cfa(" visual =~ x1 + x2 + x3 ", data = HS)))
  mB <- quiet(semPlotModel(cfa(" textual =~ x4 + x5 + x6 ", data = HS)))
  m <- quiet(semPlotModel(mA + mB)); inherits(m, "semPlotModel") && length(m@Original) == 2 })

check("T1g non-model input still errors informatively", {
  e <- tryCatch({ semPlotModel(1 + 1); NULL }, error = function(e) conditionMessage(e))
  !is.null(e) && grepl("not recognized", e) })


## ================= item02 =================

# Item 2: GLM std fallback must yield raw coefficients, not NA
d2 <- local({ set.seed(3); d <- data.frame(x = rnorm(300), z = rnorm(300));
              d$y <- rbinom(300, 1, plogis(0.5*d$x - 0.7*d$z)); d })

check("T2a binomial glm std column has non-NA slopes", {
  fitg <- glm(y ~ x + z, data = d2, family = binomial)
  m <- quiet(semPlotModel(fitg))
  slopes <- m@Pars$std[m@Pars$edge == "->"]
  sum(!is.na(slopes)) == 2 })

check("T2b binomial glm std equals raw coefficients (fallback)", {
  fitg <- glm(y ~ x + z, data = d2, family = binomial)
  m <- quiet(semPlotModel(fitg))
  co <- coef(fitg)
  sl <- m@Pars[m@Pars$edge == "->", ]
  all(abs(sl$std[match(c("x","z"), sl$lhs)] - co[c("x","z")]) < 1e-10) })

check("T2c binomial glm intercept std stays NA", {
  fitg <- glm(y ~ x + z, data = d2, family = binomial)
  m <- quiet(semPlotModel(fitg))
  is.na(m@Pars$std[m@Pars$edge == "int"]) })

if (requireNamespace("rockchalk", quietly = TRUE)) check("T2d gaussian lm std still matches rockchalk", {
  fitl <- lm(y ~ x + z, data = d2)
  m <- quiet(semPlotModel(fitl))
  sc <- coef(rockchalk::standardize(fitl))
  sl <- m@Pars[m@Pars$edge == "->", ]
  all(abs(sl$std[match(c("x","z"), sl$lhs)] - sc[c("xs","zs")]) < 1e-10) })

check("T2e binomial glm renders with what='std'", {
  fitg <- glm(y ~ x + z, data = d2, family = binomial)
  p <- quiet(semPaths(fitg, what = "std", DoNotPlot = TRUE))
  inherits(p, "qgraph") })


## ================= item03 =================

# Item 3: graphics-state hygiene.
# T3a is a strong whole-suite assertion: every earlier test used DoNotPlot=TRUE,
# so if the fixes hold, NO graphics device may exist at this point.
check("T3a no graphics device opened by any DoNotPlot call so far", is.null(dev.list()))

check("T3b par(ask) restored after an error mid-plot", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  par(ask = FALSE)
  # force an error inside the plot loop AFTER par(ask=TRUE) is set:
  try(suppressWarnings(semPaths(fit_cfa, ask = TRUE,
      edge.color = function(x) x)), silent = TRUE)
  isTRUE(par("ask") == FALSE) })

check("T3c layout restored after semCors", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  quiet(try(semCors(fit_cfa), silent = TRUE))
  # after semCors exits, a new plot must occupy the full device again:
  plot.new(); pr <- par("mfg")
  pr[3] == 1 && pr[4] == 1 })


## ================= item04 =================

# Item 4: informative errors from the importer cascade
check("T4a garbage string error names the lavaan importer failure", {
  e <- tryCatch({ semPlotModel("definitely ~~~ not // a model"); NULL },
                error = function(e) conditionMessage(e))
  !is.null(e) && grepl("lavaan", e, ignore.case = TRUE) && !identical(e, "Object not recognized as SEM model") })

check("T4b unrecognized FILE reports per-importer failures", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("this is", "just some text"), f)
  e <- tryCatch({ semPlotModel(f); NULL }, error = function(e) conditionMessage(e))
  unlink(f)
  !is.null(e) && grepl("Importer attempts failed", e) && grepl("Mplus", e) })

check("T4c valid lavaan model string still imports", {
  m <- quiet(semPlotModel(" visual =~ x1 + x2 + x3 "))
  inherits(m, "semPlotModel") })


## ================= item06 =================

# Item 6: polish bundle
check("T6a predictor named intercept_time stays a path", {
  d <- local({ set.seed(4); data.frame(intercept_time = rnorm(50), y = rnorm(50)) })
  m <- quiet(semPlotModel(lm(y ~ intercept_time, data = d)))
  sum(m@Pars$edge == "->" & m@Pars$lhs == "intercept_time") == 1 &&
    sum(m@Pars$edge == "int") == 1 })

check("T6b semCors errors informatively on lm models", {
  d <- data.frame(x = rnorm(30), y = rnorm(30))
  e <- tryCatch({ semCors(quiet(semPlotModel(lm(y ~ x, data = d)))); NULL },
                error = function(e) conditionMessage(e))
  !is.null(e) && grepl("covariance", e) })

check("T6c semSyntax warns and skips interaction terms from lm", {
  d <- local({ set.seed(5); d <- data.frame(x = rnorm(80), z = rnorm(80)); d$y <- d$x * d$z + rnorm(80); d })
  m <- quiet(semPlotModel(lm(y ~ x * z, data = d)))
  w <- NULL
  s <- withCallingHandlers(
    capture.output(res <- semSyntax(m, "lavaan")),
    warning = function(cond){ w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  any(grepl("Interaction", w)) && !grepl("par[0-9]+\\*", res) })

check("T6d main-effect syntax from lm still emitted", {
  d <- data.frame(x = rnorm(30)); d$y <- d$x + rnorm(30)
  m <- quiet(semPlotModel(lm(y ~ x, data = d)))
  s <- quiet(capture.output(res <- semSyntax(m, "lavaan")))
  grepl("y ~", res) })


## ================= item08 =================

# Item 8: lavaan modernization — composites, multilevel, key-based alignment
names(HS)[names(HS) == "x9"] -> .x9nm  # keep HS intact; use a copy
HS8 <- HS

check("T8a composite (<~) edges point INTO the composite", {
  fitc <- quiet(sem(" comp <~ 1*x1 + x2 + x3\n x9 ~ comp ", data = HS8))
  m <- quiet(semPlotModel(fitc))
  ce <- m@Pars[m@Pars$composite %in% TRUE, ]
  nrow(ce) == 3 && all(ce$edge == "~>") && all(ce$rhs == "comp") &&
    all(ce$lhs %in% c("x1","x2","x3")) && "comp" %in% m@Vars$name[!m@Vars$manifest] })

check("T8b composite model renders", {
  fitc <- quiet(sem(" comp <~ 1*x1 + x2 + x3\n x9 ~ comp ", data = HS8))
  inherits(quiet(semPaths(fitc, DoNotPlot = TRUE)), "qgraph") })

check("T8c semSyntax re-emits <~ and refits", {
  fitc <- quiet(sem(" comp <~ 1*x1 + x2 + x3\n x9 ~ comp ", data = HS8))
  s <- quiet(capture.output(res <- semSyntax(fitc, "lavaan")))
  grepl("comp <~", res) &&
    { fit2 <- quiet(sem(res, data = HS8)); inherits(fit2, "lavaan") &&
        length(coef(fit2)) == length(coef(fitc)) } })

check("T8d multilevel: Pars carry Within/Between and levels are separated", {
  fit2 <- quiet(sem(" level: 1\n fw =~ y1 + y2 + y3\n level: 2\n fb =~ y1 + y2 + y3 ",
                    data = Demo.twolevel, cluster = "cluster"))
  m <- quiet(semPlotModel(fit2))
  !is.null(m@Pars$BetweenWithin) &&
    all(sort(unique(m@Pars$BetweenWithin)) == c("Between","Within")) &&
    all(m@Pars$BetweenWithin[m@Pars$lhs == "fw" | m@Pars$rhs == "fw"] == "Within") &&
    all(m@Pars$BetweenWithin[m@Pars$lhs == "fb" | m@Pars$rhs == "fb"] == "Between") })

check("T8e multilevel renders two panels (Within and Between)", {
  fit2 <- quiet(sem(" level: 1\n fw =~ y1 + y2 + y3\n level: 2\n fb =~ y1 + y2 + y3 ",
                    data = Demo.twolevel, cluster = "cluster"))
  p <- quiet(semPaths(fit2, DoNotPlot = TRUE, ask = FALSE))
  is.list(p) && length(p) == 2 && all(sapply(p, inherits, "qgraph")) })

check("T8f multilevel block covariances named Within/Between", {
  fit2 <- quiet(sem(" level: 1\n fw =~ y1 + y2 + y3\n level: 2\n fb =~ y1 + y2 + y3 ",
                    data = Demo.twolevel, cluster = "cluster"))
  m <- quiet(semPlotModel(fit2))
  identical(names(m@ObsCovs), c("Within","Between")) &&
    identical(names(m@ImpCovs), c("Within","Between")) })

check("T8g key-based alignment: defined parameters and labels intact", {
  fitd <- quiet(sem(" visual =~ x1 + a*x2 + b*x3\n ab := a*b ", data = HS8))
  m <- quiet(semPlotModel(fitd))
  # := row removed from Pars; labeled loadings a,b kept with free par ids:
  all(!grepl(":=", m@Pars$edge)) &&
    sum(m@Pars$label %in% c("a","b")) == 2 &&
    all(m@Pars$par[m@Pars$label %in% c("a","b")] > 0) })

check("T8h ordered/threshold model still imports with thresholds", {
  d <- HS8[paste0("x",1:4)]
  for (v in names(d)) d[[v]] <- cut(d[[v]], 3, labels = FALSE, ordered_result = TRUE)
  fito <- quiet(cfa(" f =~ x1 + x2 + x3 + x4 ", data = d, ordered = names(d)))
  m <- quiet(semPlotModel(fito))
  nrow(m@Thresholds) > 0 && inherits(quiet(semPaths(fito, DoNotPlot = TRUE)), "qgraph") })

check("T8i combining lavaan and lm models now works (column union)", {
  d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
  mA <- quiet(semPlotModel(cfa(" visual =~ x1 + x2 + x3 ", data = HS8)))
  mB <- quiet(semPlotModel(lm(y ~ x, data = d)))
  m <- mA + mB
  inherits(m, "semPlotModel") && "knot" %in% names(m@Pars) && "composite" %in% names(m@Pars) })


## ================= item09 =================

# Item 9: efaList (lavaan::efa) importer
check("T9a efa() single-solution imports and renders", {
  ef <- quiet(efa(data = HS[paste0("x",1:9)], nfactors = 2))
  m <- quiet(semPlotModel(ef))
  inherits(m, "semPlotModel") && sum(!m@Vars$manifest) == 2 &&
    inherits(quiet(semPaths(ef, DoNotPlot = TRUE)), "qgraph") })

check("T9b efa() multi-solution: default is last, 'which' selects", {
  ef <- quiet(efa(data = HS[paste0("x",1:9)], nfactors = 1:2))
  m_def <- quiet(semPlotModel(ef))
  m_1 <- quiet(semPlotModel(ef, which = 1))
  m_nf2 <- quiet(semPlotModel(ef, which = "nf2"))
  sum(!m_def@Vars$manifest) == 2 && sum(!m_1@Vars$manifest) == 1 &&
    sum(!m_nf2@Vars$manifest) == 2 })

check("T9c efa() invalid 'which' errors informatively", {
  ef <- quiet(efa(data = HS[paste0("x",1:9)], nfactors = 1:2))
  e <- tryCatch({ semPlotModel(ef, which = "nf9"); NULL }, error = function(e) conditionMessage(e))
  !is.null(e) && grepl("nf", e) })


## ================= item10 =================

# Item 10: psychonetrics importer (lvm + varcov, all subtypes)
if (!(requireNamespace("psychonetrics", quietly = TRUE) && not_cran)) {
  cat("SKIP item10: psychonetrics not installed or running on CRAN\n")
} else {
  pn_L <- matrix(0, 9, 3); pn_L[1:3,1] <- pn_L[4:6,2] <- pn_L[7:9,3] <- 1
  pn_lats <- c("visual","textual","speed")
  pn_lvm_cov <- quiet(psychonetrics::runmodel(psychonetrics::lvm(
    HS[paste0("x",1:9)], lambda = pn_L, latents = pn_lats)))

  check("T10a lvm(cov) imports with expected structure", {
    m <- quiet(semPlotModel(pn_lvm_cov))
    sum(m@Pars$edge == "->") == 9 &&                             # loadings
    sum(m@Pars$edge == "<->" & m@Pars$lhs != m@Pars$rhs) == 3 && # latent covs
    sum(m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs) == 12 &&# variances
    sum(m@Pars$edge == "int") == 9 &&                            # nu (nu_eta fixed at 0, dropped)
    sum(m@Pars$edge == "->" & m@Pars$fixed) == 3 })              # marker loadings kept

  check("T10b lvm(cov) renders", {
    p <- quiet(semPaths(pn_lvm_cov, DoNotPlot = TRUE)); inherits(p, "qgraph") })

  pn_lvm_ggm <- quiet(psychonetrics::runmodel(psychonetrics::lvm(
    HS[paste0("x",1:9)], lambda = pn_L, latents = pn_lats, latent = "ggm")))

  check("T10c lvm(latent=ggm): undirected latent edges, no delta rows", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    om <- psychonetrics::getmatrix(pn_lvm_ggm, "omega_zeta")
    sum(m@Pars$edge == "--") == sum(om[upper.tri(om)] != 0) &&
    all(m@Pars$lhs[m@Pars$edge == "--"] %in% pn_lats) })

  check("T10d lvm(latent=ggm): latent self-loops equal diag(sigma_zeta)", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    sz <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma_zeta")
    sl <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs & m@Pars$lhs %in% pn_lats, ]
    nrow(sl) == 3 && all(abs(sl$est[match(pn_lats, sl$lhs)] - diag(sz)) < 1e-8) })

  check("T10e lvm(latent=ggm) delta='ignore': no latent self-loops", {
    m <- quiet(semPlotModel(pn_lvm_ggm, delta = "ignore"))
    nrow(m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs & m@Pars$lhs %in% pn_lats, ]) == 0 })

  check("T10f lvm(latent=ggm) renders with undirected edges", {
    p <- quiet(semPaths(pn_lvm_ggm, DoNotPlot = TRUE))
    inherits(p, "qgraph") && sum(!p$Edgelist$directed) >= 1 })

  pn_net <- quiet(psychonetrics::runmodel(psychonetrics::varcov(
    HS[paste0("x",1:4)], type = "ggm")))

  check("T10g varcov(ggm): pure network with variances from implied sigma", {
    m <- quiet(semPlotModel(pn_net))
    sg <- psychonetrics::getmatrix(pn_net, "sigma")
    om <- psychonetrics::getmatrix(pn_net, "omega")
    sum(m@Pars$edge == "--") == sum(om[upper.tri(om)] != 0) &&
    sum(m@Pars$edge == "~>") + sum(m@Pars$edge == "->") == 0 &&
    { sl <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs, ]
      all(abs(sl$est[match(paste0("x",1:4), sl$lhs)] - diag(sg)) < 1e-8) } })

  check("T10h varcov(ggm) renders in circle layout, all edges undirected", {
    p <- quiet(semPaths(pn_net, layout = "circle", residuals = FALSE, intercepts = FALSE, DoNotPlot = TRUE))
    inherits(p, "qgraph") && all(!p$Edgelist$directed[p$Edgelist$from != p$Edgelist$to]) })

  check("T10i ggm() wrapper equals varcov(type='ggm')", {
    w <- quiet(psychonetrics::runmodel(psychonetrics::ggm(HS[paste0("x",1:4)])))
    m1 <- quiet(semPlotModel(w)); m2 <- quiet(semPlotModel(pn_net))
    identical(dim(m1@Pars), dim(m2@Pars)) && all(abs(sort(m1@Pars$est) - sort(m2@Pars$est)) < 1e-6) })

  check("T10j multigroup lvm imports with group labels and renders 2 panels", {
    mg <- quiet(psychonetrics::runmodel(psychonetrics::lvm(
      HS, vars = paste0("x",1:9), lambda = pn_L, groups = "school", latents = pn_lats)))
    m <- quiet(semPlotModel(mg))
    p <- quiet(semPaths(mg, DoNotPlot = TRUE, ask = FALSE))
    all(c("Pasteur","Grant-White") %in% m@Pars$group) &&
      is.list(p) && length(p) == 2 && all(sapply(p, inherits, "qgraph")) })

  check("T10k lvm(latent=chol) shows implied latent covariances", {
    ch <- quiet(psychonetrics::runmodel(psychonetrics::lvm(
      HS[paste0("x",1:9)], lambda = pn_L, latents = pn_lats, latent = "chol")))
    m <- quiet(semPlotModel(ch))
    sz <- psychonetrics::getmatrix(ch, "sigma_zeta")
    lc <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs %in% pn_lats & m@Pars$rhs %in% pn_lats, ]
    nrow(lc) == 6 && !any(m@Pars$edge == "--") })

  check("T10l unsupported framework gives informative error", {
    v <- pn_net; v@model <- "var1"
    e <- tryCatch({ semPlotModel(v); NULL }, error = function(e) conditionMessage(e))
    !is.null(e) && grepl("var1", e) && grepl("supported", e) })

  # --- standardized solution (psychonetrics leaves @parameters$std empty) ---

  check("T10n lvm(latent=ggm): std is populated, not all NA", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    !all(is.na(m@Pars$std)) && !any(is.na(m@Pars$std)) })

  check("T10o lvm(latent=ggm): std loadings equal lambda * sd_eta / sd_y", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    lam <- psychonetrics::getmatrix(pn_lvm_ggm, "lambda")
    bet <- psychonetrics::getmatrix(pn_lvm_ggm, "beta")
    sz <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma_zeta")
    sg <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma")
    IB <- solve(diag(nrow(bet)) - bet)
    sdL <- sqrt(diag(IB %*% sz %*% t(IB))); sdY <- sqrt(diag(sg))
    # NB: not named 'ok' -- check() evaluates this in the global frame, where
    # 'ok' is the harness' own pass counter.
    allMatch <- TRUE
    for (i in 1:9) for (j in 1:3) {
      if (lam[i,j] == 0) next
      got <- m@Pars$std[m@Pars$edge == "->" & m@Pars$lhs == pn_lats[j] & m@Pars$rhs == paste0("x",i)]
      allMatch <- allMatch && length(got) == 1 && abs(got - lam[i,j] * sdL[j] / sdY[i]) < 1e-8
    }
    allMatch })

  check("T10p lvm(latent=ggm): omega is already standardized (std == est)", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    om <- m@Pars[m@Pars$edge == "--", ]
    nrow(om) >= 1 && all(abs(om$std - om$est) < 1e-12) })

  check("T10q lvm(latent=ggm): std residual variances are proportions in [0,1]", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    se <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma_epsilon")
    sg <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma")
    rv <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs & m@Pars$lhs %in% paste0("x",1:9), ]
    all(rv$std >= 0 & rv$std <= 1) &&
      all(abs(rv$std[match(paste0("x",1:9), rv$lhs)] - diag(se)/diag(sg)) < 1e-8) })

  check("T10r lvm(latent=ggm): latent variances standardize to 1 (no latent regressions)", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    lv <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs & m@Pars$lhs %in% pn_lats, ]
    nrow(lv) == 3 && all(abs(lv$std - 1) < 1e-8) })

  check("T10s std intercepts follow the lavaan std.all convention (est / sd_y)", {
    m <- quiet(semPlotModel(pn_lvm_ggm))
    sg <- psychonetrics::getmatrix(pn_lvm_ggm, "sigma")
    it <- m@Pars[m@Pars$edge == "int", ]
    all(abs(it$std[match(paste0("x",1:9), it$rhs)] - it$est[match(paste0("x",1:9), it$rhs)]/sqrt(diag(sg))) < 1e-8) })

  check("T10t semPaths(what='std') renders without non-finite weights", {
    w <- NULL
    p <- withCallingHandlers(quiet(semPaths(pn_lvm_ggm, "std", "est", DoNotPlot = TRUE)),
                             warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
    inherits(p, "qgraph") && is.null(w) })

  check("T10u varcov(ggm): observed variances standardize to 1", {
    m <- quiet(semPlotModel(pn_net))
    sl <- m@Pars[m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs, ]
    all(abs(sl$std - 1) < 1e-8) })

  # An empty residual network leaves every omega_epsilon fixed at zero; the
  # block must still be recognized as GGM so its delta rows are handled
  # rather than leaking through as unsupported '~/~' parameters.
  check("T10v residual='ggm' with an empty network drops no parameters", {
    rg <- quiet(psychonetrics::runmodel(psychonetrics::lvm(
      HS[paste0("x",1:9)], lambda = pn_L, latents = pn_lats, residual = "ggm")))
    w <- NULL
    m <- withCallingHandlers(semPlotModel(rg),
                             warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
    is.null(w) && !any(m@Pars$edge == "~/~") && !any(is.na(m@Pars$std)) })

  check("T10w improper solution warns and yields NA std rather than wrong values", {
    bad <- pn_lvm_ggm
    # Force a negative model-implied variance for x1 (a Heywood-like case):
    bad@modelmatrices[[1]]$sigma[1,1] <- -1
    w <- NULL
    m <- withCallingHandlers(semPlotModel(bad),
                             warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
    std_x1 <- m@Pars$std[m@Pars$edge == "->" & m@Pars$rhs == "x1"]
    !is.null(w) && grepl("improper", w) && grepl("x1", w) &&
      length(std_x1) >= 1 && all(is.na(std_x1)) })

  check("T10m ObsCovs/ImpCovs populated and symmetric", {
    m <- quiet(semPlotModel(pn_lvm_cov))
    length(m@ObsCovs) == 1 && length(m@ImpCovs) == 1 &&
      isSymmetric(unname(m@ObsCovs[[1]])) && nrow(m@ImpCovs[[1]]) == 9 })
}


## ================= item11 =================

# Item 11: undirected (GGM) support in semPaths — defExo symmetry, curveAdjacent, semSyntax guard
mk_ggm_model <- function(){
  Pars <- data.frame(label = "",
    lhs = c("x1","x1","x2","x3", paste0("x",1:4)),
    edge = c(rep("--",4), rep("<->",4)),
    rhs = c("x2","x3","x4","x4", paste0("x",1:4)),
    est = c(.3,-.2,.25,.4, rep(1,4)), std = NA, group = "", fixed = FALSE,
    par = 1:8, stringsAsFactors = FALSE)
  Vars <- data.frame(name = paste0("x",1:4), manifest = TRUE, exogenous = NA,
                     stringsAsFactors = FALSE)
  m <- new("semPlotModel", Pars = Pars, Vars = Vars, Thresholds = data.frame(),
           Computed = TRUE, ObsCovs = list(), ImpCovs = list(), Original = list())
  m
}

check("T11a defExo: pure GGM has symmetric (all-endogenous) layout roles", {
  m <- mk_ggm_model()
  e <- semPlot:::defExo(m, "tree")@Vars$exogenous
  all(e == FALSE) })

check("T11b tree layout puts all pure-GGM nodes on one level", {
  p <- quiet(semPaths(mk_ggm_model(), what = "est", layout = "tree", DoNotPlot = TRUE))
  length(unique(p$layout[,2])) == 1 })

check("T11c same-level non-adjacent -- edge is curved", {
  # 3 nodes on one level; x1--x3 passes over x2 and must curve:
  Pars <- data.frame(label = "", lhs = c("x1","x1","x2"), edge = "--",
                     rhs = c("x2","x3","x3"), est = .3, std = NA, group = "",
                     fixed = FALSE, par = 1:3, stringsAsFactors = FALSE)
  Vars <- data.frame(name = paste0("x",1:3), manifest = TRUE, exogenous = NA,
                     stringsAsFactors = FALSE)
  m <- new("semPlotModel", Pars = Pars, Vars = Vars, Thresholds = data.frame(),
           Computed = TRUE, ObsCovs = list(), ImpCovs = list(), Original = list())
  p <- quiet(semPaths(m, what = "est", layout = "tree", DoNotPlot = TRUE))
  # find the x1--x3 edge (nodes 1 and 3 in layout order) and check curvature:
  ce <- p$graphAttributes$Edges$curve
  el <- cbind(p$Edgelist$from, p$Edgelist$to)
  far <- which((el[,1] == 1 & el[,2] == 3) | (el[,1] == 3 & el[,2] == 1))
  length(far) == 1 && abs(ce[far]) > 0 })

check("T11d semSyntax warns and skips -- edges", {
  m <- mk_ggm_model()
  w <- NULL
  s <- withCallingHandlers(
    capture.output(res <- semSyntax(m, "lavaan")),
    warning = function(cond){ w <<- conditionMessage(cond); invokeRestart("muffleWarning") })
  !is.null(w) && grepl("Undirected", w) && !grepl("--", res) })

check("T11e factanal loadings model unaffected by defExo change", {
  fa <- factanal(HS[paste0("x",1:9)], factors = 3)
  m <- quiet(semPlotModel(fa$loadings))
  e <- semPlot:::defExo(m, "tree")@Vars$exogenous
  # Historical behavior (verified pre-change): the two role loops mark everything
  # exogenous, and the all-exo reset flips everything to endogenous — the MIMIC
  # branch never fires for loadings models under either the old or new edge sets.
  all(e == FALSE) &&
    inherits(quiet(semPaths(m, what = "est", DoNotPlot = TRUE)), "qgraph") })


## ================= item14 =================
# Item 14: importer registration API + model validator
check("T14a registered importer is dispatched by semPlotModel/semPaths", {
  toy <- function(object, ...){
    Pars <- data.frame(label = "", lhs = "A", edge = "~>", rhs = "B", est = 0.5,
                       std = NA, group = "", fixed = FALSE, par = 1,
                       stringsAsFactors = FALSE)
    Vars <- data.frame(name = c("A","B"), manifest = TRUE, exogenous = NA,
                       stringsAsFactors = FALSE)
    m <- new("semPlotModel", Pars = Pars, Vars = Vars, Thresholds = data.frame(),
             Computed = TRUE, ObsCovs = list(), ImpCovs = list(), Original = list())
    m
  }
  registerSemPlotImporter("myToyModel", toy)
  obj <- structure(list(), class = "myToyModel")
  m <- semPlotModel(obj)
  p <- quiet(semPaths(obj, DoNotPlot = TRUE))
  inherits(m, "semPlotModel") && nrow(m@Pars) == 1 && inherits(p, "qgraph") })

check("T14b validator catches unknown edge types", {
  Pars <- data.frame(label = "", lhs = "A", edge = "=>", rhs = "B", est = 1,
                     std = NA, group = "", fixed = FALSE, par = 1, stringsAsFactors = FALSE)
  Vars <- data.frame(name = c("A","B"), manifest = TRUE, exogenous = NA, stringsAsFactors = FALSE)
  bad <- new("semPlotModel", Pars = Pars, Vars = Vars, Thresholds = data.frame(),
             Computed = TRUE, ObsCovs = list(), ImpCovs = list(), Original = list())
  e <- tryCatch({ validateSemPlotModel(bad); NULL }, error = function(e) conditionMessage(e))
  !is.null(e) && grepl("Unknown edge", e) })

check("T14c validator catches Pars referring to unknown variables", {
  Pars <- data.frame(label = "", lhs = "A", edge = "~>", rhs = "C", est = 1,
                     std = NA, group = "", fixed = FALSE, par = 1, stringsAsFactors = FALSE)
  Vars <- data.frame(name = c("A","B"), manifest = TRUE, exogenous = NA, stringsAsFactors = FALSE)
  bad <- new("semPlotModel", Pars = Pars, Vars = Vars, Thresholds = data.frame(),
             Computed = TRUE, ObsCovs = list(), ImpCovs = list(), Original = list())
  e <- tryCatch({ validateSemPlotModel(bad); NULL }, error = function(e) conditionMessage(e))
  !is.null(e) && grepl("not in Vars", e) })

check("T14d validator passes a real imported model", {
  m <- quiet(semPlotModel(fit_cfa))
  identical(validateSemPlotModel(m), m) })

## ================= item15a =================
# Item 15a: psych fa()/omega() importers
if (!requireNamespace("psych", quietly = TRUE)) {
  cat("SKIP item15a: psych not installed\n")
} else {
  fa3 <- quiet(psych::fa(HS[paste0("x",1:9)], nfactors = 3, rotate = "oblimin", fm = "ml"))

  check("T15a1 fa() imports: loadings, Phi, uniquenesses", {
    m <- quiet(semPlotModel(fa3))
    sum(m@Pars$edge == "->") == 27 &&
      sum(m@Pars$edge == "<->" & m@Pars$lhs != m@Pars$rhs) == 3 &&
      sum(m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs) == 9 &&
      sum(!m@Vars$manifest) == 3 })

  check("T15a2 fa() renders", {
    inherits(quiet(semPaths(fa3, DoNotPlot = TRUE)), "qgraph") })

  check("T15a3 fa() varimax: no factor correlations", {
    fv <- quiet(psych::fa(HS[paste0("x",1:9)], nfactors = 2, rotate = "varimax", fm = "ml"))
    m <- quiet(semPlotModel(fv))
    sum(m@Pars$edge == "<->" & m@Pars$lhs != m@Pars$rhs) == 0 })

  check("T15a4 fa() cut argument drops small loadings", {
    m0 <- quiet(semPlotModel(fa3))
    m3 <- quiet(semPlotModel(fa3, cut = 0.3))
    sum(m3@Pars$edge == "->") < sum(m0@Pars$edge == "->") &&
      all(abs(m3@Pars$est[m3@Pars$edge == "->"]) >= 0.3) })

  check("T15a5 omega() imports bifactor with g first and renders", {
    om <- quiet(psych::omega(HS[paste0("x",1:9)], nfactors = 3, plot = FALSE))
    m <- quiet(semPlotModel(om, cut = 0.01))
    lats <- m@Vars$name[!m@Vars$manifest]
    gload <- m@Pars[m@Pars$edge == "->" & m@Pars$lhs == "g", ]
    lats[1] == "g" && nrow(gload) >= 8 &&
      sum(m@Pars$edge == "<->" & m@Pars$lhs == m@Pars$rhs) == 9 &&
      sum(m@Pars$edge == "<->" & m@Pars$lhs != m@Pars$rhs) == 0 &&
      inherits(quiet(semPaths(om, bifactor = "g", layout = "tree2", DoNotPlot = TRUE)), "qgraph") })

  check("T15a6 principal() still works via delegation", {
    pc <- quiet(psych::principal(HS[paste0("x",1:9)], nfactors = 2))
    m <- quiet(semPlotModel(pc))
    inherits(m, "semPlotModel") })
}

## ================= item15b =================
# Item 15b: piecewiseSEM psem importer
if (!requireNamespace("piecewiseSEM", quietly = TRUE)) {
  cat("SKIP item15b: piecewiseSEM not installed\n")
} else {
  data(keeley, package = "piecewiseSEM")

  check("T15b1 psem imports directed paths matching coefs()", {
    mod <- quiet(piecewiseSEM::psem(
      lm(rich ~ cover, data = keeley),
      lm(cover ~ firesev, data = keeley),
      lm(firesev ~ age, data = keeley),
      data = keeley))
    co <- quiet(piecewiseSEM::coefs(mod))
    m <- quiet(semPlotModel(mod))
    dir <- m@Pars[m@Pars$edge == "~>", ]
    nrow(dir) == 3 &&
      all(abs(sort(dir$est) - sort(as.numeric(co$Estimate))) < 1e-10) &&
      all(!is.na(dir$std)) })

  check("T15b2 psem renders", {
    mod <- quiet(piecewiseSEM::psem(
      lm(rich ~ cover, data = keeley),
      lm(cover ~ firesev, data = keeley),
      data = keeley))
    inherits(quiet(semPaths(mod, DoNotPlot = TRUE)), "qgraph") })

  check("T15b3 correlated errors become <-> edges", {
    `%~~%` <- piecewiseSEM::`%~~%`
    mod <- quiet(piecewiseSEM::psem(
      lm(rich ~ firesev, data = keeley),
      lm(cover ~ firesev, data = keeley),
      rich %~~% cover,
      data = keeley))
    m <- quiet(semPlotModel(mod))
    ce <- m@Pars[m@Pars$edge == "<->", ]
    nrow(ce) == 1 && sort(c(ce$lhs, ce$rhs))[1] == "cover" &&
      !grepl("~~", paste(ce$lhs, ce$rhs)) })

  check("T15b4 glm component supported", {
    k2 <- keeley; k2$rich01 <- as.numeric(k2$rich > stats::median(k2$rich))
    mod <- quiet(piecewiseSEM::psem(
      glm(rich01 ~ cover, family = binomial, data = k2),
      lm(cover ~ firesev, data = k2),
      data = k2))
    m <- quiet(semPlotModel(mod))
    sum(m@Pars$edge == "~>") == 2 })
}

## ================= item13 =================
# Item 13: lavaan.mi pooled multiple-imputation importer
if (!requireNamespace("lavaan.mi", quietly = TRUE)) {
  cat("SKIP item13: lavaan.mi not installed\n")
} else {
  mi_fit <- local({
    set.seed(11)
    d <- HS[paste0("x",1:6)]
    d[sample(nrow(d), 30), "x1"] <- NA
    imps <- lapply(1:3, function(i){ di <- d
      for (v in names(di)){ na <- is.na(di[[v]])
        di[[v]][na] <- mean(di[[v]], na.rm = TRUE) + rnorm(sum(na), 0, sd(di[[v]], na.rm = TRUE)) }
      di })
    quiet(lavaan.mi::cfa.mi(" visual =~ x1+x2+x3\n textual =~ x4+x5+x6 ", data = imps))
  })

  check("T13a lavaan.mi imports with pooled estimates", {
    m <- quiet(semPlotModel(mi_fit))
    pe <- quiet(lavaan.mi::parameterEstimates.mi(mi_fit))
    inherits(m, "semPlotModel") && nrow(m@Pars) == nrow(pe) &&
      all(abs(sort(m@Pars$est) - sort(pe$est)) < 1e-10) })

  check("T13b lavaan.mi std column pools standardizedSolution.mi", {
    m <- quiet(semPlotModel(mi_fit))
    ss <- quiet(lavaan.mi::standardizedSolution.mi(mi_fit))
    ld <- m@Pars[m@Pars$edge == "->", ]
    rows_ok <- sapply(seq_len(nrow(ld)), function(i){
      row <- ss[ss$lhs == ld$lhs[i] & ss$op == "=~" & ss$rhs == ld$rhs[i], ]
      nrow(row) == 1 && abs(row$est.std - ld$std[i]) < 1e-10 })
    all(rows_ok) })

  check("T13c lavaan.mi renders and semCors errors informatively", {
    p <- quiet(semPaths(mi_fit, DoNotPlot = TRUE))
    e <- tryCatch({ semCors(quiet(semPlotModel(mi_fit))); NULL },
                  error = function(e) conditionMessage(e))
    inherits(p, "qgraph") && !is.null(e) && grepl("covariance", e) })

  check("T13d lavaan.mi edge structure matches complete-data fit", {
    fitc <- quiet(cfa(" visual =~ x1+x2+x3\n textual =~ x4+x5+x6 ", data = HS))
    mc <- quiet(semPlotModel(fitc)); mm <- quiet(semPlotModel(mi_fit))
    identical(table(mc@Pars$edge), table(mm@Pars$edge)) })
}

## ================= item15c =================
# Item 15c: umx compatibility (umxRAM returns MxRAMModel; OpenMx importer handles it)
not_cran0 <- function() !exists("not_cran") || isTRUE(not_cran)
if (!(requireNamespace("umx", quietly = TRUE) && not_cran0())) {
  cat("SKIP item15c: umx not installed or CRAN mode\n")
} else {
  # umx initializes its options in .onAttach, so attach it (as its users do):
  suppressMessages(require(umx, quietly = TRUE))
  check("T15c1 umxRAM model imports and renders via OpenMx importer", {
    HSu <- HS[paste0("x",1:3)]
    m1 <- quiet(umxRAM("m1", data = HSu,
      umxPath(from = "vis", to = paste0("x",1:3)),
      umxPath(var = paste0("x",1:3)),
      umxPath(var = "vis", fixedAt = 1),
      umxPath(means = paste0("x",1:3)),
      autoRun = TRUE))
    m <- quiet(semPlotModel(m1))
    p <- quiet(semPaths(m1, DoNotPlot = TRUE))
    inherits(m, "semPlotModel") && sum(m@Pars$edge == "->") == 3 && inherits(p, "qgraph") })
}

## ================= item15d =================
# Item 15d: metaSEM stage-2 (wls) importer via meta2semPlot
not_cran2 <- function() !exists("not_cran") || isTRUE(not_cran)
if (!(requireNamespace("metaSEM", quietly = TRUE) && not_cran2())) {
  cat("SKIP item15d: metaSEM not installed or CRAN mode\n")
} else {
  ms_fit2 <- tryCatch({
    f1 <- quiet(metaSEM::tssem1(metaSEM::Becker92$data, metaSEM::Becker92$n, method = "FEM"))
    A <- metaSEM::create.mxMatrix(c(0, "0.2*Spatial2Math", "0.2*Verbal2Math",
                                    0, 0, 0,
                                    0, 0, 0),
                                  type = "Full", ncol = 3, nrow = 3, byrow = TRUE, name = "A")
    S <- metaSEM::create.mxMatrix(c("0.2*ErrMath",
                                    0, 1,
                                    0, "0.2*SpatialVerbal", 1),
                                  type = "Symm", byrow = TRUE, ncol = 3, nrow = 3, name = "S")
    quiet(metaSEM::tssem2(f1, Amatrix = A, Smatrix = S, diag.constraints = FALSE))
  }, error = function(e) NULL)

  if (is.null(ms_fit2)) {
    cat("SKIP item15d: could not fit the metaSEM example\n")
  } else {
    check("T15d1 metaSEM wls imports via meta2semPlot and renders", {
      m <- quiet(semPlotModel(ms_fit2))
      inherits(m, "semPlotModel") && nrow(m@Pars) >= 4 &&
        inherits(quiet(semPaths(ms_fit2, DoNotPlot = TRUE)), "qgraph") })

    check("T15d2 metaSEM stage-1 object errors informatively", {
      f1 <- quiet(metaSEM::tssem1(metaSEM::Becker92$data, metaSEM::Becker92$n, method = "FEM"))
      e <- tryCatch({ semPlotModel(f1); NULL }, error = function(e) conditionMessage(e))
      !is.null(e) && grepl("stage-1", e) })
  }
}

## ================= item12 =================
# Item 12: blavaan support without the class-overwrite hack
not_cran3 <- function() !exists("not_cran") || isTRUE(not_cran)
if (!(requireNamespace("blavaan", quietly = TRUE) && not_cran3())) {
  cat("SKIP item12: blavaan not installed or CRAN mode\n")
} else {
  suppressMessages(require(blavaan, quietly = TRUE))
  fit_blav <- local({
    f <- function(){
      sink(tempfile()); on.exit(sink(), add = TRUE)
      quiet(blavaan::bcfa(" visual =~ x1 + x2 + x3 ", data = HS,
                          burnin = 100, sample = 100))
    }
    tryCatch(f(), error = function(e) { cat("note: blavaan fit error:", conditionMessage(e), "\n"); NULL })
  })
  if (is.null(fit_blav)) {
    cat("SKIP item12: blavaan example fit failed\n")
  } else {
    check("T12a blavaan imports with class intact and posterior std", {
      m <- quiet(semPlotModel(fit_blav))
      is(m@Original[[1]], "blavaan") &&
        sum(!is.na(m@Pars$std)) == nrow(m@Pars) &&
        inherits(m, "semPlotModel") })
    check("T12b blavaan matches lavaan edge structure and renders", {
      fitl <- quiet(cfa(" visual =~ x1 + x2 + x3 ", data = HS))
      ml <- quiet(semPlotModel(fitl)); mb <- quiet(semPlotModel(fit_blav))
      identical(table(ml@Pars$edge), table(mb@Pars$edge)) &&
        inherits(quiet(semPaths(fit_blav, DoNotPlot = TRUE)), "qgraph") })
  }
}

## ================= item15e =================
# Item 15e: PLS-SEM importers (seminr, cSEM)
not_cran4 <- function() !exists("not_cran") || isTRUE(not_cran)
if (!(requireNamespace("seminr", quietly = TRUE) && not_cran4())) {
  cat("SKIP item15e (seminr): not installed or CRAN mode\n")
} else {
  sem_fit <- quiet({
    mm <- seminr::constructs(
      seminr::composite("Image", seminr::multi_items("IMAG", 1:5)),
      seminr::composite("Satisfaction", seminr::multi_items("CUSA", 1:3), weights = seminr::mode_B),
      seminr::reflective("Loyalty", seminr::multi_items("CUSL", 1:3)))
    sm <- seminr::relationships(
      seminr::paths(from = "Image", to = c("Satisfaction","Loyalty")),
      seminr::paths(from = "Satisfaction", to = "Loyalty"))
    seminr::estimate_pls(data = seminr::mobi, measurement_model = mm, structural_model = sm)
  })
  check("T15e1 seminr: modes map to loading vs weight arrows; renders", {
    m <- quiet(semPlotModel(sem_fit))
    modeB <- m@Pars[m@Pars$edge == "~>" & m@Pars$rhs == "Satisfaction" & m@Pars$lhs %in% paste0("CUSA",1:3), ]
    refl <- m@Pars[m@Pars$edge == "->" & m@Pars$lhs %in% c("Image","Loyalty"), ]
    struct <- m@Pars[m@Pars$edge == "~>" & m@Pars$lhs %in% c("Image","Satisfaction") & m@Pars$rhs %in% c("Satisfaction","Loyalty"), ]
    nrow(modeB) == 3 && nrow(refl) == 8 && nrow(struct) == 3 &&
      inherits(quiet(semPaths(sem_fit, DoNotPlot = TRUE)), "qgraph") })
}
if (!(requireNamespace("cSEM", quietly = TRUE) && not_cran4())) {
  cat("SKIP item15e (cSEM): not installed or CRAN mode\n")
} else {
  csem_fit <- quiet({
    model <- "
      EXPE ~ IMAG
      SAT ~ EXPE
      IMAG =~ imag1 + imag2 + imag3
      EXPE =~ expe1 + expe2 + expe3
      SAT <~ sat1 + sat2 + sat3
    "
    cSEM::csem(.data = cSEM::satisfaction[, c(paste0("imag",1:3), paste0("expe",1:3), paste0("sat",1:3))],
               .model = model)
  })
  check("T15e2 cSEM: composites get weights, common factors loadings; renders", {
    m <- quiet(semPlotModel(csem_fit))
    watSAT <- m@Pars[m@Pars$edge == "~>" & m@Pars$rhs == "SAT" & grepl("^sat", m@Pars$lhs), ]
    loads <- m@Pars[m@Pars$edge == "->" & m@Pars$lhs %in% c("IMAG","EXPE"), ]
    noSATload <- !any(m@Pars$edge == "->" & m@Pars$lhs == "SAT")
    struct <- m@Pars[m@Pars$edge == "~>" & m@Pars$lhs %in% c("IMAG","EXPE") & m@Pars$rhs %in% c("EXPE","SAT"), ]
    nrow(watSAT) == 3 && nrow(loads) == 6 && noSATload && nrow(struct) == 2 &&
      inherits(quiet(semPaths(csem_fit, DoNotPlot = TRUE)), "qgraph") })
}

cat("\n==== RESULT:", ok, "passed,", fail, "failed ====\n")
if (fail > 0) stop("semPlot regression tests failed: ", paste(fails, collapse = "; "))
