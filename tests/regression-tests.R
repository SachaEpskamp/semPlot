# Regression tests for semPlot (plain stopifnot-style, no test framework).
# Structure: baseline pipeline tests, then per-fix regression blocks matching
# the audit items in fable_audit_semPlot.md. Long-running blocks (psychonetrics
# model fitting) are skipped on CRAN (set NOT_CRAN=true to run them).
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


cat("\n==== RESULT:", ok, "passed,", fail, "failed ====\n")
if (fail > 0) stop("semPlot regression tests failed: ", paste(fails, collapse = "; "))
