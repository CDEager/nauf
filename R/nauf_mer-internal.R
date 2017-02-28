
#' @importFrom stringr str_split
nauf_mkBlist <- function(x, frloc, drop.unused.levels = TRUE) {
  # code taken from lme4:::mkBlist and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored
  mtfe <- attr(frloc, "terms")
  mefcfe <- attr(mtfe, "mefc")
  drop.unused.levels <- TRUE
  
  frloc <- lme4::factorize(x, frloc)
  if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac), 
    list(fac = x[[3]])), frloc), error = function(e) NULL))) 
    stop("couldn't evaluate grouping factor ", deparse(x[[3]]), 
      " within model frame:", " try adding grouping factor to data ", 
      "frame explicitly if possible", call. = FALSE)
  if (all(is.na(ff))) {
    stop("Invalid grouping factor specification, ", deparse(x[[3]]),
      call. = FALSE)
  }
  if (drop.unused.levels) ff <- droplevels(ff)
  nl <- length(levels(ff))
  
  re_form <- eval(substitute(~foo, list(foo = x[[2]])))
  mtre <- attr(nauf_model_frame(re_form, frloc[!is.na(ff), , drop = FALSE]),
    "terms")
  mf <- stats::model.frame(re_form, frloc, na.action = "na.pass")
  mefcre <- attr(mtre, "mefc")
  retl <- attr(mtre, "term.labels")
  for (f in names(mefcre)[names(mefcre) %in% names(mefcfe)]) {
    if (mefcfe[[f]]$ordered) {
      mefcre[[f]] <- mefcfe[[f]]
    } else if (length(mefcfe[[f]]$levels) != length(mefcre[[f]]$levels)) {
      mefcre[[f]]$ccrena <- TRUE
      mf[, f] <- factor(mf[, f], levels = mefcre[[f]]$levels, ordered = FALSE)
      contrasts(mf[, f]) <- mefcre[[f]]$contrasts
      cn <- which(colnames(mf) == f)
      sep <- "_"
      cnm <- paste("ccrena", f, sep = sep)
      # this while loop shouldn't need to execute
      while (cnm %in% colnames(frloc)) {
        sep <- paste(sep, "_", sep = "")
        cnm <- paste("ccrena", f, sep = sep)
      }
      colnames(mf)[cn] <- cnm
      tnms <- stringr::str_split(retl, ":")
      tnms <- lapply(tnms, function(u) {
          u[u == f] <- cnm
          u <- paste(u, collapse = ":")
          u
        })
      names(mefcre)[names(mefcre) == f] <- cnm
      retl <- tnms
    } else {
      mefcre[[f]]$ccrena <- FALSE
    }
  }
  re_form <- stats::terms(stats::as.formula(paste(
    ifelse(attr(mtre, "intercept") == 1, "~1+", "~0+"),
    paste(retl, collapse = "+"))))
  attr(mf, "terms") <- NULL
  mf <- nauf_model_frame(re_form, mf)
  attr(attr(mf, "terms"), "mefc") <- mefcre
  attr(attr(mf, "terms"), "re_hasna") <- anyNA(ff)
  mm <- nauf_model_matrix(mf)
  
  sm <- Matrix::fac2sparse(ff, to = "d",
    drop.unused.levels = drop.unused.levels)
  sm <- Matrix::KhatriRao(sm, t(mm))
  dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))
  return(list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm),
    mt = attr(mf, "terms")))
}


nauf_mkReTrms <- function(bars, fr, drop.unused.levels = TRUE) {
  # code taken from lme4::mkReTrms and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored
  drop.unused.levels <- TRUE
  if (!length(bars)) {
    stop("No random effects terms specified in formula", call. = FALSE)
  }
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr, 
    "data.frame"))
  names(bars) <- lme4:::barnames(bars)
  term.names <- vapply(bars, lme4:::safeDeparse, "")
  blist <- lapply(bars, nauf_mkBlist, fr, drop.unused.levels)
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (any(diff(nl) > 0)) {
    ord <- rev(order(nl))
    blist <- blist[ord]
    nl <- nl[ord]
    term.names <- term.names[ord]
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(Matrix::rBind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1)) / 2)
  nb <- nc * nl
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)", 
      sum(nb), q))
  }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  Lambdat <- Matrix::t(do.call(Matrix::sparseMatrix, do.call(Matrix::rBind,
    lapply(seq_along(blist), function(i) {
      mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
      dd <- diag(nc[i])
      ltri <- lower.tri(dd, diag = TRUE)
      ii <- row(dd)[ltri]
      jj <- col(dd)[ltri]
      data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[, 
        jj]) + boff[i], x = as.double(rep.int(seq_along(ii), 
        rep.int(nl[i], length(ii))) + thoff[i]))
    }))))
  thet <- numeric(sum(nth))
  ll <- list(Zt = Matrix::drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x), 
    Gp = unname(c(0L, cumsum(nb))))
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(Matrix::diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else {
    asgn <- seq_along(fl)
  }
  names(fl) <- ufn
  fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll$nauf_mt <- lapply(blist, function(n) n$mt)
  return(nauf_on(ll))
}


nauf_lFormula <- function(formula, data = NULL, REML = TRUE,
                          subset, weights, na.action, offset, contrasts = NULL,
                          control = lme4::lmerControl(), ...) {
  # code taken from lme4::lFormula and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored
  control <- control$checkControl
  mf <- mc <- match.call()

  ignoreArgs <- c("start","verbose","devFunOnly","control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("lmer"), l...))
  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(nauf:::nauf_glFormula)
    if (missing(control)) mc[["control"]] <- lme4::glmerControl()
    return(eval(mc, parent.frame()))
  }

  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4:::checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  formula <- stats::as.formula(formula, env = denv)
  formula[[length(formula)]] <- lme4::expandDoubleVerts(lme4:::RHSForm(formula))
  mc$formula <- formula

  m <- match(c("data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model_frame)
  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x=.(i))))) {
      assign(i,get(i,parent.frame()),environment(fr.form))
    }
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  attr(fr,"formula") <- formula
  attr(fr,"offset") <- mf$offset
  n <- nrow(fr)

  reTrms <- nauf_mkReTrms(lme4::findbars(lme4:::RHSForm(formula)), fr)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ",
      "please use ",
      shQuote("na.action='na.omit'"),
      " or ",
      shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e6)

  fixedform <- formula
  fixedform[[length(fixedform)]] <- lme4::nobars(lme4:::RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr,"terms"), "predvars.fixed") <- attr(attr(fixedfr, "terms"),
    "predvars")

  ranform <- formula
  ranform[[length(ranform)]] <- lme4::subbars(lme4:::RHSForm(
    lme4:::reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), "predvars")

  X <- nauf_model_matrix(lme4::nobars(formula), data)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-7)
  if(is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)

  return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}

nauf_lmer <- function(formula, data = NULL, REML = TRUE,
                      control = lme4::lmerControl(), 
                      start = NULL, verbose = 0L, subset, weights, na.action,
                      offset, contrasts = NULL, devFunOnly = FALSE, ...) {
  # code taken from lme4::lmer and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored
  mc <- match.call()
  mc$na.action <- "na.pass"
  mc$contrats <- NULL
  mcout <- mc
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control)) 
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead", 
      immediate. = TRUE)
    control <- do.call(lme4::lmerControl, control)
  }
  if (!is.null(list(...)[["family"]])) {
    warning("calling lmer with 'family' is deprecated; please use glmer() instead")
    mc[[1]] <- quote(nauf:::nauf_glmer)
    if (missCtrl) 
      mc$control <- lme4::glmerControl()
    return(eval(mc, parent.frame(1L)))
  }
  mc$control <- control
  mc[[1]] <- quote(nauf:::nauf_lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  devfun <- do.call(lme4::mkLmerDevfun, c(lmod, list(start = start, 
    verbose = verbose, control = control)))
  if (devFunOnly) {
    return(devfun)
  }
  opt <- if (control$optimizer == "none") {
    list(par = NA, fval = NA, conv = 1000, message = "no optimization")
  } else {
    lme4::optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge, 
      boundary.tol = control$boundary.tol, control = control$optCtrl, 
      verbose = verbose, start = start, calc.derivs = control$calc.derivs, 
      use.last.params = control$use.last.params)
  }
  cc <- lme4:::checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
    lbound = environment(devfun)$lower)
  return(lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = nauf_off(lmod$fr), 
    mc = mcout, lme4conv = cc))
}


nauf_glFormula <- function(formula, data = NULL, family = gaussian,
                      subset, weights, na.action, offset,
                      contrasts = NULL, start, mustart, etastart,
                      control = lme4::glmerControl(), ...) {
  # code taken from lme4::glFormula and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored

  control <- control$checkControl
  mf <- mc <- match.call()

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    mc[[1]] <- quote(nauf:::nauf_lFormula)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  if (family$family %in% c("quasibinomial", "quasipoisson", "quasi")) {
    stop('"quasi" families cannot be used in glmer')
  }
  
  ignoreArgs <- c("start","verbose","devFunOnly","optimizer", "control", "nAGQ")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("glmer"), l...))

  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])

  denv <- lme4:::checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  mc$formula <- formula <- stats::as.formula(formula, env = denv)

  m <- match(c("data", "subset", "weights", "na.action", "offset",
    "mustart", "etastart"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf_model_frame)
  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x=.(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- fr.form
  fr <- eval(mf, parent.frame())
  attr(fr,"formula") <- formula
  attr(fr,"offset") <- mf$offset
  if (!missing(start) && is.list(start)) {
      attr(fr,"start") <- start$fixef
  }
  n <- nrow(fr)
  
  reTrms <- nauf_mkReTrms(lme4::findbars(lme4:::RHSForm(formula)), fr)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e6, allow.n = TRUE)

  fixedform <- formula
  fixedform[[length(fixedform)]] <- lme4::nobars(lme4:::RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr,"terms"),"predvars.fixed") <- attr(attr(fixedfr, "terms"),
    "predvars")

  ranform <- formula
  ranform[[length(ranform)]] <- lme4::subbars(lme4:::RHSForm(
    lme4:::reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), "predvars")

  X <- nauf_model_matrix(lme4::nobars(formula), data)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-7)
  if(is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)

  return(list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}


nauf_glmer <- function(formula, data = NULL, family = gaussian,
                       control = lme4::glmerControl(), start = NULL,
                       verbose = 0L, nAGQ = 1L, subset, weights, na.action, 
                       offset, contrasts = NULL, mustart, etastart,
                       devFunOnly = FALSE, ...) {
  # code taken from lme4::glmer and altered as little as possible
  # to implement nauf methods. drop.unused.levels is ignored
  if (!inherits(control, "glmerControl")) {
    if (!is.list(control)) {
      stop("'control' is not a list; use glmerControl()")
    }
    msg <- "Use control=glmerControl(..) instead of passing a list"
    if (length(cl <- class(control))) {
      msg <- paste(msg, "of class", dQuote(cl[1]))
    }
    warning(msg, immediate. = TRUE)
    control <- do.call(lme4::glmerControl, control)
  }
  mc <- mcout <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if (is.function(family)) {
    family <- family()
  }
  if (isTRUE(all.equal(family, gaussian()))) {
      warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;", 
          " please call lmer() directly")
      mc[[1]] <- quote(nauf:::nauf_lmer)
      mc["family"] <- NULL
      return(eval(mc, parent.frame()))
  }
  mc[[1]] <- quote(nauf:::nauf_glFormula)
  glmod <- eval(mc, parent.frame(1L))
  mcout$formula <- glmod$formula
  glmod$formula <- NULL
  nAGQinit <- if (control$nAGQ0initStep) {
    0L
  } else {
    1L
  }
  devfun <- do.call(lme4::mkGlmerDevfun, c(glmod, list(verbose = verbose, 
    control = control, nAGQ = nAGQinit)))
  if (nAGQ == 0 && devFunOnly) {
    return(devfun)
  }
  if (is.list(start)) {
    start.bad <- setdiff(names(start), c("theta", "fixef"))
    if (length(start.bad) > 0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s", 
        paste(start.bad, collapse = ", "), shQuote("theta"), 
        shQuote("fixef")), call. = FALSE)
    }
    if (!is.null(start$fixef) && nAGQ == 0) {
      stop("should not specify both start$fixef and nAGQ==0")
    }
  }
  if (control$nAGQ0initStep) {
    opt <- lme4::optimizeGlmer(devfun, optimizer = control$optimizer[[1]], 
      restart_edge = if (nAGQ == 0) 
        control$restart_edge
      else FALSE, boundary.tol = if (nAGQ == 0) 
        control$boundary.tol
      else 0, control = control$optCtrl, start = start, 
      nAGQ = 0, verbose = verbose, calc.derivs = FALSE)
  }
  if (nAGQ > 0L) {
    start <- lme4:::updateStart(start, theta = opt$par)
    devfun <- lme4::updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    if (devFunOnly) {
      return(devfun)
    }
    opt <- lme4::optimizeGlmer(devfun, optimizer = control$optimizer[[2]], 
      restart_edge = control$restart_edge, boundary.tol = control$boundary.tol, 
      control = control$optCtrl, start = start, nAGQ = nAGQ, 
      verbose = verbose, stage = 2, calc.derivs = control$calc.derivs, 
      use.last.params = control$use.last.params)
  }
  cc <- if (!control$calc.derivs) {
    NULL
  } else {
    if (verbose > 10) {
      cat("checking convergence\n")
    }
    lme4:::checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
      lbound = environment(devfun)$lower)
  }
  return(lme4::mkMerMod(environment(devfun), opt, glmod$reTrms,
    fr = nauf_off(glmod$fr), mc = mcout, lme4conv = cc))
}


nauf_glmer_nb <- function(..., interval = log(th) + c(-3, 3),
                     tol = 5e-5, verbose = FALSE, nb.control = NULL,
                     initCtrl = list(limit = 20, eps = 2 * tol, trace = verbose)) {
  dotE <- as.list(substitute(E(...))[-1])

  mc <- match.call()
  mc[[1]] <- quote(nauf_glmer)
  mc$family <- quote(stats::poisson)
  mc$verbose <- (verbose >= 2)
  g0 <- eval(mc, parent.frame(1L))

  th <- lme4:::est_theta(g0, limit = initCtrl$limit,
    eps = initCtrl$eps, trace = initCtrl$trace)
  
  if (verbose) cat("th := est_theta(glmer(..)) =", format(th))

  mc$family <- bquote(MASS::negative.binomial(theta=.(th)))
  g1 <- eval(mc, parent.frame(1L))

  if (verbose) cat(" --> dev.= -2*logLik(.) =", format(-2 * logLik(g1)),"\n")
  if ("data" %in% names(g1@call)) {
    if (!is.null(dotE[["data"]])) {
      g1@call[["data"]] <- dotE[["data"]]
    }
  } else {
    warning("no 'data = *' in glmer.nb() call ... Not much is guaranteed")
  }
  other.args <- c("verbose", "control")
  for (a in other.args) {
    if (a %in% names(g1@call)) {
      g1@call[[a]] <- dotE[[a]]
    }
  }
  
  return(lme4:::optTheta(g1, interval = interval, tol = tol, verbose = verbose,
    control = c(eval.parent(g1@call$control), nb.control)))
}

