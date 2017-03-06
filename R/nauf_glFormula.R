

#' @export
nauf_glmer.nb <- function(..., interval = log(th) + c(-3, 3), tol = 5e-05,
                          verbose = FALSE, nb.control = NULL,
                          initCtrl = list(limit = 20, eps = 2 * tol,
                          trace = verbose)) {
  dotE <- as.list(substitute(E(...))[-1])
  mc <- match.call()
  
  mc[[1]] <- quote(nauf::nauf_glmer)
  mc$family <- quote(poisson)
  mc$verbose <- (verbose >= 2)
  g0 <- eval(mc, parent.frame(1L))
  th <- lme4:::est_theta(g0, limit = initCtrl$limit, eps = initCtrl$eps, 
    trace = initCtrl$trace)
  if (verbose) cat("th := est_theta(glmer(..)) =", format(th))
  
  mc$family <- bquote(MASS::negative.binomial(theta = .(th)))
  g1 <- eval(mc, parent.frame(1L))
  if (verbose) {
    cat(" --> dev.= -2*logLik(.) =",
      format(-2 * lme4:::logLik.merMod(g1)), "\n")
  }
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
  
  return(nauf.glmerMod(lme4:::optTheta(g1, interval = interval, tol = tol,
    verbose = verbose, control = c(eval.parent(g1@call$control), nb.control))))
}


#' @export
nauf_lFormula <- function(formula, data = NULL, REML = TRUE, subset, weights, 
                          na.action, offset, contrasts = NULL,
                          control = lme4::lmerControl(), ...) {
  control <- control$checkControl
  mf <- mc <- match.call()
  
  ignoreArgs <- c("start", "verbose", "devFunOnly", "control")
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
  formula[[length(formula)]] <- lme4::expandDoubleVerts(
    formula[[length(formula)]])
  mc$formula <- formula
  
  m <- match(c("data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf:::nauf_model.frame)
  fr <- eval(mf, parent.frame())
  
  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)
  
  reTrms <- nauf_mkReTrms(fr)
  n <- nrow(fr)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ", "please use ", 
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  
  return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula, 
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}


#' @export
nauf_glFormula <- function(formula, data = NULL, family = gaussian, subset,
                           weights, na.action, offset, contrasts = NULL,
                           mustart, etastart, control = lme4::glmerControl(),
                           ...) {
  control <- control$checkControl
  mf <- mc <- match.call()
  
  family <- get_family(family)
  if (isTRUE(all.equal(family, gaussian()))) {
    mc[[1]] <- quote(nauf:::nauf_lFormula)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  if (family$family %in% c("quasibinomial", "quasipoisson", "quasi")) {
    stop("\"quasi\" families cannot be used in glmer")
  }
  
  ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer",
    "control", "nAGQ")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4:::checkArgs, c(list("glmer"), l...))
  
  cstr <- "check.formula.LHS"
  lme4:::checkCtrlLevels(cstr, control[[cstr]])
  formula[[length(formula)]] <- lme4::expandDoubleVerts(
    formula[[length(formula)]])
  mc$formula <- formula
  
  m <- match(c("data", "subset", "weights", "na.action", "offset", "mustart",
    "etastart"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf:::nauf_model.frame)
  fr <- eval(mf, parent.frame())
  
  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4:::chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4:::checkScaleX(X, kind = scaleX.chk)
  
  reTrms <- nauf_mkReTrms(fr)
  n <- nrow(fr)
  wmsgNlev <- lme4:::checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4:::checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ", "please use ", 
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4:::checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  
  return(list(fr = fr, X = X, reTrms = reTrms, family = family,
    formula = formula, wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims,
    Zrank = wmsgZrank)))
}

