

#' @export
nauf_lFormula <- function(formula, data = NULL, REML = TRUE, subset, weights, 
                          na.action = na.pass, offset, contrasts = NULL,
                          control = lme4::lmerControl(), ncs_scale = NULL,
                          ...) {
  # based on lme4::lFormula
  control <- control$checkControl
  mf <- mc <- match.call()
  
  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }
  
  ignoreArgs <- c("start", "verbose", "devFunOnly", "control")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4_checkArgs, c(list("lmer"), l...))
  
  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(nauf::nauf_glFormula)
    if (missing(control)) mc[["control"]] <- lme4::glmerControl()
    return(eval(mc, parent.frame()))
  }
  
  cstr <- "check.formula.LHS"
  lme4_checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4_checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  formula <- stats::as.formula(formula, env = denv)
  lme4_RHSForm(formula) <- lme4::expandDoubleVerts(lme4_RHSForm(formula))
  mc$formula <- formula
  m <- match(c("data", "subset", "weights", "offset", "ncs_scale"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  
  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- formula
  fr <- eval(mf, parent.frame())
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)
  
  reTrms <- nauf_mkReTrms(fr)
  wmsgNlev <- lme4_checkNlevels(reTrms$flist, n = n, control)
  wmsgZdims <- lme4_checkZdims(reTrms$Ztlist, n = n, control, allow.n = FALSE)
  if (anyNA(reTrms$Zt)) {
    stop("NA in Z (random-effects model matrix): ", "please use ", 
      shQuote("na.action='na.omit'"), " or ", shQuote("na.action='na.exclude'"))
  }
  wmsgZrank <- lme4_checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06)
  
  mf[[1L]] <- quote(stats::model.frame)
  fixedform <- formula
  lme4_RHSForm(fixedform) <- lme4::nobars(lme4_RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr,
    "terms"), "predvars")
    
  ranform <- formula
  lme4_RHSForm(ranform) <- lme4::subbars(lme4_RHSForm(lme4_reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), 
    "predvars")
  
  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4_chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4_checkScaleX(X, kind = scaleX.chk)
  
  return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula, 
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}


#' @export
nauf_glFormula <- function(formula, data = NULL, family = gaussian, subset,
                           weights, na.action = na.pass, offset,
                           contrasts = NULL, mustart, etastart,
                           control = lme4::glmerControl(), ncs_scale = NULL,
                           ...) {
  # based on lme4::glFormula
  control <- control$checkControl
  mf <- mc <- match.call()
  
  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }
  
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if (is.function(family)) {
    family <- family()
  }
  if (isTRUE(all.equal(family, stats::gaussian()))) {
    mc[[1]] <- quote(nauf::nauf_lFormula)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  if (family$family %in% c("quasibinomial", "quasipoisson", "quasi")) {
    stop("\"quasi\" families cannot be used in glmer")
  }
  
  ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer", "control",
    "nAGQ")
  l... <- list(...)
  l... <- l...[!names(l...) %in% ignoreArgs]
  do.call(lme4_checkArgs, c(list("glmer"), l...))

  cstr <- "check.formula.LHS"
  lme4_checkCtrlLevels(cstr, control[[cstr]])
  denv <- lme4_checkFormulaData(formula, data,
    checkLHS = control$check.formula.LHS == "stop")
  mc$formula <- formula <- stats::as.formula(formula, env = denv)
  m <- match(c("data", "subset", "weights", "offset", "ncs_scale",
    "mustart", "etastart"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(nauf::nauf_model.frame)
    
  fr.form <- lme4::subbars(formula)
  environment(fr.form) <- environment(formula)
  for (i in c("weights", "offset")) {
    if (!eval(bquote(missing(x = .(i))))) {
      assign(i, get(i, parent.frame()), environment(fr.form))
    }
  }
  mf$formula <- formula
  fr <- eval(mf, parent.frame())
  attr(fr, "formula") <- formula
  attr(fr, "offset") <- mf$offset
  n <- nrow(fr)
  
  reTrms <- nauf_mkReTrms(fr)
  wmsgNlev <- lme4_checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
  wmsgZdims <- lme4_checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
  wmsgZrank <- lme4_checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06, 
    allow.n = TRUE)
        
  mf[[1L]] <- quote(stats::model.frame)
  fixedform <- formula
  lme4_RHSForm(fixedform) <- lme4::nobars(lme4_RHSForm(fixedform))
  mf$formula <- fixedform
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr,
    "terms"), "predvars")

  ranform <- formula
  lme4_RHSForm(ranform) <- lme4::subbars(lme4_RHSForm(lme4_reOnly(formula)))
  mf$formula <- ranform
  ranfr <- eval(mf, parent.frame())
  attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), "predvars")

  X <- nauf_model.matrix(fr)
  if (is.null(rankX.chk <- control[["check.rankX"]])) {
    rankX.chk <- eval(formals(lme4::lmerControl)[["check.rankX"]])[[1]]
  }
  X <- lme4_chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
  if (is.null(scaleX.chk <- control[["check.scaleX"]])) {
    scaleX.chk <- eval(formals(lme4::lmerControl)[["check.scaleX"]])[[1]]
  }
  X <- lme4_checkScaleX(X, kind = scaleX.chk)
  
  return(list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula, 
    wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)))
}

