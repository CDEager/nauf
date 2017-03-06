

#' @export
nauf_lmer <- function(formula, data = NULL, REML = TRUE,
                      control = lme4::lmerControl(), start = NULL, verbose = 0L,
                      subset, weights, na.action, offset, contrasts = NULL,
                      devFunOnly = FALSE, ...) {
  mc <- mcout <- match.call()
  
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control)) {
      stop("'control' is not a list; use lmerControl()")
    }
    warning("passing control as list is deprecated: please use lmerControl() instead", 
      immediate. = TRUE)
    control <- do.call(lmerControl, control)
  }
  
  if (!is.null(list(...)[["family"]])) {
    mc[[1]] <- quote(nauf::nauf_glmer)
    if (missCtrl) {
      mc$control <- lme4::glmerControl()
    }
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
  
  if (control$optimizer == "none") {
    opt <- list(par = NA, fval = NA, conv = 1000, message = "no optimization")
  } else {
    opt <- lme4::optimizeLmer(devfun, optimizer = control$optimizer,
      restart_edge = control$restart_edge, boundary.tol = control$boundary.tol,
      control = control$optCtrl, verbose = verbose, start = start,
      calc.derivs = control$calc.derivs,
      use.last.params = control$use.last.params)
  }
  
  cc <- lme4:::checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
    lbound = environment(devfun)$lower)
  
  return(nauf.lmerMod(lme4::mkMerMod(environment(devfun), opt, lmod$reTrms,
    fr = lmod$fr, mc = mcout, lme4conv = cc)))
}


#' @export
nauf_glmer <- function(formula, data = NULL, family = gaussian,
                       control = lme4::glmerControl(), start = NULL,
                       verbose = 0L, nAGQ = 1L, subset, weights, na.action,
                       offset, contrasts = NULL, mustart, etastart,
                       devFunOnly = FALSE, ...) {
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
  family <- get_family(family)
  if (isTRUE(all.equal(family, gaussian()))) {
    mc[[1]] <- quote(nauf::nauf_lmer)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  
  mc[[1]] <- quote(nauf:::nauf_glFormula)
  glmod <- eval(mc, parent.frame(1L))
  mcout$formula <- glmod$formula
  glmod$formula <- NULL
  
  if (control$nAGQ0initStep) {
    nAGQinit <- 0L
  } else {
    nAGQinit <- 1L
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
      restart_edge = if (nAGQ == 0) control$restart_edge else FALSE,
      boundary.tol = if (nAGQ == 0) control$boundary.tol else 0,
      control = control$optCtrl, start = start, nAGQ = 0, verbose = verbose,
      calc.derivs = FALSE)
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
  
  if (!control$calc.derivs) {
    cc <- NULL
  } else {
    if (verbose > 10) cat("checking convergence\n")
    cc <- lme4:::checkConv(attr(opt, "derivs"), opt$par,
      ctrl = control$checkConv, lbound = environment(devfun)$lower)
  }
  
  return(nauf.glmerMod(mkMerMod(environment(devfun), opt, glmod$reTrms,
    fr = glmod$fr, mc = mcout, lme4conv = cc)))
}

