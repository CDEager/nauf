

#' Fit a mixed effects regression using \code{nauf} contrasts.
#'
#' The mixed effects regression functions \code{nauf_lmer}, \code{nauf_glmer.nb},
#' and \code{nauf_glmer} fit linear, negative binomial, and other generalized
#' linear mixed effects models, respectively, impelementing
#' \code{\link{nauf_contrasts}}.
#'
#' \code{nauf_lmer}, \code{nauf_glmer}, and \code{nauf_glmer.nb} are based on
#' the \code{lme4} functions \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
#' and \code{\link[lme4]{glmer.nb}}, respectively, but implement
#' \code{\link{nauf_contrasts}}.  The \code{nauf} functions have all the same
#' arguments as the functions they are based on, but additionally
#' \code{ncs_scale}, which is passed to \code{\link{nauf_model.frame}}.  Other
#' than \code{ncs_scale}, the arguments have the same functions as they do in
#' the functions they are based on.  The default values for \code{na.action} and
#' \code{contrasts} cannot be changed (see \code{\link{nauf_model.frame}}).
#'
#' @param formula,data,subset,weights,offset,control,start,devFunOnly,verbose,REML,family,etastart,mustart,nAGQ,interval,tol,nb.control,initCtrl,... See
#'   \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}, and
#'   \code{\link[lme4]{glmer.nb}}.  Note that many arguments that are passed to
#'   \code{nauf_glmer} are passed to \code{nauf_glmer.nb} through \code{...},
#'   including \code{ncs_scale}.
#' @param na.action,contrasts Changes to the default values for these arguments
#'   are ignored with a warning.  See \code{\link{nauf_model.frame}}.
#' @param ncs_scale A positive number to be passed as the \code{scale} argument
#'   to \code{\link[standardize]{named_contr_sum}} for all unordered factors.
#'   See \code{\link{nauf_model.frame}}.  For \code{nauf_glmer.nb},
#'   \code{ncs_scale} is passed through \code{...}.
#'
#' @return A fitted model of class \code{\linkS4class{nauf.lmerMod}}
#'   (\code{nauf_lmer}) or \code{\linkS4class{nauf.glmerMod}} (\code{nauf_glmer}
#'   and \code{nauf_glmer.nb}).
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' sdat <- standardize(intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item), dat)
#' mod <- nauf_lmer(sdat$formula, sdat$data)
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors; and \code{\link[lme4]{lmer}},
#'   \code{\link[lme4]{glmer}}, and \code{\link[lme4]{glmer.nb}} for argument
#'   definitions.
#'
#' @export
nauf_glmer <- function(formula, data = NULL, family = gaussian,
                       control = lme4::glmerControl(), start = NULL,
                       verbose = 0L, nAGQ = 1L, subset, weights,
                       na.action = na.pass, offset, contrasts = NULL, mustart,
                       etastart, devFunOnly = FALSE,
                       ncs_scale = attr(formula, "standardized.scale"), ...) {
  # based on lme4::glmer
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

  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  mc$contrasts <- NULL
  mc$na.action <- na.pass

  if (is.linear(family <- get_family(family))) {
    warning("calling nauf_glmer with family=gaussian (identity link) as a ",
      "shortcut to lmer() is deprecated; please call nauf_lmer directly")
    mc[[1]] <- quote(nauf::nauf_lmer)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  } else if (is.character(family) && family == "negbin") {
    stop("To fit a negative binomial model, use nauf_glmer.nb")
  }

  mc[[1]] <- quote(nauf::nauf_glFormula)
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
  if (nAGQ == 0 && devFunOnly) return(devfun)

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
    start <- lme4_updateStart(start, theta = opt$par)
    devfun <- lme4::updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    if (devFunOnly) return(devfun)
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
    cc <- lme4_checkConv(attr(opt, "derivs"), opt$par,
      ctrl = control$checkConv, lbound = environment(devfun)$lower)
  }

  return(nauf.glmerMod(lme4::mkMerMod(environment(devfun), opt, glmod$reTrms,
    fr = glmod$fr, mc = mcout, lme4conv = cc)))
}


#' @rdname nauf_glmer
#' @export
nauf_lmer <- function(formula, data = NULL, REML = TRUE,
                      control = lme4::lmerControl(), start = NULL, verbose = 0L,
                      subset, weights, na.action = na.pass, offset,
                      contrasts = NULL, devFunOnly = FALSE,
                      ncs_scale = attr(formula, "standardized.scale"), ...) {
  # based on lme4::lmer
  mc <- mcout <- match.call()
  missCtrl <- missing(control)

  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  mc$contrasts <- NULL
  mc$na.action <- na.pass

  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control)) stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
      immediate. = TRUE)
    control <- do.call(lme4::lmerControl, control)
  }

  if (!is.null(list(...)[["family"]])) {
    warning("calling nauf_lmer with 'family' is deprecated; please use ",
      "nauf_glmer instead")
    mc[[1]] <- quote(nauf::nauf_glmer)
    if (missCtrl) mc$control <- lme4::glmerControl()
    return(eval(mc, parent.frame(1L)))
  }

  mc$control <- control
  mc[[1]] <- quote(nauf::nauf_lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  devfun <- do.call(lme4::mkLmerDevfun, c(lmod, list(start = start,
    verbose = verbose, control = control)))
  if (devFunOnly) return(devfun)

  if (control$optimizer == "none") {
    opt <- list(par = NA, fval = NA, conv = 1000, message = "no optimization")
  } else {
    opt <- lme4::optimizeLmer(devfun, optimizer = control$optimizer,
      restart_edge = control$restart_edge, boundary.tol = control$boundary.tol,
      control = control$optCtrl, verbose = verbose, start = start,
      calc.derivs = control$calc.derivs,
      use.last.params = control$use.last.params)
  }

  cc <- lme4_checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
    lbound = environment(devfun)$lower)

  return(nauf.lmerMod(lme4::mkMerMod(environment(devfun), opt, lmod$reTrms,
    fr = lmod$fr, mc = mcout, lme4conv = cc)))
}


#' @rdname nauf_glmer
#' @export
nauf_glmer.nb <- function (..., interval = log(th) + c(-3, 3), tol = 5e-05,
                           verbose = FALSE, nb.control = NULL,
                           initCtrl = list(limit = 20, eps = 2 * tol,
                           trace = verbose)) {
  # based on lme4::glmer.nb
  dotE <- as.list(substitute(E(...))[-1])
  mc <- match.call()
  mc[[1]] <- quote(nauf::nauf_glmer)
  mc$family <- quote(stats::poisson)
  mc$verbose <- (verbose >= 2)
  g0 <- eval(mc, parent.frame(1L))

  th <- lme4_est_theta(g0, limit = initCtrl$limit, eps = initCtrl$eps,
    trace = initCtrl$trace)
  if (verbose) cat("th := est_theta(glmer(..)) =", format(th))

  mc$family <- bquote(MASS::negative.binomial(theta = .(th)))
  g1 <- eval(mc, parent.frame(1L))
  if (verbose) cat(" --> dev.= -2*logLik(.) =", format(-2 * logLik(g1)), "\n")
  if ("data" %in% names(g1@call)) {
    if (!is.null(dotE[["data"]])) {
      g1@call[["data"]] <- dotE[["data"]]
    }
  } else {
    warning("no 'data = *' in nauf_glmer.nb() call ... Not much is guaranteed")
  }

  other.args <- c("verbose", "control")
  for (a in other.args) {
    if (a %in% names(g1@call)) {
      g1@call[[a]] <- dotE[[a]]
    }
  }

  return(nauf.glmerMod(lme4_optTheta(g1, interval = interval, tol = tol,
    verbose = verbose, control = c(eval.parent(g1@call$control), nb.control))))
}

