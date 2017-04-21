

#' Create a model frame and fixed and random effects model matrices using \code{nauf} contrasts.
#'
#' The same as the \code{lme4} \code{\link[lme4]{modular}} functions
#' \code{glFormula} and \code{lFormula}, but implementing
#' \code{\link{nauf_contrasts}}.  \code{nauf_lFormula} is used for linear mixed
#' effects regressions (i.e. those that would be fit with
#' \code{\link{nauf_lmer}}) and \code{nauf_glFormula} is used for genarlized
#' linear mixed effects regressions (i.e. those that would be fit with
#' \code{\link{nauf_glmer}} or \code{\link{nauf_glmer.nb}}).  Both of the
#' functions contain a call to \code{nauf_mkReTrms}, which serves the same
#' purpose as the \code{lme4} function \code{\link[lme4]{mkReTrms}}, but with
#' \code{\link{nauf_contrasts}}, and, while \code{\link[lme4]{mkReTrms}} is
#' exported by \code{lme4}, \code{nauf_mkReTrms} is an internal function in the
#' \code{nauf} package.
#'
#' @param formula,data,REML,subset,weights,offset,control,mustart,etastart,...
#'   See \code{\link[lme4]{glFormula}}.
#' @param na.action,contrasts Changes from default values are ignored.  See
#'   \code{\link{nauf_model.frame}}.
#' @param ncs_scale A positive number to be passed as the \code{scale} argument
#'   to \code{\link[standardize]{named_contr_sum}} for all unordered factors.
#'   See \code{\link{nauf_model.frame}}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{fr}{The model frame (with class \code{\linkS4class{nauf.frame}}).
#'     See \code{\link{nauf_model.frame}}.}
#'   \item{X}{The fixed effects model matrix with \code{\link{nauf_contrasts}}
#'     applied. See \code{\link{nauf_model.matix}}.}
#'   \item{reTrms}{A list containing the random effects model matrix and other
#'     information about the random effects structure.  The elements of the list
#'     have the same structure as that returned by \code{\link[lme4]{mkReTrms}},
#'     but incorportating \code{\link{nauf_contrasts}}.}
#'   \item{REML}{(\code{nauf_lFormula} only): A logical indicating if restricted
#'     maximum likelihood was used (Copy of argument).}
#' }
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' dat_form <- intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item)
#' sdat <- standardize(dat_form, dat)
#' lmod <- nauf_lFormula(sdat$formula, sdat$data)
#'
#' vless <- droplevels(subset(dat, voicing == "Voiceless"))
#' vless$fully_voiced <- vless$vdur == 0
#' vless_form <- fully_voiced ~ dialect * spont +
#'   (1 + spont | speaker) + (1 + dialect | item)
#' svless <- standardize(vless_form, vless, family = binomial)
#' glmod <- nauf_glFormula(svless$formula, svless$data, family = binomial)
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors; \code{\link{nauf_model.frame}} and
#'   \code{\link{nauf_model.matrix}} for the creation of the \code{fr} and
#'   \code{X} elements of the returned list, respectively; and
#'   \code{\link{nauf_lmer}}, \code{\link{nauf_glmer.nb}}, and
#'   \code{\link{nauf_glmer}} for fitting mixed effects regressions with gaussian,
#'   negative binomial, and all other families, respectively.
#'
#' @export
nauf_glFormula <- function(formula, data = NULL, family = gaussian, subset,
                           weights, na.action = na.pass, offset,
                           contrasts = NULL, mustart, etastart,
                           control = lme4::glmerControl(),
                           ncs_scale = attr(formula, "standardized.scale"),
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


#' @rdname nauf_glFormula
#' @export
nauf_lFormula <- function(formula, data = NULL, REML = TRUE, subset, weights,
                          na.action = na.pass, offset, contrasts = NULL,
                          control = lme4::lmerControl(),
                          ncs_scale = attr(formula, "standardized.scale"),
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

