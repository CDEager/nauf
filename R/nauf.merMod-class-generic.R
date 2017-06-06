

#' Class for fitted mixed effects models with \code{nauf} contrasts.
#'
#' Models fit with \code{\link{nauf_lmer}} have class \code{nauf.lmerMod}
#' (inheriting from \code{\linkS4class{lmerMod}}) and models fit with
#' \code{\link{nauf_glmer}} and \code{\link{nauf_glmer.nb}} have class
#' \code{nauf.glmerMod} (inheriting from \code{\linkS4class{glmerMod}}).
#'
#' @slot rsp,Gp,call,frame,flist,cnms,lower,theta,beta,u,devcomp,pp,optinfo See
#'   \code{\linkS4class{merMod}}.
#'
#' @section Generic Methods:
#' There are S3 methods specific to the
#' \code{nauf.lmerMod} and \code{nauf.glmerMod} classes for the generic
#' functions \code{\link[=predict.nauf.merMod]{predict}} and
#' \code{\link[=anova.nauf.merMod]{anova}}, which behave differently from the
#' methods for \code{\linkS4class{merMod}} objects. Other methods for generic
#' functions from the \code{stats}, \code{MASS}, and \code{lme4} packages
#' should work as they would for \code{\linkS4class{merMod}} objects, with the
#' exception of the method for the \code{\link[stats]{simulate}} function,
#' for which the \code{nauf} method has more limited options than
#' \code{\link[lme4]{simulate.merMod}}, and is not meant for end user use
#' at this time (the method is necessary for the \code{PB} method to work
#' for \code{\link[=anova.nauf.merMod]{anova}}). If you encounter a generic
#' function in these packages which does not function properly, please report
#' the issue at \url{https://github.com/CDEager/nauf/issues}.
#'
#' @seealso \code{\link{nauf_glmer}}, \code{\link{nauf_contrasts}}, and
#'   \code{\linkS4class{merMod}}.
#'
#' @name nauf.merMod
NULL


#' @rdname nauf.merMod
nauf.lmerMod <- setClass("nauf.lmerMod", contains = "lmerMod")


#' @rdname nauf.merMod
nauf.glmerMod <- setClass("nauf.glmerMod", contains = "glmerMod")


formula.nauf.merMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                ...) {
  # based on lme4_formula.merMod
  if (missing(fixed.only) && random.only) fixed.only <- FALSE
  if (fixed.only && random.only) {
    stop("can't specify 'only fixed' and 'only random' terms")
  }

  if (is.null(form <- attr(x@frame, "formula"))) {
    if (!grepl("lmer$", deparse(getCall(x)[[1]]))) {
      stop("can't find formula stored in model frame or call")
    }
    form <- stats::as.formula(formula(getCall(x), ...))
  }

  if (fixed.only) form <- lme4_getFixedFormula(form)
  if (random.only) form <- lme4_reOnly(form, response = TRUE)

  first_class(form) <- "nauf.formula"
  attr(form, "standardized.scale") <- attr(attr(x@frame, "formula"),
    "standardized.scale")
  return(form)
}


#' @export
formula.nauf.lmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                 ...) {
  return(formula.nauf.merMod(x, fixed.only, random.only, ...))
}


#' @export
formula.nauf.glmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                  ...) {
  return(formula.nauf.merMod(x, fixed.only, random.only, ...))
}


terms.nauf.merMod <- function(x, fixed.only = TRUE, random.only = FALSE, ...) {
  # based on lme4_terms.merMod
  if (missing(fixed.only) && random.only) fixed.only <- FALSE
  if (fixed.only && random.only) {
    stop("can't specify 'only fixed' and 'only random' terms")
  }

  tt <- attr(x@frame, "terms")

  if (fixed.only) {
    tt <- terms.formula(formula(x, fixed.only = TRUE))
    attr(tt, "predvars") <- attr(terms(x@frame), "predvars.fixed")
  }
  if (random.only) {
    tt <- terms.formula(subbars(formula(x, random.only = TRUE)))
    attr(tt, "predvars") <- attr(terms(x@frame), "predvars.random")
  }

  first_class(tt) <- "nauf.terms"
  nauf.info(tt) <- nauf.info(x@frame)

  return(tt)
}


#' @export
terms.nauf.lmerMod <- function(x, fixed.only = TRUE, random.only = FALSE,
                               ...) {
  return(terms.nauf.merMod(x, fixed.only, random.only, ...))
}


#' @export
terms.nauf.glmerMod <- function(x, fixed.only = TRUE, random.only = FALSE,
                                ...) {
  return(terms.nauf.merMod(x, fixed.only, random.only, ...))
}


#' Predictions from a mixed effects \code{nauf} model at new data values.
#'
#' The \code{\link[stats]{predict}} method for
#' \code{\linkS4class{nauf.glmerMod}} and \code{\linkS4class{nauf.lmerMod}}
#' objects (the results of \code{\link{nauf_glmer}}, \code{\link{nauf_glmer.nb}},
#' and \code{\link{nauf_glmer}}).  It is based on
#' \code{\link[lme4]{predict.merMod}}, but currently some options are not
#' supported for \code{nauf} models.
#'
#' @param object A \code{\linkS4class{nauf.lmerMod}} or
#'   \code{\linkS4class{nauf.glmerMod}}.
#' @param newdata A data frame to make predictions on.
#' @param newparams,terms,allow.new.levels Changes to default values are not
#'   currently supported and result in an error.
#' @param re.form Formula for random effects to condition on.  Currently, only
#'   \code{NULL} (the default, indicating conditioning on \emph{all} random
#'   effects in the model) and \code{NA} or \code{~ 0} (indicating to use only
#'   the fixed effects in the predictions) are supported (i.e. you cannot
#'   currently condition on a subset of the random effects).
#' @param ReForm,REForm,REform Older versions of \code{re.form} in \code{lme4}
#'   which are now deprecated.
#' @param type Whether the predictions should be transformed with the inverse
#'   link function.
#' @param na.action Changes from default of \code{na.pass} are ignored with a
#'   warning.
#' @param ... Additional parameters (currently unused and ignored with a
#'   warning).
#'
#' @return A numeric vector of predicted values.
#'
#' @examples
#' \dontrun{
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' sdat <- standardize(intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item), dat)
#'
#' mod <- nauf_lmer(sdat$formula, sdat$data)
#' fit <- predict(mod)  # fitted values
#' preds <- predict(mod, sdat$data)  # predict same data using all ranef
#' preds_fe <- predict(mod, sdat$data, re.form = NA)  # only use fixef
#'
#' isTRUE(all.equal(fit, preds))  # TRUE
#' isTRUE(all.equal(preds, preds_fe))  # FALSE
#' }
#'
#' @seealso \code{\link[lme4]{predict.merMod}}, \code{\link{nauf_lmer}},
#'   \code{\link{nauf_glmer}}, \code{\link{nauf_glmer.nb}},
#'   \code{\linkS4class{nauf.lmerMod}}, and \code{\linkS4class{nauf.glmerMod}}.
#'
#' @name predict.nauf.merMod
NULL


predict.nauf.merMod <- function(object, newdata = NULL, newparams = NULL,
                                re.form = NULL, ReForm, REForm, REform,
                                terms = NULL, type = c("link", "response"),
                                allow.new.levels = FALSE, na.action = na.pass,
                                ...) {
  # based on lme4_predict.merMod but fewer options
  # notably cannot specify new ranef structure currently
  # and cannot have new levels or new paramms

  type <- match.arg(type)
  if (allow.new.levels) {
    stop("'allow.new.levels' not currently supported; must be FALSE")
  }
  if (!is.null(terms)) {
    stop("terms functionality for predict not yet implemented")
  }
  if (!is.null(newparams)) {
    stop("'newparams' not currently supported; must be NULL")
  }
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }
  re <- condition_on_re(re.form, ReForm, REForm, REform)
  if (length(list(...)) > 0) warning("unused arguments ignored")

  if ((is.null(newdata) && re)) {
    if (lme4::isLMM(object)) {
      pred <- stats::na.omit(fitted(object))
    } else {
      pred <- switch(type, response = object@resp$mu, link = object@resp$eta)
      if (is.null(nm <- rownames(model.frame(object)))) nm <- seq_along(pred)
      names(pred) <- nm
    }

  } else {
    X <- lme4::getME(object, "X")
    X.col.dropped <- attr(X, "col.dropped")
    if (is.null(newdata)) {
      offset <- model.offset(model.frame(object))
      if (is.null(offset)) offset <- 0

    } else {
      if (re) {
        Terms <- attr(object@frame, "terms")
      } else {
        Terms <- terms(object, fixed.only = TRUE)
      }
      mfnew <- suppressWarnings(model.frame(stats::delete.response(Terms),
        newdata))
      attr(mfnew, "formula") <- attr(object@frame, "formula")
      X <- nauf_mm(mfnew)
      offset <- 0
      tt <- terms(object)
      if (!is.null(off.num <- attr(tt, "offset"))) {
        for (i in off.num) {
          offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
        }
      }
      if (is.numeric(X.col.dropped) && length(X.col.dropped)) {
        X <- X[, -X.col.dropped, drop = FALSE]
      }
    }

    pred <- drop(X %*% lme4::fixef(object))
    pred <- pred + offset

    if (re) {
      z <- Matrix::t(nauf_mkReTrms(mfnew, levasgn(object))$Zt)
      b <- as.vector(lme4::getME(object, "b"))
      pred <- pred + base::drop(as(z %*% b, "matrix"))
    }

    if (lme4::isGLMM(object) && type == "response") {
      pred <- object@resp$family$linkinv(pred)
    }
  }

  return(pred)
}


#' @rdname predict.nauf.merMod
#' @export
predict.nauf.glmerMod <- function(object, newdata = NULL, newparams = NULL,
                                  re.form = NULL, ReForm, REForm, REform,
                                  terms = NULL, type = c("link", "response"),
                                  allow.new.levels = FALSE, na.action = na.pass,
                                  ...) {
  return(predict.nauf.merMod(object, newdata, newparams, re.form, ReForm, REForm,
    REform, terms, type, allow.new.levels, na.action, ...))
}


#' @rdname predict.nauf.merMod
#' @export
predict.nauf.lmerMod <- function(object, newdata = NULL, newparams = NULL,
                                 re.form = NULL, ReForm, REForm, REform,
                                 terms = NULL, type = c("link", "response"),
                                 allow.new.levels = FALSE, na.action = na.pass,
                                 ...) {
  return(predict.nauf.merMod(object, newdata, newparams, re.form, ReForm, REForm,
    REform, terms, type, allow.new.levels, na.action, ...))
}



#' Type III anovas for mixed effects \code{nauf} models.
#'
#' Obtain an anova table for a \code{\linkS4class{nauf.lmerMod}} or
#' \code{\linkS4class{nauf.glmerMod}} model.  Currently only Type III tests
#' are supported.
#'
#' There are six methods of p-value calculation which are supported:
#'
#' \describe{
#'   \item{lme4}{The default method. See \code{\link[lme4]{anova.merMod}}.}
#'   \item{S}{\emph{nauf.lmerMod models only}. Computes F-tests using the
#'     Satterthwaite approximation of denominator degrees of freedom,
#'     implemented with \code{\link[lmerTest]{calcSatterth}}.}
#'   \item{KR}{\emph{nauf.lmerMod models only}. If \code{object} was fit with
#'     maximum likelihood (\code{ML}), then the model is refit with restricted
#'     maximum likelihood (\code{REML}) first. Then computes F-tests using the
#'     Kenward-Roger approximation of denominator degrees of freedom,
#'     implemented with \code{\link[car]{Anova.merMod}}.}
#'   \item{nested-KR}{\emph{nauf.lmerMod models only}. If \code{object} was fit
#'     with maximum likelihood (\code{ML}), then the model is refit with
#'     restricted maximum likelihood (\code{REML}) first. Then for each fixed
#'     effects term, a restricted nested model is fit lacking only that fixed
#'     effects term, and F-tests are computed using the Kenward-Roger
#'     approximation of denominator degrees of freedom,
#'     implemented with \code{\link[pbkrtest]{KRmodcomp}}.  The full model and
#'     restricted models are returned along with the anova table, similar to
#'     \code{\link[afex]{mixed}}.}
#'   \item{LRT}{If \code{object} is a \code{\linkS4class{nauf.lmerMod}} fit with
#'     \code{REML}, it is first refit with \code{ML}.  Then restricted models
#'     are fit as in the \code{nested-KR} method, and Chi-squared (likelihood
#'     ratio) tests are computed.  The full model and restricted models are
#'     returned along with the anova table, similar to
#'     \code{\link[afex]{mixed}}.}
#'   \item{PB}{If \code{object} is a \code{\linkS4class{nauf.lmerMod}} fit with
#'     \code{REML}, it is first refit with \code{ML}. Then likelihood ratios are
#'     computed as for the \code{LRT} method, and p-values for the likelihood
#'     ratios are computed using parametric bootstrapping, implemented with
#'     \code{\link[pbkrtest]{PBmodcomp}}.  The full model and restricted models
#'     are returned along with the anova table, similar to
#'     \code{\link[afex]{mixed}}.}
#' }
#'
#' @param object A \code{\linkS4class{nauf.lmerMod}} or
#'   \code{\linkS4class{nauf.glmerMod}}.
#' @param ... Additional \code{nauf} models for the \code{lme4} method. See
#'   \code{\link[lme4]{anova.merMod}}.
#' @param refit For the \code{lme4} method, a logical indicating whether
#'   \code{\linkS4class{nauf.lmerMod}} models fit with \code{REML} should be
#'   refit with \code{ML} prior to comparison with models in \code{...}; default
#'   \code{TRUE}. See \code{\link[lme4]{anova.merMod}}.
#' @param model.names For the \code{lme4} method, character vectors of model
#'   names to be used in the anova table. See \code{\link[lme4]{anova.merMod}}.
#' @param method The method for calculating p-values.  See 'Details'.
#' @param test_intercept For all methods besides \code{lme4}, whether a test
#'   should be performed for the intercept term (default \code{FALSE}).
#' @param args_test For methods \code{nested-KR} and \code{PB}, an optional
#'   named list of arguments to be passed to \code{\link[pbkrtest]{KRmodcomp}}
#'   and \code{\link[pbkrtest]{PBmodcomp}}, respectively.
#'
#' @return The object returned depends on the method, and has class
#'   \code{nauf.mer.anova}.  For the \code{lme4},
#'   \code{S}, and \code{KR} methods, it is an \code{\link[stats]{anova}} table.
#'   For the \code{nested-KR}, \code{PB}, and \code{LRT} methods,
#'   a list with the anova table and restricted models is returned (similar
#'   to the output of \code{\link[afex]{mixed}}).
#'
#' @examples
#' dat <- droplevels(subset(plosives, voicing == "Voiceless"))
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' sdat <- standardize(cdur ~ dialect * spont + (1 | speaker) + (1 | item), dat)
#'
#' mod <- nauf_lmer(sdat$formula, sdat$data)
#'
#' \dontrun{
#' # lme4 method anova table
#' anova(mod)
#'
#' # anova table using Satterthwaite approximation
#' anova(mod, method = "S")
#'
#' # anova table using Kenward-Roger approximation
#' anova(mod, method = "KR")
#'
#' # list with restricted models and Kenward-Roger table
#' anova(mod, method = "nested-KR")
#'
#' # list with restricted models and parametric bootstrap table
#' # model is first refit with maximum likelihood
#' anova(mod, method = "PB")
#'
#' # list with restricted models and likelihood ratio test table
#' # model is first refit with maximum likelihood
#' anova(mod, method = "LRT")
#' }
#'
#' @seealso \code{\linkS4class{nauf.lmerMod}} and
#'   \code{\linkS4class{nauf.glmerMod}} classes;
#'   \code{\link[lme4]{anova.merMod}} for the \code{lme4} method;
#'   \code{\link[car]{Anova.merMod}} for the \code{KR} method;
#'   \code{\link[lmerTest]{calcSatterth}} for the \code{S} method;
#'   \code{\link[afex]{mixed}} for the \code{nested-KR}, \code{LRT}, and
#'   \code{PB} methods; \code{\link[pbkrtest]{KRmodcomp}} for the
#'   \code{nested-KR} method; \code{\link[pbkrtest]{PBmodcomp}} for the
#'   \code{PB} method.
#'
#' @name anova.nauf.merMod
NULL


anova.nauf.merMod <- function(object, ..., refit = TRUE, model.names = NULL,
                               method = c("lme4", "S", "KR", "LRT", "PB", "nested-KR"),
                              test_intercept = FALSE, args_test = NULL) {
  # drawing from afex::mixed
  mc <- match.call()
  method <- match.arg(method)
  dots <- list(...)
  linear <- lme4::isLMM(object)
  
  if (method %in% c("S", "KR", "nested-KR") && !linear) {
    stop("Method ", method, " is only valid for linear models")
  }
  if ("type" %in% names(dots) && !(dots$type %in% c("3", "III"))) {
    stop("Currently only Type III tests are supported")
  }
  if (dotmods <- length(dots)) dotmods <- sum(sapply(dots, is.nauf.model))
  
  if (method == "lme4" || dotmods) {
    if (method != "lme4") {
      stop("Multiple nauf models supplied. Must use method 'lme4'")
    }
    out <- NextMethod("anova")
    first_class(out) <- "nauf.mer.anova"
    return(out)
  }

  reml <- method %in% c("KR", "nested-KR")
  if (method != "S" && reml != lme4::isREML(object)) {
    call <- getCall(object)
    call[["REML"]] <- reml
    cat("Refitting model with REML =", reml, "since method =", method, "\n")
    object <- eval(call, parent.frame())
  }

  if (method == "KR") {
    out <- car::Anova(object, type = 3, test.statistic = "F")
    if (!test_intercept) out <- out[-1, , drop = FALSE]
    first_class(out) <- "nauf.mer.anova"
    return(out)
  }
  
  fenr <- stats::delete.response(terms(object))
  fe <- c(if (test_intercept) "(Intercept)", attr(fenr, "term.labels"))
  a1 <- as.integer(!test_intercept)
  asgn <- attr(X <- lme4::getME(object, "X"), "assign")
  
  if (method == "S") {
    L <- lapply(msplit(diag(length(asgn)), asgn)[paste(seq_along(fe) + a1 - 1)],
      drop)
    out <- data.frame(t(sapply(L, lmerTest::calcSatterth,
      model = object))[, c(4, 1:3), drop = FALSE])
    dimnames(out) <- list(fe, c("num Df", "den Df", "F", "Pr(>F)"))
    
    return(structure(out, class = c("nauf.mer.anova", "anova", "data.frame"),
      heading = c(paste0("Mixed Model Anova Table (Type III tests, ", method,
      "-method)\n"), paste0("Model: ", deparse(getCall(object)$formula)),
      paste0("Data: ", getCall(object)[["data"]])),
      sig_symbols = c(" +", " *", " **", " ***")))
  }
  
  mcout <- getCall(object)
  glmod <- list(fr = object@frame, X = X, reTrms = nauf_mkReTrms(object@frame))
  if (linear) {
    glmod$REML <- reml
  } else {
    glmod$family <- get_family(object)
  }
  
  getargs <- c("start", "verbose", "control", if (!linear) "nAGQ")
  defaults <- formals(if (linear) nauf_lmer else nauf_glmer)[getargs]
  vals <- list()
  for (a in getargs) {
    if (a %in% names(mcout)) {
      vals[[a]] <- eval(mcout[[a]], parent.frame())
    } else {
      vals[[a]] <- eval(defaults[[a]], parent.frame())
    }
  }
  
  fits <- do.call(fit_nested_models,
    c(rsa_nlist(fe, a1, glmod, X, asgn, mcout), vals))
  
  lik_full <- as.numeric(logLik(object))
  lik_rest <- as.numeric(vapply(fits, logLik, 0))
  if (length(better <- which(lik_rest > lik_full))) {
    warning(paste0("Following nested model(s) ",
      "provide better fit than full model:\n  ", add_quotes(fe[better])))
  }

  tests <- setNames(vector("list", length(fe)), fe)
  
  
  cat("Obtaining", length(fe), "p-values [")
  if (method == "nested-KR") {
    for (ft in 1:length(fe)) {
      tests[[ft]] <- do.call(pbkrtest::KRmodcomp,
        args = c(largeModel = object, smallModel = fits[[ft]], args_test))
      cat(".")
    }

    anova_table <- data.frame(t(vapply(tests, function(x)
      unlist(x[["test"]][1, ]), unlist(tests[[1]][["test"]][1, ]))))
    dimnames(anova_table) <- list(fe, c("F", "num Df", "den Df", "F.scaling",
      "Pr(>F)"))
    anova_table <- anova_table[, c("num Df", "den Df", "F.scaling", "F",
      "Pr(>F)")]
    anova_tab_addition <- NULL

  } else if (method == "PB") {
    for (ft in 1:length(fe)) {
      tests[[ft]] <- do.call(pbkrtest::PBmodcomp,
        args = c(largeModel = object, smallModel = fits[[ft]], args_test))
      cat(".")
    }

    anova_table <- data.frame(t(vapply(tests, function(x)
      unlist(x[["test"]][2, ]), unlist(tests[[1]][["test"]][2, ]))))
    anova_table <- anova_table[, -2]
    LRT <- vapply(tests, function(x) unlist(x[["test"]][1, ]),
      unlist(tests[[1]][["test"]][1, ]))
    row.names(LRT) <- stringr::str_c(row.names(LRT), ".LRT")
    anova_table <- cbind(anova_table, t(LRT))
    rownames(anova_table) <- fe
    anova_table <- anova_table[, c("stat", "df.LRT", "p.value.LRT",
      "p.value")]
    colnames(anova_table) <- c("Chisq", "Chi Df", "Pr(>Chisq)", "Pr(>PB)")
    if (!is.null(args_test$nsim)) {
      nsim <- args_test$nsim
    } else {
      nsim <- formals(pbkrtest::PBmodcomp)$nsim
    }
    anova_tab_addition <- paste("Based on", nsim, "simulations")

  } else if (method == "LRT") {
    for (ft in 1:length(fe)) {
      tests[[ft]] <- anova(object, fits[[ft]], model.names = paste(1:2))
      cat(".")
    }
  
    anova_table <- data.frame(
      Df = vapply(tests, function(x) x[["Df"]][1], 0),
      Chisq = vapply(tests, function(x) x[["Chisq"]][2], 0),
      `Chi df` = vapply(tests, function(x) x[["Chi Df"]][2], 0),
      `Pr(>Chisq)` = vapply(tests, function(x) x[["Pr(>Chisq)"]][2], 0),
      stringsAsFactors = FALSE, check.names = FALSE, row.names = fe
    )
    anova_tab_addition <- paste0("Df full model: ", tests[[1]][["Df"]][2])
  }
  cat("]\n")

  anova_table <- structure(anova_table, class = c("anova", "data.frame"),
    heading = c(paste0("Mixed Model Anova Table (Type III tests, ",
      method, "-method)\n"), paste0("Model: ", deparse(getCall(object)$formula)),
      paste0("Data: ", getCall(object)[["data"]]), anova_tab_addition),
    sig_symbols = c(" +", " *", " **", " ***"))

  return(structure(list(anova_table = anova_table, full_model = object,
    restricted_models = fits, tests = tests), class = "nauf.mer.anova",
    method = method, type = "III"))
}


fit_nested_models <- function(fe, a1, glmod, X, asgn, mcout, start, verbose,
                              control, ...) {
  UseMethod("fit_nested_models")
}


fit_nested_models.nauf.lmerMod <- function(fe, a1, glmod, X, asgn, mcout, start,
                                           verbose, control, ...) {
  fits <- setNames(vector("list", length(fe)), fe)
  
  if (!inherits(control, "lmerControl")) {
    control <- do.call(lme4::lmerControl, control)
  }
  
  formula <- deparse(mcout[["formula"]])
  devfun_args <- c(glmod, rsa_nlist(start, verbose, control))
  opt_args <- c(list(devfun = NULL), control[c("optimizer", "restart_edge",
    "boundary.tol")], rsa_nlist(control = control$optCtrl, verbose, start),
    control[c("calc.derivs", "use.last.params")])
  
  cat("Fitting", length(fe), "nested nauf_lmer models [")
  for (ft in 1:length(fe)) {
    mcout[["formula"]] <- stats::formula(paste(formula, "-", fe[ft]))
    
    devfun_args$X <- X[, -which(asgn == (ft - 1 + a1)), drop = FALSE]
    opt_args$devfun <- do.call(lme4::mkLmerDevfun, devfun_args)
    
    opt <- do.call(lme4::optimizeLmer, opt_args)
    cc <- lme4_checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
      lbound = environment(opt_args$devfun)$lower)
      
    fits[[ft]] <- nauf.lmerMod(lme4::mkMerMod(environment(opt_args$devfun),
      opt, glmod$reTrms, fr = glmod$fr, mc = mcout, lme4conv = cc))
      
    cat(".")
  }
  cat("]\n")
  
  return(fits)
}


fit_nested_models.nauf.glmerMod <- function(fe, a1, glmod, X, asgn, mcout,
                                            start, verbose, control, nAGQ,
                                            ...) {
  fits <- setNames(vector("list", length(fe)), fe)
  
  if (!inherits(control, "glmerControl")) {
    control <- do.call(lme4::glmerControl, control)
  }
  if (control$nAGQ0initStep) {
    nAGQinit <- 0L
  } else {
    nAGQinit <- 1L
  }
  
  oc.start <- start
  fe.start <- is.list(start) && ("fixef" %in% names(start))
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
  
  formula <- deparse(mcout[["formula"]])
  devfun_args <- c(glmod, rsa_nlist(verbose, control, nAGQ = nAGQinit))
  opt_args1 <- list(optimizer = control$optimizer[[1]],
    restart_edge = if (nAGQ == 0) control$restart_edge else FALSE,
    boundary.tol = if (nAGQ == 0) control$boundary.tol else 0,
    control = control$optCtrl, start = NULL, nAGQ = 0, verbose = verbose,
    calc.derivs = FALSE)
  opt_args2 <- list(optimizer = control$optimizer[[2]],
    restart_edge = control$restart_edge, boundary.tol = control$boundary.tol,
    control = control$optCtrl, start = NULL, nAGQ = nAGQ, verbose = verbose,
    stage = 2, calc.derivs = control$calc.derivs,
    use.last.params = control$use.last.params)
  
  
  cat("Fitting", length(fe), "nested nauf_glmer models [")
  for (ft in 1:length(fe)) {
    mcout[["formula"]] <- stats::formula(paste(formula, "-", fe[ft]))
    devfun_args$X <- X[, -which(asgn == (ft - 1 + a1)), drop = FALSE]
    devfun <- do.call(lme4::mkGlmerDevfun, devfun_args)
    
    if (fe.start) {
      start <- oc.start
      start$fixef <- start[-which(asgn == (ft - 1 + a1))]
    }
    
    if (control$nAGQ0initStep) {
      opt_args1$start <- start
      opt <- do.call(lme4::optimizeGlmer, c(list(devfun = devfun), opt_args1))
      if (nAGQ > 0L) start <- lme4_updateStart(start, theta = opt$par)
    }
    
    if (nAGQ > 0L) {
      opt_args2$start <- start
      devfun <- lme4::updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
      opt <- do.call(lme4::optimizeGlmer, c(list(devfun = devfun), opt_args2))
    }

    if (!control$calc.derivs) {
      cc <- NULL
    } else {
      if (verbose > 10) cat("checking convergence\n")
      cc <- lme4_checkConv(attr(opt, "derivs"), opt$par,
        ctrl = control$checkConv, lbound = environment(devfun)$lower)
    }
      
    fits[[ft]] <- nauf.glmerMod(lme4::mkMerMod(environment(devfun), opt,
      glmod$reTrms, fr = glmod$fr, mc = mcout, lme4conv = cc))
      
    cat(".")
  }
  cat("]\n")
  
  return(fits)
}


#' @rdname anova.nauf.merMod
#'
#' @importFrom pbkrtest KRmodcomp PBmodcomp
#' @importFrom lmerTest calcSatterth
#' @importFrom car Anova
#'
#' @export
anova.nauf.lmerMod <- function(object, ..., refit = TRUE, model.names = NULL,
                               method = c("lme4", "S", "KR", "LRT", "PB", "nested-KR"),
                               test_intercept = FALSE, args_test = NULL) {
  method <- match.arg(method)
  return(anova.nauf.merMod(object, ..., refit = refit, model.names = model.names,
    method = method, test_intercept = test_intercept, args_test = args_test))
}


#' @rdname anova.nauf.merMod
#'
#' @importFrom pbkrtest PBmodcomp
#'
#' @export
anova.nauf.glmerMod <- function(object, ..., refit = TRUE, model.names = NULL,
                                method = c("lme4", "LRT", "PB"),
                                test_intercept = FALSE, args_test = NULL) {
  method <- match.arg(method)
  return(anova.nauf.merMod(object, ..., refit = refit, model.names = model.names,
    method = method, test_intercept = test_intercept, args_test = args_test))
}


#' @export
print.nauf.mer.anova <- function(x, ...) {
  if (!inherits(x, "list")) {
    class(x) <- class(x)[-1]
    print(x, ...)
    invisible(NULL)

  } else {
    # based on afex:::print.mixed
    lik_full <- as.numeric(logLik(x[["full_model"]]))
    lik_rest <- as.numeric(vapply(x[["restricted_models"]], logLik, 0))
    if (length(better <- which(lik_rest > lik_full))) {
      warning(paste0("Following nested model(s) ",
        "provide better fit than full model:\n  ",
        add_quotes(names(x[["restricted_models"]])[better])))
    }
    afex_get_mixed_warnings(x)
    print(x$anova_table)
    invisible(NULL)
  }
}


simulate.nauf.merMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                                  re.form = NA, ReForm, REForm, REform,
                                  newdata = NULL, newparams = NULL,
                                  family = NULL, allow.new.levels = FALSE,
                                  na.action = na.pass, ...) {
  dots <- list(...)
  if (is.null(dots$weights)) {
    if (is.null(newdata)) {
      weights <- object@resp$weights
    } else {
      weights <- rep(1, nrow(newdata))
    }
  }

  if (missing(object)) {
    stop("Currently, method where 'object' is missing is not supported")
  }
  stopifnot((nsim <- as.integer(nsim[1])) > 0, is(object, "merMod"))
  if (!is.null(newparams)) {
    stop("'newparams' not currently supported")
  }

  re.form.miss <- missing(re.form)
  re.form <- lme4_reFormHack(re.form, ReForm, REForm, REform)
  if (!missing(use.u)) {
    if (!re.form.miss) {
      stop("should specify only one of ", sQuote("use.u"),
        " and ", sQuote("re.form"))
    }
    if (use.u) {
      re.form <- NULL
    } else {
      re.form <- NA
    }
  }

  if (sim_new_re <- !is.null(re.form)) {
    re.form <- NA
  }
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  RNGstate <- .Random.seed
  sigma <- sigma(object)

  etapred <- predict(object, newdata = newdata, re.form = re.form, type = "link")
  n <- length(etapred)

  if (sim_new_re) {
    if (is.null(newdata)) {
      newRE <- nauf_mkReTrms(object@frame)
    } else {
      Terms <- attr(object@frame, "terms")
      mfnew <- suppressWarnings(model.frame(stats::delete.response(Terms),
        newdata))
      attr(mfnew, "formula") <- attr(object@frame, "formula")
      lvs <- lapply(object@flist, levels)
      newRE <- nauf_mkReTrms(mfnew, lvs)
    }
    U <- Matrix::t(newRE$Lambdat %*% newRE$Zt)
    u <- rnorm(ncol(U) * nsim)
    sim.reff <- as(U %*% matrix(u, ncol = nsim), "matrix")
  } else {
    sim.reff <- 0
  }

  if (lme4::isLMM(object)) {
    val <- etapred + sigma * (sim.reff + matrix(rnorm(n * nsim), ncol = nsim))

  } else if (lme4::isGLMM(object)) {
    etasim <- etapred + sim.reff
    family <- object@resp$family
    if (grepl("^Negative ?Binomial", family$family, ignore.case = TRUE)) {
      family$family <- "negative.binomial"
    }
    musim <- family$linkinv(etasim)

    if (is.null(sfun <- lme4_simfunList()[[family$family]]) &&
    is.null(family$simulate)) {
      stop("simulation not implemented for family", family$family)
    }

    val <- sfun(object, nsim = 1, ftd = rep_len(musim, n * nsim), wts = weights)

    if (family$family == "binomial" && is.matrix(r <- model.response(
    object@frame))) {
      val <- lapply(split(val[[1]], gl(nsim, n, 2 * nsim * n)),
        matrix, ncol = 2, dimnames = list(NULL, colnames(r)))

    } else if (family$family == "binomial" && is.factor(val[[1]])) {
      val <- split(val[[1]], gl(nsim, n))

    } else {
      val <- split(val, gl(nsim, n))
    }

  } else {
    stop("simulate method for NLMMs not yet implemented")
  }

  if (!is.list(val)) {
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
  } else {
    class(val) <- "data.frame"
  }
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  f <- fitted(object)
  nm <- names(f)[!is.na(f)]
  if (length(nm) == 0) {
    nm <- as.character(seq(n))
  } else if (!is.null(newdata)) {
    nm <- rownames(newdata)
  }
  row.names(val) <- nm

  return(structure(val, na.action = na.pass, seed = RNGstate))
}


#' @export
simulate.nauf.lmerMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                                  re.form = NA, ReForm, REForm, REform,
                                  newdata = NULL, newparams = NULL,
                                  family = NULL, allow.new.levels = FALSE,
                                  na.action = na.pass, ...) {
  return(simulate.nauf.merMod(object, nsim = nsim, seed = seed, use.u = use.u,
    re.form = re.form, ReForm = ReForm, REForm = REForm, REform = REform,
    newdata = newdata, newparams = newparams, family = family,
    allow.new.levels = allow.new.levels, na.action = na.action, ...))
}


#' @export
simulate.nauf.glmerMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                                   re.form = NA, ReForm, REForm, REform,
                                   newdata = NULL, newparams = NULL,
                                   family = NULL, allow.new.levels = FALSE,
                                   na.action = na.pass, ...) {
  return(simulate.nauf.merMod(object, nsim = nsim, seed = seed, use.u = use.u,
    re.form = re.form, ReForm = ReForm, REForm = REForm, REform = REform,
    newdata = newdata, newparams = newparams, family = family,
    allow.new.levels = allow.new.levels, na.action = na.action, ...))
}

