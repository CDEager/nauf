

###### is ######

#' @export
is.nauf.formula <- function(object) {
  return(inherits(object, "nauf.formula"))
}


#' @export
is.nauf.terms <- function(object) {
  return(inherits(object, "nauf.terms"))
}


#' @export
is.nauf.frame <- function(object) {
  return(inherits(object, "nauf.frame"))
}


#' @export
is.nauf.glm <- function(object) {
  return(inherits(object, "nauf.glm"))
}


#' @export
is.nauf.lmerMod <- function(object) {
  return(inherits(object, "nauf.lmerMod"))
}


#' @export
is.nauf.glmerMod <- function(object) {
  return(inherits(object, "nauf.glmerMod"))
}


###### formula ######

formula_nauf_merMod <- function(x, fixed.only = FALSE, random.only = FALSE,
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
  
  class(form) <- c("nauf.formula", "formula")
  
  return(form)
}


#' @export
formula.nauf.lmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                 ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::formula_nauf_merMod)
  return(eval(mc, parent.frame()))
}


#' @export
formula.nauf.glmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                  ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::formula_nauf_merMod)
  return(eval(mc, parent.frame()))
}


###### terms ######

terms_nauf_merMod <- function(x, fixed.only = TRUE, random.only = FALSE, ...) {
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
  
  attr(tt, "nauf.info") <- attr(terms(x@frame), "nauf.info")
  class(tt) <- c("nauf.terms", "terms", "formula")
  
  return(tt)
}


#' @export
terms.nauf.lmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                               ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::terms_nauf_merMod)
  return(eval(mc, parent.frame()))
}


#' @export
terms.nauf.glmerMod <- function(x, fixed.only = FALSE, random.only = FALSE,
                                ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::terms_nauf_merMod)
  return(eval(mc, parent.frame()))
}


###### model.frame ######

#' @export
model.frame.nauf.formula <- function(formula, data = NULL, subset = NULL,
                                   na.action = na.pass,
                                   drop.unused.levels = TRUE, xlev = NULL,
                                   ncs_scale = NULL, ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf::nauf_model.frame)
  return(eval(mc, parent.frame()))
}


#' @export
model.frame.nauf.terms <- function(formula, data = NULL, subset = NULL,
                                   na.action = na.pass,
                                   drop.unused.levels = TRUE, xlev = NULL,
                                   ncs_scale = NULL, ...) {
  mc <- match.call()
  info <- attr(formula, "nauf.info")
  class(formula) <- c("terms", "formula")
  ncs <- info$ncs_scale
  mc[[1]] <- quote(stats::model.frame)
  mc$formula <- formula
  mc$na.action <- na.pass
  mc$drop.unused.levels <- TRUE
  mc$xlev <- NULL
  mc["ncs_scale"] <- NULL
  mf <- eval(mc, parent.frame())
  attr(attr(mf, "terms"), "nauf.info") <- info
  class(attr(mf, "terms")) <- c("nauf.terms", "terms", "formula")
  attr(mf, "formula") <- stats::formula(attr(mf, "terms"))
  class(attr(mf, "formula")) <- c("nauf.formula", "formula")
  class(mf) <- c("data.frame", "nauf.frame")
    
  for (j in colnames(mf)) {
    if (j %in% names(info$uf)) {
      mf[[j]] <- standardize::fac_and_contr(mf[[j]], ordered = FALSE,
        levels = info$uf[[j]][[1]],
        contrasts = named_contr_sum(info$uf[[j]][[1]], ncs))
    } else if (j %in% names(info$of)) {
      mf[[j]] <- standardize::fac_and_contr(mf[[j]], ordered = TRUE,
        levels = info$of[[j]]$levels, contrasts = info$of[[j]]$contrasts)
    } else if (j %in% names(info$groups)) {
      mf[[j]] <- factor(mf[[j]], levels = info$groups[[j]])
    }
  }

  return(mf)
}


###### model.matrix ######

#' @export
model.matrix.nauf.terms <- function(object, data = environment(object),
                                    contrasts.arg = NULL, xlev = NULL, ...) {
  if (!is.nauf.frame(data)) {
    data <- model.frame(object, data)
  }
  return(nauf_mm(data))
}


###### predict ######

predict_nauf_merMod <- function(object, newdata = NULL, newparams = NULL,
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
  re.form <- lme4_reFormHack(re.form, ReForm, REForm, REform)
  if (!is.null(re.form) && !is.na(re.form) && !isTRUE(all.equal(re.form, ~ 0))) {
    stop("New random effects structures not currently supported")
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
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  if ((is.null(newdata) && is.null(re.form) && is.null(newparams))) {
    if (lme4_isLMM(object) || lme4_isNLMM(object)) {
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
      if (is.null(re.form)) {
        Terms <- attr(object@frame, "terms")
      } else {
        RHS <- stats::formula(substitute(~R, list(R = lme4_RHSForm(
          formula(object, fixed.only = TRUE)))))
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
      if (is.numeric(X.col.dropped) && length(X.col.dropped) > 0) {
        X <- X[, -X.col.dropped, drop = FALSE]
      }
    }

    pred <- drop(X %*% lme4::fixef(object))
    pred <- pred + offset
    
    if (is.null(re.form)) {
      lvs <- lapply(object@flist, levels)
      z <- Matrix::t(nauf_mkReTrms(mfnew, lvs)$Zt)
      b <- as.vector(lme4::getME(object, "b"))
      pred <- pred + base::drop(as(z %*% b, "matrix"))
    }
    
    if (lme4_isGLMM(object) && type == "response") {
      pred <- object@resp$family$linkinv(pred)
    }
  }

  return(pred)
}


#' @export
predict.nauf.lmerMod <- function(object, newdata = NULL, newparams = NULL,
                                 re.form = NULL, ReForm, REForm, REform,
                                 terms = NULL, type = c("link", "response"), 
                                 allow.new.levels = FALSE, na.action = na.pass,
                                 ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::predict_nauf_merMod)
  return(eval(mc, parent.frame()))
}


#' @export
predict.nauf.glmerMod <- function(object, newdata = NULL, newparams = NULL,
                                  re.form = NULL, ReForm, REForm, REform,
                                  terms = NULL, type = c("link", "response"), 
                                  allow.new.levels = FALSE, na.action = na.pass,
                                  ...) {
  mc <- match.call()
  mc[[1]] <- quote(nauf:::predict_nauf_merMod)
  return(eval(mc, parent.frame()))
}


###### pmmeans ######


###### lsmeans ######

