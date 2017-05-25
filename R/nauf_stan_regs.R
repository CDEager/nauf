

#' Fit a Bayesian mixed effects regression with \code{nauf} contrasts.
#'
#' ADD DESCRIPTION
#'
#' ADD DETAILS
#'
#' @param formula,data,family,subset,weights,na.action,offset,contrasts,ncs_scale See
#'   \code{\link{nauf_model.frame}} and \code{\link{nauf_glFormula}}.
#' @param ... Further arguments to be passed to \code{\link[rstan]{sampling}}.
#'   See \code{\link[rstanarm]{stan_glmer}} for details.
#' @param link,prior,prior_intercept,prior_aux,prior_covariance,prior_PD,adapt_delta,QR,sparse See
#'   \code{\link[rstanarm]{stan_glmer}}.
#' @param algorithm Changes from the default \code{"sampling"} result in an
#'   error.  Only MCMC is currently supported.
#'
#' @return A \code{\link{nauf.stanreg}} object.
#'
#' @examples
#' \dontrun{
#' dat <- fricatives
#' # add NA stuff
#' mod <- nauf_stan_lmer(f, d)
#' }
#' 
#' @seealso \code{\link{nauf_contrasts}} for a description of the treatment of 
#'   \code{NA} values, \code{\link[rstanarm]{stan_glmer}} for a description of 
#'   the priors, and the documentation for \code{Stan} and the \code{rstan} and 
#'   \code{rstanarm} packages for algorithmic details.
#'
#' @export
nauf_stan_glmer <- function(formula, data = NULL, family = gaussian, subset,
                            weights, na.action = na.pass, offset, contrasts = NULL,
                            ncs_scale = attr(formula, "standardized.scale"),
                            ..., prior = rstanarm::normal(),
                            prior_intercept = rstanarm::normal(),
                            prior_aux = rstanarm::cauchy(0, 5),
                            prior_covariance = rstanarm::decov(), prior_PD = FALSE,
                            algorithm = "sampling", adapt_delta = NULL, QR = FALSE,
                            sparse = FALSE) {
  # based on rstanarm::stan_glmer
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  
  if (algorithm != "sampling") {
    stop("Only algorithm = \"sampling\" is supported.")
  }
  
  family <- get_family(family)
  
  mc[[1]] <- quote(nauf::nauf_glFormula)
  mc$control <- lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop", check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  mc <- mc[c(1, which(names(mc) %in% setdiff(intersect(names(mc),
    names(formals(nauf::nauf_glFormula))), "...")))]
  glmod <- eval(mc, parent.frame())
  
  X <- glmod$X
  y <- glmod$fr[, as.character(glmod$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L) y <- as.vector(y)
  
  offset <- model.offset(glmod$fr)
  if (is.null(offset)) offset <- double(0)
  if (missing(weights) || is.null(weights)) {
    weights <- double(0)
  } else {
    if (!is.numeric(weights)) {
      stop("'weights' must be a numeric vector.", call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("Negative weights are not allowed.", call. = FALSE)
    }
    if (length(weights) != nrow(glmod$fr)) {
      stop("'weights' must have exactly one element per row of 'data'",
        call. = FALSE)
    }
  }

  if (is.null(prior)) prior <- list()
  if (is.null(prior_intercept)) prior_intercept <- list()
  if (is.null(prior_aux)) prior_aux <- list()
  if (is.null(prior_covariance)) {
    stop("'prior_covariance' can't be NULL.", call. = FALSE)
  }
  
  group <- glmod$reTrms
  group$decov <- prior_covariance
  
  stanfit <- rstanarm::stan_glm.fit(x = X, y = y, weights = weights,
    offset = offset, family = family, prior = prior,
    prior_intercept = prior_intercept, prior_aux = prior_aux,
    prior_PD = prior_PD, algorithm = algorithm, adapt_delta = adapt_delta,
    group = group, QR = QR, sparse = sparse, ...)
    
  Z <- rsa_pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms,
    flist = group$flist)$Z
  colnames(Z) <- rsa_b_names(names(stanfit), value = TRUE)
  
  fit <- rsa_nlist(stanfit, family, formula, offset, weights,
    x = if (getRversion() < "3.2.0") Matrix::cBind(X, Z) else cbind2(X, Z),
    y = y, data, call, terms = NULL, model = NULL, na.action = na.pass,
    contrasts, algorithm, glmod)
  
  out <- rsa_stanreg(fit)
  add_class(out) <- list("nauf.stanreg", "lmerMod")
  out$modeling_function <- "nauf_stan_glmer"
  
  return(out)
}


#' @rdname nauf_stan_glmer
#'
#' @export
nauf_stan_lmer <- function(formula, data = NULL, subset, weights,
                           na.action = na.pass, offset, contrasts = NULL,
                           ncs_scale = attr(formula, "standardized.scale"),
                           ..., prior = rstanarm::normal(),
                           prior_intercept = rstanarm::normal(),
                           prior_aux = rstanarm::cauchy(0, 5),
                           prior_covariance = rstanarm::decov(), prior_PD = FALSE,
                           algorithm = "sampling", adapt_delta = NULL,
                           QR = FALSE) {
  # based on rstanarm::stan_lmer
  if ("family" %in% names(list(...))) stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call)) names(call)[2L] <- "formula"
  mc[[1L]] <- quote(nauf::nauf_stan_glmer)
  mc$REML <- NULL
  mc$family <- "gaussian"
  
  out <- eval(mc, parent.frame())
  out$call <- call
  out$modeling_function <- "nauf_stan_lmer"
  
  return(out)
}


#' @rdname nauf_stan_glmer
#'
#' @export
nauf_stan_glmer.nb <- function(formula, data = NULL, subset, weights,
                               na.action = na.pass, offset, contrasts = NULL,
                               link = "log",
                               ncs_scale = attr(formula, "standardized.scale"),
                               ..., prior = rstanarm::normal(),
                               prior_intercept = rstanarm::normal(),
                               prior_aux = rstanarm::cauchy(0, 5),
                               prior_covariance = rstanarm::decov(),
                               prior_PD = FALSE, algorithm = "sampling",
                               adapt_delta = NULL, QR = FALSE) {
  # based on rstanarm::stan_glmer.nb
  if ("family" %in% names(list(...))) stop("'family' should not be specified.")
  mc <- call <- match.call(expand.dots = TRUE)
  if (!"formula" %in% names(call)) names(call)[2L] <- "formula"
  mc[[1]] <- quote(nauf::nauf_stan_glmer)
  mc$REML <- mc$link <- NULL
  mc$family <- rstanarm::neg_binomial_2(link = link)
  
  out <- eval(mc, parent.frame())
  out$call <- call
  out$modeling_function <- "nauf_stan_glmer.nb"
  
  return(out)
}


#' Fit a Bayesian fixed effects regression with \code{nauf} contrasts.
#'
#' ADD DESCRIPTION
#'
#' ADD DETAILS
#'
#' @param formula,data,family,subset,weights,na.action,offset,contrasts,model,x,y,ncs_scale See
#'   \code{\link{nauf_model.frame}}.
#' @param ... Further arguments to be passed to \code{\link[rstan]{sampling}}.
#'   See \code{\link[rstanarm]{stan_glmer}} for details.
#' @param link,prior,prior_intercept,prior_aux,prior_PD,adapt_delta,QR,sparse See
#'   \code{\link[rstanarm]{stan_glm}}.
#' @param algorithm Changes from the default \code{"sampling"} result in an
#'   error.  Only MCMC is currently supported.
#'
#' @return A \code{\link{nauf.stanreg}} object.
#'
#' @examples
#' \dontrun{
#' dat <- fricatives
#' # add NA stuff
#' mod <- nauf_stan_lm(f, d)
#' }
#' 
#' @seealso \code{\link{nauf_contrasts}} for a description of the treatment of 
#'   \code{NA} values, \code{\link[rstanarm]{stan_glm}} for a description of 
#'   the priors, and the documentation for \code{Stan} and the \code{rstan} and 
#'   \code{rstanarm} packages for algorithmic details.
#'
#' @export
nauf_stan_glm <- function(formula, family = gaussian(), data = NULL, weights, 
                          subset, na.action = na.pass, offset = NULL,
                          model = TRUE, x = TRUE, y = TRUE, contrasts = NULL,
                          ncs_scale = attr(formula, "standardized.scale"), ...,
                          prior = rstanarm::normal(),
                          prior_intercept = rstanarm::normal(), 
                          prior_aux = rstanarm::cauchy(0, 5), prior_PD = FALSE,
                          algorithm = "sampling", adapt_delta = NULL, 
                          QR = FALSE, sparse = FALSE) {
  # based on stan_glm
  if (algorithm != "sampling") {
    stop("Only algorithm = \"sampling\" is supported.")
  }
  if (length(lme4::findbars(formula))) {
    stop("'formula' cannot contain '|'.  Mixed effects models should be fit ",
      "with nauf_stan_(g)lmer(.nb)")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  
  if (!model) warning("Ignoring 'model'; must be TRUE")
  if (!x) warning("Ignoring 'x'; must be TRUE")
  if (!y) warning("Ignoring 'y'; must be TRUE")
  
  family <- get_family(family)
  
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
      "offset", "ncs_scale"), table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  mf <- eval(mf, parent.frame())
  
  rsa_check_constant_vars(mf)
  
  mt <- attr(mf, "terms")
  Y <- rsa_array1D_check(model.response(mf, type = "any"))
  X <- nauf_mm(mf)
  
  offset <- as.vector(stats::model.offset(mf))
  if (is.null(offset)) {
    offset <- double(0)
  } else if (length(o) != nrow(mf)) {
    stop("'offset' must have exactly one element per row of 'data'")
  }
  
  weights <- as.vector(stats::model.weights(mf))
  if (is.null(weights)) {
    weights <- double(0)
  } else {
    if (!is.numeric(weights)) {
      stop("'weights' must be a numeric vector.", call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("Negative weights are not allowed.", call. = FALSE)
    }
    if (length(weights) != nrow(mf)) {
      stop("'weights' must have exactly one element per row of 'data'",
        call. = FALSE)
    }
  }
  
  if (rsa_binom_y_prop(Y, family, weights)) {
    y1 <- as.integer(as.vector(Y) * weights)
    Y <- cbind(y1, y0 = weights - y1)
    weights <- double(0)
  }
  
  stanfit <- rstanarm::stan_glm.fit(x = X, y = Y, weights = weights,
    offset = offset, family = family, prior = prior,
    prior_intercept = prior_intercept, prior_aux = prior_aux,
    prior_PD = prior_PD, algorithm = algorithm, adapt_delta = adapt_delta,
    QR = QR, sparse = sparse, ...)
      
  fit <- rsa_nlist(stanfit, algorithm, family, formula, data, offset, 
    weights, x = X, y = Y, model = mf, terms = mt, call,
    na.action = na.pass, contrasts = NULL,
    modeling_function = "nauf_stan_glm")
    
  out <- rsa_stanreg(fit)
  out$xlevels <- .getXlevels(mt, mf)
  first_class(out) <- "nauf.stanreg"
  
  return(out)
}


#' @rdname nauf_stan_glm
#'
#' @export
nauf_stan_lm <- function(formula, data = NULL, subset, weights,
                         na.action = na.pass, model = TRUE, x = TRUE, y = TRUE,
                         singular.ok = TRUE, contrasts = NULL, offset,
                         ncs_scale = attr(formula, "standardized.scale"), ...,
                         prior = rstanarm::R2(stop("'location' must be specified")), 
                         prior_intercept = NULL, prior_PD = FALSE,
                         algorithm = "sampling", adapt_delta = NULL) {
  # based on stan_lm
  if (algorithm != "sampling") {
    stop("Only algorithm = \"sampling\" is supported.")
  }
  if (length(lme4::findbars(formula))) {
    stop("'formula' cannot contain '|'.  Mixed effects models should be fit ",
      "with nauf_stan_(g)lmer(.nb)")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  
  if (!model) warning("Ignoring 'model'; must be TRUE")
  if (!x) warning("Ignoring 'x'; must be TRUE")
  if (!y) warning("Ignoring 'y'; must be TRUE")
  
  call <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
      "offset", "ncs_scale"), table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  modelframe <- eval(mf, parent.frame())
  
  mt <- attr(modelframe, "terms")
  Y <- stats::model.response(modelframe, "numeric")
  X <- nauf_mm(modelframe)
  w <- as.vector(stats::model.weights(modelframe))
  offset <- as.vector(stats::model.offset(modelframe))
  
  stanfit <- rstanarm::stan_lm.wfit(y = Y, x = X, w, offset,
    singular.ok = singular.ok, prior = prior, prior_intercept = prior_intercept,
    prior_PD = prior_PD, algorithm = algorithm, adapt_delta = adapt_delta, ...)

  fit <- rsa_nlist(stanfit, family = gaussian(), formula, offset, 
    weights = w, x = X[, intersect(colnames(X), dimnames(stanfit)[[3]]),
    drop = FALSE], y = Y, data, prior.info = prior, algorithm, call, terms = mt, 
    model = modelframe, na.action = na.pass, contrasts = NULL,
    modeling_function = "nauf_stan_lm")
    
  out <- rsa_stanreg(fit)
  out$xlevels <- .getXlevels(mt, modelframe)
  first_class(out) <- "nauf.stanreg"
  
  return(out)
}


#' @rdname nauf_stan_glm
#'
#' @export
nauf_stan_glm.nb <- function(formula, data = NULL, weights, subset,
                             na.action = na.pass, offset = NULL, model = TRUE,
                             x = TRUE, y = TRUE, contrasts = NULL, link = "log", 
                             ncs_scale = attr(formula, "standardized.scale"),
                             ..., prior = rstanarm::normal(),
                             prior_intercept = rstanarm::normal(),
                             prior_aux = rstanarm::cauchy(0, 5), prior_PD = FALSE,
                             algorithm = "sampling", adapt_delta = NULL,
                             QR = FALSE) {
  # based on stan_glm.nb
  if ("family" %in% names(list(...))) stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call)) names(call)[2L] <- "formula"
  mc[[1L]] <- quote(nauf::nauf_stan_glm)
  mc$link <- NULL
  mc$family <- rstanarm::neg_binomial_2(link = link)
  out <- eval(mc, parent.frame())
  out$call <- call
  out$modeling_function <- "nauf_stan_glm.nb"
  return(out)
}

