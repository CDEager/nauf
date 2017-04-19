

#' @export
nauf_lm <- function(formula, data = NULL, subset, weights, na.action = na.pass,
                    method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE,
                    singular.ok = TRUE, contrasts = NULL, offset,
                    ncs_scale = NULL, ...) {
  # based on stats::lm
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  if (!model) warning("Ignoring 'model'; must be TRUE")
  if (!x) warning("Ignoring 'x'; must be TRUE")
  if (!y) warning("Ignoring 'y'; must be TRUE")
  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  m <- match(c("formula", "data", "subset", "weights", "offset", "ncs_scale"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  mf <- eval(mf, parent.frame())

  if (method == "model.frame") return(mf)

  if (method != "qr") {
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
      method), domain = NA)
  }

  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")

  w <- as.vector(stats::model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")

  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset) && length(offset) != NROW(y)) {
    stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
      length(offset), NROW(y)), domain = NA)
  }

  if (stats::is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 3) else numeric(),
      residuals = y, fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if (!is.null(w)) sum(w != 0) else if (is.matrix(y)) nrow(y)
      else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }

  } else {
    x <- nauf_model.matrix(mf)
    if (is.null(w)) {
      z <- stats::lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...)
    } else {
      z <- stats::lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
        ...)
    }
  }

  class(z) <- c("nauf.glm", if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- NULL
  z$xlevels <- stats::.getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  z$x <- x
  z$y <- y
  if (!qr) z$qr <- NULL

  return(z)
}


#' @export
nauf_glm <- function(formula, family = gaussian, data = NULL, weights, subset,
                     na.action = na.pass, start = NULL, etastart, mustart,
                     offset, control = list(...), model = TRUE,
                     method = "glm.fit", x = TRUE, y = TRUE, contrasts = NULL,
                     ncs_scale = NULL, ...) {
  # based on stats::glm
  
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)

  if (!model) warning("Ignoring 'model'; must be TRUE")
  if (!x) warning("Ignoring 'x'; must be TRUE")
  if (!y) warning("Ignoring 'y'; must be TRUE")
  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  m <- match(c("formula", "data", "subset", "weights", "ncs_scale", "etastart",
    "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  mf <- eval(mf, parent.frame())

  if (identical(method, "model.frame")) return(mf)
  if (!is.character(method) && !is.function(method)) {
    stop("invalid 'method' argument")
  }
  if (identical(method, "glm.fit")) {
    control <- do.call(stats::glm.control, control)
  }

  mt <- attr(mf, "terms")
  Y <- stats::model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) names(Y) <- nm
  }

  if (!stats::is.empty.model(mt)) {
    X <- nauf_model.matrix(mf)
  } else {
    X <- matrix(, NROW(Y), 0L)
  }

  weights <- as.vector(stats::model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }

  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset) && length(offset) != NROW(Y)) {
    stop(gettextf("number of offsets is %d should equal %d (number of observations)",
      length(offset), NROW(Y)), domain = NA)
  }

  mustart <- stats::model.extract(mf, "mustart")
  etastart <- stats::model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
    x = X, y = Y, weights = weights, start = start, etastart = etastart,
    mustart = mustart, offset = offset, family = family, control = control,
    intercept = attr(mt, "intercept") > 0L))

  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
      x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights,
      offset = offset, family = family, control = control, intercept = TRUE))
    if (!fit2$converged) {
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    }
    fit$null.deviance <- fit2$deviance
  }

  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  fit$x <- X
  fit <- c(fit, list(call = call, formula = formula, terms = mt, data = data,
    offset = offset, control = control, method = method, contrasts = NULL,
    xlevels = .getXlevels(mt, mf)))
  class(fit) <- c("nauf.glm", fit$class, c("glm", "lm"))

  return(fit)
}


#' @export
nauf_glm.nb <- function(formula, data = NULL, weights, subset,
                        na.action = na.pass, start = NULL, etastart, mustart,
                        control = stats::glm.control(...), method = "glm.fit",
                        model = TRUE, x = TRUE, y = TRUE, contrasts = NULL,
                        ncs_scale = NULL, ..., init.theta, link = log) {
  # based on MASS::glm.nb
  
  loglik <- function(n, th, mu, y, w) {
    sum(w * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
      log(mu + (y == 0)) - (th + y) * log(th + mu)))
  }

  link <- substitute(link)
  if (missing(init.theta)) {
    fam0 <- do.call(stats::poisson, list(link = link))
  } else {
    fam0 <- do.call(MASS::negative.binomial, list(theta = init.theta,
      link = link))
  }

  if (!model) warning("Ignoring 'model'; must be TRUE")
  if (!x) warning("Ignoring 'x'; must be TRUE")
  if (!y) warning("Ignoring 'y'; must be TRUE")
  if (!is.null(contrasts)) warning("Ignoring 'contrasts'; must be NULL")
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action'; must be na.pass")
  }

  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "ncs_scale", "etastart",
    "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  mf <- eval.parent(mf)

  Terms <- attr(mf, "terms")
  if (method == "model.frame") return(mf)

  Y <- stats::model.response(mf, "numeric")
  if (!stats::is.empty.model(Terms)) {
    X <- nauf_model.matrix(mf)
  } else {
    X <- matrix(, NROW(Y), 0)
  }

  w <- stats::model.weights(mf)
  if (!length(w)) {
    w <- rep(1, nrow(mf))
  } else if (any(w < 0)) {
    stop("negative weights not allowed")
  }

  offset <- stats::model.offset(mf)
  mustart <- stats::model.extract(mf, "mustart")
  etastart <- stats::model.extract(mf, "etastart")
  n <- length(Y)

  if (!missing(method)) {
    if (!exists(method, mode = "function")) {
      stop(gettextf("unimplemented method: %s", sQuote(method)), domain = NA)
    }
    glm.fitter <- get(method)
  } else {
      method <- "glm.fit"
      glm.fitter <- stats::glm.fit
  }

  if (control$trace > 1) message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
    mustart = mustart, offset = offset, family = fam0,
    control = list(maxit = control$maxit, epsilon = control$epsilon,
    trace = control$trace > 1), intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")

  mu <- fit$fitted.values
  th <- as.vector(MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
    trace = control$trace > 2))
  if (control$trace > 1) {
    message(gettextf("Initial value for 'theta': %f", signif(th)), domain = NA)
  }
  fam <- do.call(MASS::negative.binomial, list(theta = th, link = link))

  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1

  while ((iter <- iter + 1) <= control$maxit &&
  (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, offset = offset,
      family = fam, control = list(maxit = control$maxit,
      epsilon = control$epsilon, trace = control$trace > 1),
      intercept = attr(Terms, "intercept") > 0)
    t0 <- th
    th <- MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
      trace = control$trace > 2)
    fam <- do.call(MASS::negative.binomial, list(theta = th, link = link))
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w)
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", iter, signif(th),
        signif(Dev)), domain = NA)
    }
  }

  if (!is.null(attr(th, "warn"))) fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }

  if (length(offset) && attr(Terms, "intercept")) {
    if (length(Terms)) {
      null.deviance <- glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
        offset = offset, family = fam, control = list(maxit = control$maxit,
        epsilon = control$epsilon, trace = control$trace > 1),
        intercept = TRUE)$deviance
    } else {
      null.deviance <- fit$deviance
    }
    fit$null.deviance <- null.deviance
  }

  class(fit) <- c("nauf.glm", "negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  fit$x <- X
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- NULL
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset

  return(fit)
}

