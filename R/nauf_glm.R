## these functions are adapted from stats::lm, stats::glm, and MASS::glm.nb


nauf_lm <- function(formula, data, subset, weights, na.action = "na.pass",
                    method = "qr", model = TRUE, x = FALSE, y = FALSE,
                    qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset,
                    sumcoef = 1, ...) {
  cl <- match.call()
  
  if (!model) warning("forcing model = TRUE")
  model <- TRUE
  
  if (!is.null(contrasts)) warning("forcing contrasts = NULL")
  contrasts <- NULL
  
  if (na.action != "na.pass") warning("forcing na.action = 'na.pass'")
  na.action <- "na.pass"
  
  ret.x <- x
  ret.y <- y
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
    "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf$sumcoef <- sumcoef
  mf[[1L]] <- quote(nauf::nauf_model.frame)
  mf <- eval(mf, parent.frame())
  
  if (method == "model.frame") {
    return(mf)
  }
  else if (method != "qr") {
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
      method), domain = NA)
  }
  
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  w <- as.vector(stats::model.weights(mf))
  
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) {
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
        length(offset), NROW(y)), domain = NA)
    }
  }
  
  if (stats::is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 3) else numeric(),
      residuals = y, fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if (!is.null(w)) sum(w != 0)
                    else if (is.matrix(y)) nrow(y)
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
      z <- stats::lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, ...)
    }
  }
  
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- "na.pass"
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  if (ret.x) z$x <- x
  if (ret.y) z$y <- y
  if (!qr) z$qr <- NULL
  
  return(nauf.glm(z))
}


nauf_glm <- function(formula, data, family = gaussian, weights, subset, 
                     na.action = "na.pass", start = NULL, etastart, mustart,
                     offset, control = list(...), model = TRUE,
                     method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
                     sumcoef = 1, ...) {
  call <- match.call()
  family <- get_family(family)
  if (isTRUE(all.equal(family, gaussian()))) {
    warning("calling nauf_lm")
    
    
  } else if (is.character(family)) {
    warning("calling nauf_glm.nb")
    
    
  }
  
  if (!model) warning("forcing model = TRUE")
  model <- TRUE
  
  if (!is.null(contrasts)) warning("forcing contrasts = NULL")
  contrasts <- NULL
  
  if (na.action != "na.pass") warning("forcing na.action = 'na.pass'")
  na.action <- "na.pass"
  
  if (missing(data)) stop("must supply a data.frame")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
      "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$sumcoef <- sumcoef
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
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) names(Y) <- nm
  }
  if (!is.empty.model(mt)) {
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
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
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
  fit$na.action <- "na.pass"
  if (x) fit$x <- X
  if (!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt, 
    data = data, offset = offset, control = control, method = method, 
    contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  
  return(nauf.glm(fit))
}


nauf_glm.nb <- function(formula, data, weights, subset, na.action = "na.pass",
                        start = NULL, etastart, mustart,
                        control = stats::glm.control(...), method = "glm.fit", 
                        model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,
                        sumcoef = 1, ..., init.theta, link = log) {

  if (!model) warning("forcing model = TRUE")
  model <- TRUE
  
  if (!is.null(contrasts)) warning("forcing contrasts = NULL")
  contrasts <- NULL
  
  if (na.action != "na.pass") warning("forcing na.action = 'na.pass'")
  na.action <- "na.pass"

  loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + y) - lgamma(th) -
    lgamma(y + 1) + th * log(th) + y * log(mu + (y == 0)) - (th + y) *
    log(th + mu)))
    
  link <- substitute(link)
  if (missing(init.theta)) {
    fam0 <- do.call(stats::poisson, list(link = link))
  } else {
    fam0 <- do.call(MASS::negative.binomial, list(theta = init.theta,
      link = link))
  }
  
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
      "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf$sumcoef <- sumcoef
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
      stop(gettextf("unimplemented method: %s", sQuote(method)), 
        domain = NA)
    }
    glm.fitter <- get(method)
  } else {
      method <- "glm.fit"
      glm.fitter <- stats::glm.fit
  }
  
  if (control$trace > 1) message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart, 
    mustart = mustart, offset = offset, family = fam0, control = list(
    maxit = control$maxit, epsilon = control$epsilon,
    trace = control$trace > 1), intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
      trace = control$trace > 2))
  if (control$trace > 1) {
    message(gettextf("Initial value for 'theta': %f", signif(th)), 
      domain = NA)
  }
  
  fam <- do.call(MASS::negative.binomial, list(theta = th, link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  
  while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - Lm)/d1 +
  abs(del)/d2) > control$epsilon) {
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
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
        iter, signif(th), signif(Dev)), domain = NA)
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
  
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  fit$model <- mf
  fit$na.action <- "na.pass"
  if (x) fit$x <- X
  if (!y) fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  
  return(nauf.glm(fit))
}

