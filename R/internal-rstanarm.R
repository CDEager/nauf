

## this file contains internal functions taken from rstanarm v2.15.3
## which are altered as little as possible and preceded by 'rsa_'


#' @importFrom Matrix rBind Matrix
rsa_pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    if (grepl("^Xr", names(p)[i])) next
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])),
      paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    if (grepl("^Xr", names(p)[i])) next
    if (getRversion() < "3.2.0") {
      Ztlist[[i]] <- Matrix::rBind(Ztlist[[i]], Matrix::Matrix(0, nrow = p[i],
        ncol = n, sparse = TRUE))
    } else {
      Ztlist[[i]] <- rbind2(Ztlist[[i]], Matrix::Matrix(0, nrow = p[i],
        ncol = n, sparse = TRUE))
    }
  }
  
  Z <- Matrix::t(do.call(rbind, args = Ztlist))
  
  return(rsa_nlist(Z, cnms, flist))
}


rsa_b_names <- function(x, ...) {
  return(grep("^b\\[", x, ...))
}


rsa_nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  
  if (no_names <- is.null(names(out))) {
    has_name <- FALSE
  } else {
    has_name <- nzchar(names(out))
  }
  
  if (all(has_name)) return(out)
  
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }
  
  return(out)
}


#' @importFrom utils packageVersion
rsa_stanreg <- function(object) {
  opt <- object$algorithm == "optimizing"
  mer <- !is.null(object$glmod)
  stanfit <- object$stanfit
  family <- object$family
  y <- object$y
  x <- object$x
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  if (is.matrix(y)) {
    ynames <- rownames(y)
  } else {
    ynames <- names(y)
  }
  
  if (is_betareg <- family$family == "beta") {
    family_phi <- object$family_phi
    z <- object$z
    nvars_z <- ncol(z)
  }
  
  if (opt) {
    stanmat <- stanfit$theta_tilde
    probs <- c(0.025, 0.975)
    stan_summary <- cbind(Median = apply(stanmat, 2L, median), 
      MAD_SD = apply(stanmat, 2L, mad), t(apply(stanmat, 2L, quantile, probs)))
    xnms <- colnames(x)
    covmat <- cov(stanmat)[xnms, xnms]
    coefs <- apply(stanmat[, xnms, drop = FALSE], 2L, median)
    ses <- apply(stanmat[, xnms, drop = FALSE], 2L, mad)
    rank <- qr(x, tol = .Machine$double.eps, LAPACK = TRUE)$rank
    df.residual <- nobs - sum(object$weights == 0) - rank
    if (is_betareg) {
      if (length(colnames(z)) == 1) {
        coefs_z <- apply(stanmat[, grepl("(phi)", colnames(stanmat),
          fixed = TRUE), drop = FALSE], 2L, median)
      } else {
        coefs_z <- apply(stanmat[, paste0("(phi)_", colnames(z)), drop = FALSE],
          2L, median)
      }
    }

  } else {
    stan_summary <- rsa_make_stan_summary(stanfit)
    coefs <- stan_summary[1:nvars, rsa_select_median(object$algorithm)]
      
    if (is_betareg) {
      coefs_z <- stan_summary[(nvars + 1):(nvars + nvars_z), 
        rsa_select_median(object$algorithm)]
      if (length(coefs_z) == 1L) {
        names(coefs_z) <- rownames(stan_summary)[nvars + 1]
      }
    }
      
    if (length(coefs) == 1L) names(coefs) <- rownames(stan_summary)[1L]
    if (is_betareg) {
      stanmat <- as.matrix(stanfit)[, c(names(coefs), names(coefs_z)),
        drop = FALSE]
      colnames(stanmat) <- c(names(coefs), names(coefs_z))
    } else {
      stanmat <- as.matrix(stanfit)[, 1:nvars, drop = FALSE]
      colnames(stanmat) <- colnames(x)
    }
      
    ses <- apply(stanmat, 2L, mad)
    if (mer) {
      mark <- sum(sapply(object$stanfit@par_dims[c("alpha", "beta")], prod))
      stanmat <- stanmat[, 1:mark, drop = FALSE]
    }
    covmat <- cov(stanmat)
    if (object$algorithm == "sampling") {
      rsa_check_rhats(stan_summary[, "Rhat"])
    }
  }
  
  eta <- rsa_linear_predictor.default(coefs, x, object$offset)
  mu <- family$linkinv(eta)
  if (NCOL(y) == 2L) {
    residuals <- y[, 1L] / rowSums(y) - mu
  } else {
    if (is.factor(y)) {
      ytmp <- rsa_fac2bin(y)
    } else {
      ytmp <- y
    }
    residuals <- ytmp - mu
  }
  names(eta) <- names(mu) <- names(residuals) <- ynames
  
  if (is_betareg) {
    eta_z <- rsa_linear_predictor.default(coefs_z, z, object$offset)
    phi <- family_phi$linkinv(eta_z)
  }
  
  out <- rsa_nlist(coefficients = rsa_unpad_reTrms.default(coefs),
    ses = rsa_unpad_reTrms.default(ses), 
    fitted.values = mu, linear.predictors = eta, residuals, 
    df.residual = if (opt) df.residual else NA_integer_, covmat, y, x,
    model = object$model, data = object$data, family,
    offset = if (any(object$offset != 0)) object$offset else NULL,
    weights = object$weights, prior.weights = object$weights, 
    contrasts = object$contrasts, na.action = object$na.action, 
    formula = object$formula, terms = object$terms,
    prior.info = attr(stanfit, "prior.info"), algorithm = object$algorithm,
    stan_summary, stanfit = if (opt) stanfit$stanfit else stanfit,
    rstan_version = utils::packageVersion("rstan"), 
    call = object$call, modeling_function = object$modeling_function)
      
  if (opt) out$asymptotic_sampling_dist <- stanmat
  if (mer) out$glmod <- object$glmod
  if (is_betareg) {
    out$coefficients <- rsa_unpad_reTrms.default(c(coefs, coefs_z))
    out$z <- z
    out$family_phi <- family_phi
    out$eta_z <- eta_z
    out$phi <- phi
  }
  
  return(structure(out, class = c("stanreg", "glm", "lm")))
}


#' @importFrom rstan summary
rsa_make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs)/2
  probs <- sort(c(0.5, qq, 1 - qq))
  return(rstan::summary(stanfit, probs = probs, digits = 10)$summary)
}


rsa_select_median <- function(algorithm) {
  switch(algorithm, sampling = "50%", meanfield = "50%", fullrank = "50%", 
    optimizing = "Median",
    stop("Bug found (incorrect algorithm name passed to select_median)", 
    call. = FALSE))
}


rsa_check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp) rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]
  if (any(rhats > threshold, na.rm = TRUE)) {
    warning("Markov chains did not converge! Do not analyze results!", 
      call. = FALSE, noBreaks. = TRUE)
  }
}


rsa_linear_predictor.default <- function(beta, x, offset = NULL) {
  if (is.matrix(beta)) return(rsa_linear_predictor.matrix(beta, x, offset))
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (length(offset)) eta <- eta + offset
  return(eta)
}


rsa_linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset)) eta <- sweep(eta, 2L, offset, `+`)
  return(eta)
}


rsa_fac2bin <- function(y) {
  if (!is.factor(y)) {
    stop("Bug found: non-factor as input to fac2bin.", call. = FALSE)
  }
  if (!identical(nlevels(y), 2L)) {
    stop("Bug found: factor with nlevels != 2 as input to fac2bin.",
      call. = FALSE)
  }
  return(as.integer(y != levels(y)[1L]))
}


rsa_unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x)) return(rsa_unpad_reTrms.array(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  return(x[keep])
}


rsa_unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  
  if (ndim > 3) stop("'x' should be a matrix or 3-D array")
  if (columns) {
    nms <- dimnames(x)[[ndim]]
  } else {
    nms <- rownames(x)
  }

  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (ndim == 2) {
    if (columns) {
      x_keep <- x[, keep, drop = FALSE]
    } else {
      x_keep <- x[keep, , drop = FALSE]
    }
  } else {
    if (columns) {
      x_keep <- x[, , keep, drop = FALSE]
    } else {
      x_keep <- x[keep, , , drop = FALSE]
    }
  }
  
  return(x_keep)
}


rsa_check_constant_vars <- function(mf) {
  # this one is altered slightly to not treat NA as a value
  is.const <- function(x) length(sort(unique(x))) < 2
  
  nocheck <- c(if (NCOL(mf[, 1]) == 2) colnames(mf)[1], "(weights)",
    "(offset)", "(Intercept)")
  sel <- !colnames(mf) %in% nocheck
  is_constant <- apply(mf[, sel, drop = FALSE], 2, is.const)

  if (any(is_constant)) {
    stop("Constant variable(s) found: ", paste(names(is_constant)[is_constant],
      collapse = ", "), call. = FALSE)
  }
  
  invisible(NULL)
}


rsa_array1D_check <- function(y) {
  if (length(dim(y)) == 1L) {
    nms <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nms)) names(y) <- nms
  }
  return(y)
}


rsa_binom_y_prop <- function(y, family, weights) {
  if (family$family != "binomial") return(FALSE)
  yprop <- NCOL(y) == 1L && is.numeric(y) && any(y > 0 & y < 1) &&
    !any(y < 0 | y > 1)
  if (!yprop) return(FALSE)
  wtrials <- !identical(weights, double(0)) && all(weights > 0) &&
    all(abs(weights - round(weights)) < .Machine$double.eps^0.5)
  return(isTRUE(wtrials))
}


rsa_model_has_weights <- function(x) {
  wts <- x[["weights"]]
  return(length(wts) && !all(wts == wts[1]))
}


rsa_null_or_zero <- function(x) {
  return(isTRUE(is.null(x) || all(x == 0)))
}


rsa_posterior_sample_size <- function(x) {
  return(sum(x$stanfit@sim$n_save - x$stanfit@sim$warmup2))
}




rsa_pp_fun <- function(object) {
  return(get(paste0("rsa_pp_", object$family$family), mode = "function"))
}


rsa_pp_gaussian <- function(mu, sigma) {
  t(sapply(1:nrow(mu), function(s) {
    rnorm(ncol(mu), mu[s, ], sigma[s])
  }))
}


rsa_pp_binomial <- function(mu, trials) {
  t(sapply(1:nrow(mu), function(s) {
    rbinom(ncol(mu), size = trials, prob = mu[s, ])
  }))
}


rsa_pp_poisson <- function(mu) {
  t(sapply(1:nrow(mu), function(s) {
    rpois(ncol(mu), mu[s, ])
  }))
}


rsa_pp_neg_binomial_2 <- function(mu, size) {
  t(sapply(1:nrow(mu), function(s) {
    rnbinom(ncol(mu), size = size[s], mu = mu[s, ])
  }))
}


rsa_pp_inverse.gaussian <- function(mu, lambda) {
  t(sapply(1:nrow(mu), function(s) {
    rsa_rinvGauss(ncol(mu), mu = mu[s, ], lambda = lambda[s])
  }))
}


rsa_pp_Gamma <- function(mu, shape) {
  t(sapply(1:nrow(mu), function(s) {
    rgamma(ncol(mu), shape = shape[s], rate = shape[s]/mu[s, ])
  }))
}


rsa_rinvGauss <- function(n, mu, lambda) {
  mu2 <- mu^2
  y <- rnorm(n)^2
  z <- runif(n)
  tmp <- (mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * y^2))
  x <- mu + tmp/(2 * lambda)
  ifelse(z <= (mu/(mu + x)), x, mu2/x)
}




rsa_ll_fun <- function(object) {
  return(get(paste0("rsa_ll_", object$family$family, "_i"), mode = "function"))
}


rsa_ll_gaussian_i <- function(i, data, draws) {
  val <- dnorm(data$y, mean = rsa_mu(data, draws), sd = draws$sigma, log = TRUE)
  rsa_weignted(val, data$weights)
}


rsa_ll_binomial_i <- function(i, data, draws) {
  val <- dbinom(data$y, size = data$trials, prob = rsa_mu(data, draws), log = TRUE)
  rsa_weignted(val, data$weights)
}


rsa_ll_poisson_i <- function(i, data, draws) {
  val <- dpois(data$y, lambda = rsa_mu(data, draws), log = TRUE)
  rsa_weignted(val, data$weights)
}


rsa_ll_neg_binomial_2_i <- function(i, data, draws) {
  val <- dnbinom(data$y, size = draws$size, mu = rsa_mu(data, draws), log = TRUE)
  rsa_weignted(val, data$weights)
}


rsa_ll_inverse.gaussian_i <- function(i, data, draws) {
  mu <- rsa_mu(data, draws)
  val <- 0.5 * log(draws$lambda/(2 * pi)) - 1.5 * log(data$y) - 
    0.5 * draws$lambda * (data$y - mu)^2/(data$y * mu^2)
  rsa_weignted(val, data$weights)
}


rsa_ll_Gamma_i <- function(i, data, draws) {
  val <- dgamma(data$y, shape = draws$shape, rate = draws$shape/rsa_mu(data, 
    draws), log = TRUE)
  rsa_weignted(val, data$weights)
}


rsa_mu <- function(data, draws) {
  if (is.matrix(draws$beta)) {
    eta <- as.vector(rsa_linear_predictor.matrix(draws$beta, rsa_xdata(data),
      data$offset))
  } else {
    eta <- as.vector(rsa_linear_predictor.default(draws$beta, rsa_xdata(data),
      data$offset))
  }
  return(draws$f$linkinv(eta))
}


rsa_xdata <- function(data) {
  sel <- c("y", "weights", "offset", "trials")
  return(data[, -which(colnames(data) %in% sel)])
}


rsa_weighted <- function(val, w) {
  if (is.null(w)) {
    return(val)
  }
  return(val * w)
}

