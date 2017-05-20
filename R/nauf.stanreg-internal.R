

nauf_pp_data <- function(object, newdata = NULL, re.form = NULL, offset = NULL, 
                         ...) {
  if (is.nauf.stanmer(object)) {
    out <- nauf_pp_data_mer(object, newdata = newdata, re.form = re.form, ...)
    if (!is.null(offset)) out$offset <- offset
    
    return(out)
  }
  
  return(nauf_pp_data_fer(object, newdata = newdata, offset = offset, ...))
}


nauf_pp_data_mer <- function(object, newdata, re.form, ...) {
  x <- nauf_pp_data_mer_x(object, newdata, ...)
  z <- nauf_pp_data_mer_z(object, newdata, re.form, ...)
  
  offset <- model.offset(model.frame(object))
  if (!missing(newdata) && (!is.null(offset) || !is.null(object$call$offset))) {
    offset <- try(eval(object$call$offset, newdata), silent = TRUE)
    if (!is.numeric(offset)) offset <- NULL
  }
  
  return(rsa_nlist(x, offset = offset, Zt = z$Zt, Z_names = z$Z_names))
}


nauf_pp_data_mer_x <- function(object, newdata, ...) {
  if (is.null(newdata)) return(rstanarm::get_x(object))
  return(model.matrix(stats::delete.response(terms(object)), newdata))
}


nauf_pp_data_mer_z <- function(object, newdata, re.form = NULL,
                               allow.new.levels = TRUE, na.action = na.pass) {
  if (!condition_on_re(re.form)) return(list())
  if (is.null(newdata)) return(list(Zt = object$glmod$reTrms$Zt))
  
  mt <- stats::delete.response(terms(object, fixed.only = FALSE))
  attr(mt, "nauf.info")$allow.new.levels <- allow.new.levels
  mf <- model.frame(mt, newdata)
  # add to model.frame to keep track of new levels?
  
  # not forming Z_names (should only be necessary when new levels are tracked)
  return(list(Zt = nauf_mkReTrms(mf, lapply(object$glmod$flist,
    function(f) c(levels(f), "_NEW_")))$Zt))
}


nauf_pp_data_fer <- function(object, newdata = NULL, offset = NULL, ...) {
  if (is.null(newdata)) {
    x <- rstanarm::get_x(object)
    if (is.null(offset)) {
      if (is.null(offset <- object$offset)) {
        offset <- rep(0, nrow(x))
      }
    }
    
  } else {
    x <- model.matrix(stats::delete.response(terms(object)), newdata)
    offset <- nauf_pp_data_offset(object, newdata, offset)
  }
  
  return(rsa_nlist(x, offset))
}


nauf_pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    if (is.null(offset)) {
      if (is.null(offset <- object$offset)) {
        offset <- model.offset(model.frame(object))
      }
    }
    
  } else {
    if (!is.null(offset)) {
      stopifnot(length(offset) == nrow(newdata))
    } else {
      if (!is.null(object$call$offset) || !rsa_null_or_zero(object$offset) || 
        !rsa_null_or_zero(model.offset(model.frame(object)))) {
        warning("'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", call. = FALSE)
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  
  return(offset)
}


nauf_pp_eta <- function(object, data, draws = NULL) {
  x <- data$x
  S <- rsa_posterior_sample_size(object)
  if (is.null(draws)) draws <- S
  if (draws > S) stop("'draws' should be <= posterior sample size (", S, ").")
  some_draws <- isTRUE(draws < S)
  if (some_draws) samp <- sample(S, draws)
  
  if (is.null(data$Zt)) {
    stanmat <- as.matrix(object)
    beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
    if (some_draws) beta <- beta[samp, , drop = FALSE]
    eta <- rsa_linear_predictor.matrix(beta, x, data$offset)
    
  } else {
    stanmat <- as.matrix(object$stanfit)
    beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
    if (some_draws) beta <- beta[samp, , drop = FALSE]
    eta <- rsa_linear_predictor.matrix(beta, x, data$offset)
    b <- stanmat[, grepl("^b\\[", colnames(stanmat)), drop = FALSE]
    if (some_draws) b <- b[samp, , drop = FALSE]
    eta <- eta + as.matrix(b %*% data$Zt)
  }
  
  return(rsa_nlist(eta, stanmat))
}


nauf_pp_args <- function(object, data) {
  stopifnot(is.nauf.stanreg(object), is.matrix(data$stanmat))
  
  args <- list(mu = object$family$linkinv(data$eta))
  if (length(disp <- disp_name(object$family$family))) {
    args[[disp[1]]] <- data$stanmat[, disp[2]]
  }
  
  return(args)
}


disp_name <- function(fname) {
  if (fname == "gaussian") {
    return(c("sigma", "sigma"))
  } else if (fname == "Gamma") {
    return(c("shape", "shape"))
  } else if (fname == "inverse.gaussian") {
    return(c("lambda", "lambda"))
  } else if (fname == "neg_binomial_2") {
    return(c("size", "reciprocal_dispersion"))
  }
  return(NULL)
}


nauf_ll_args <- function(object, newdata = NULL, offset = NULL,
                         reloo_or_kfold = FALSE, ...) {
  f <- family(object)
  draws <- rsa_nlist(f)
  has_newdata <- !is.null(newdata)
  
  if (has_newdata && reloo_or_kfold && !is.nauf.stanmer(object)) {
    # this should be fine as long as reloo and kfold are fixe
    dots <- list(...)
    x <- dots$newx
    stanmat <- dots$stanmat
    y <- eval(formula(object)[[2L]], newdata)
    
  } else if (has_newdata) {
    ppdat <- nauf_pp_data(object, as.data.frame(newdata), offset = offset)
    tmp <- nauf_pp_eta(object, ppdat)
    eta <- tmp$eta
    stanmat <- tmp$stanmat
    x <- ppdat$x
    y <- eval(formula(object)[[2L]], newdata)
    
  } else {
    stanmat <- as.matrix(object)
    x <- rstanarm::get_x(object)
    y <- rstanarm::get_y(object)
  }
  
  fname <- f$family
  if (fname == "binomial") {
    data <- data.frame(y, x)
  } else {
    trials <- 1
    if (is.factor(y)) y <- rsa_fac2bin(y)
    stopifnot(all(y %in% c(0, 1)))
    data <- data.frame(y, trials, x)
  }
  
  draws$beta <- stanmat[, seq_len(ncol(x)), drop = FALSE]
  if (length(disp <- disp_name(fname))) {
    draws[[disp[1]]] <- stanmat[, disp[2]]
  }  
  
  if (has_newdata) {
    data$offset <- offset
  } else {
    data$offset <- object$offset
  }
  if (rsa_model_has_weights(object)) data$weights <- object$weights
  
  if (is.nauf.stanmer(object)) {
    b <- stanmat[, b_names(colnames(stanmat)), drop = FALSE]
    if (has_newdata) {
      if (is.null(ppdat$Zt)) {
        z <- matrix(NA, nrow = nrow(x), ncol = 0)
      } else {
        z <- t(ppdat$Zt)
      }
    } else {
      z <- rstanarm::get_z(object)
    }
    data <- cbind(data, as.matrix(z))
    draws$beta <- cbind(draws$beta, b)
  }
  
  return(rsa_nlist(data, draws, S = NROW(draws$beta), N = nrow(data)))
}

