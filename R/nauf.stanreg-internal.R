

## functions in this file are based on the corresponding functions in
## rstanarm v2.15.3


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
  
  if (is.null(newdata)) {
    mf <- object$glmod$fr
  } else {
    mt <- stats::delete.response(terms(object, fixed.only = FALSE))
    nauf.info(mt, "allow.new.levels") <- allow.new.levels
    mf <- model.frame(mt, newdata)
  }
  
  # not forming Z_names (should only be necessary when new levels are tracked)
  return(list(Zt = nauf_mkReTrms(mf, levasgn(object, TRUE))$Zt))
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
  S <- rsa_posterior_sample_size(object)
  if (is.null(draws)) draws <- S
  if (draws > S) stop("'draws' should be <= posterior sample size (", S, ").")
  if (some_draws <- isTRUE(draws < S)) samp <- sample(S, draws)
  
  if (is.null(data$Zt)) {
    stanmat <- as.matrix(object)
    beta <- stanmat[, seq_len(ncol(data$x)), drop = FALSE]
    if (some_draws) beta <- beta[samp, , drop = FALSE]
    eta <- rsa_linear_predictor.matrix(beta, data$x, data$offset)
    
  } else {
    stanmat <- as.matrix(object$stanfit)
    beta <- stanmat[, seq_len(ncol(data$x)), drop = FALSE]
    if (some_draws) beta <- beta[samp, , drop = FALSE]
    eta <- rsa_linear_predictor.matrix(beta, data$x, data$offset)
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
  if (fname != "binomial") {
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
    b <- stanmat[, rsa_b_names(colnames(stanmat)), drop = FALSE]
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


#' @importFrom loo psislw
nauf_loo_weights <- function(object, lw, log = FALSE, ...) {
  if (!missing(lw)) {
    stopifnot(is.matrix(lw))
  } else {
    message("Running PSIS to compute weights...")
    psis <- loo::psislw(llfun = rsa_ll_fun(object),
      llargs = nauf_ll_args(object), ...)
    lw <- psis[["lw_smooth"]]
  }
  if (log) return(lw)
  return(exp(lw))
}


#' @importFrom bayesplot available_ppc
rsa_ppc_function_name <- function(fun = character()) {
  if (!length(fun))  stop("Plotting function not specified.", call. = FALSE)
  if (identical(substr(fun, 1, 5), "mcmc_")) {
    stop("For 'mcmc_' functions use the 'plot' ", "method instead of 'pp_check'.", 
      call. = FALSE)
  }
  if (!identical(substr(fun, 1, 4), "ppc_")) fun <- paste0("ppc_", fun)
  if (!fun %in% bayesplot::available_ppc()) {
    stop(fun, " is not a valid PPC function name.",
      " Use bayesplot::available_ppc() for a list of available PPC functions.")
  }
  return(fun)
}


rsa_ppc_y_and_yrep <- function(object, nreps = NULL, seed = NULL,
                               binned_resid_plot = FALSE) {
  y <- rstanarm::get_y(object)
  
  if (binned_resid_plot) {
    yrep <- posterior_linpred(object, transform = TRUE)
    yrep <- yrep[1:nreps, , drop = FALSE]
  } else {
    yrep <- posterior_predict(object, draws = nreps, seed = seed)
  }
  
  if (object$family$family == "binomial") {
    if (NCOL(y) == 2L) {
      trials <- rowSums(y)
      y <- y[, 1L]/trials
      if (!binned_resid_plot) yrep <- sweep(yrep, 2L, trials, "/")
    } else if (is.factor(y)) {
      y <- rsa_fac2bin(y)
    }
  }
  
  return(rsa_nlist(y, yrep))
}


rsa_set_nreps <- function(nreps = NULL, fun = character()) {
  fun <- sub("ppc_", "", fun)
  
  switch(fun,
    dens_overlay = nreps %ORifNULL% 50,
    ecdf_overlay = nreps %ORifNULL% 50,
    hist = nreps %ORifNULL% 8,
    dens = nreps %ORifNULL% 8,
    boxplot = nreps %ORifNULL% 8,
    freqpoly = nreps %ORifNULL% 8,
    freqpoly_grouped = nreps %ORifNULL% 3,
    violin_grouped = nreps, 
    error_binned = nreps %ORifNULL% 3,
    error_hist = nreps %ORifNULL% 3,
    error_hist_grouped = nreps %ORifNULL% 3,
    error_scatter = nreps %ORifNULL% 3,
    error_scatter_avg = nreps,
    error_scatter_avg_vs_x = nreps,
    scatter = nreps %ORifNULL% 3,
    scatter_avg = nreps,
    scatter_avg_grouped = nreps, 
    stat = rsa_ignore_nreps(nreps),
    stat_2d = rsa_ignore_nreps(nreps), 
    stat_grouped = rsa_ignore_nreps(nreps),
    stat_freqpoly_grouped = rsa_ignore_nreps(nreps), 
    intervals = rsa_ignore_nreps(nreps),
    intervals_grouped = rsa_ignore_nreps(nreps), 
    ribbon = rsa_ignore_nreps(nreps),
    ribbon_grouped = rsa_ignore_nreps(nreps), 
    rootogram = nreps,
    bars = nreps,
    bars_grouped = nreps, 
    loo_pit = rsa_ignore_nreps(nreps),
    loo_intervals = rsa_ignore_nreps(nreps), 
    loo_ribbon = rsa_ignore_nreps(nreps),
    stop("Plotting function not supported. ", 
      "(If the plotting function is included in the output from ", 
      "bayesplot::available_ppc() then it should be available via pp_check ", 
      "and this error is probably a bug.)")
  )
}


`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) return(b)
  return(a)
}


rsa_ignore_nreps <- function(nreps) {
  if (!is.null(nreps)) warning("'nreps' is ignored for this PPC", call. = FALSE)
  return(NULL)
}


nauf_ppc_args <- function(object, y, yrep, fun, ...) {
  funname <- fun
  fun <- match.fun(fun)
  dots <- list(...)
  dots[["y"]] <- as.numeric(y)
  dots[["yrep"]] <- yrep
  argnames <- names(formals(fun))
  
  if ("group" %in% argnames) {
    groupvar <- dots[["group"]] %ORifNULL% stop("This PPC requires the 'group'",
      " argument.", call. = FALSE)
    dots[["group"]] <- rsa_ppc_groupvar(object, groupvar)
  }
  
  if ("x" %in% argnames) {
    xvar <- dots[["x"]]
    if (!is.null(xvar)) {
      dots[["x"]] <- rsa_ppc_xvar(object, xvar)
    } else {
      if (funname %in% c("ppc_intervals", "ppc_ribbon")) {
        message("'x' not specified in '...'. Using x=1:length(y).")
        dots[["x"]] <- seq_along(y)
      } else {
        stop("This PPC requires the 'x' argument.", call. = FALSE)
      }
    }
  }
  
  if ("lw" %in% argnames && is.null(dots[["lw"]])) {
    dots[["lw"]] <- nauf_loo_weights(object, log = TRUE)
  }
  
  return(dots)
}


rsa_ppc_groupvar <- function(object, var = NULL) {
  if (is.null(var) || !is.character(var)) return(var)
  mf <- model.frame(object)
  vars <- colnames(mf)
  if (var %in% vars) return(mf[, var])
  stop("Variable '", var, "' not found in model frame. ")
}


rsa_ppc_xvar <- function(object, var = NULL) {
  if (is.null(var) || !is.character(var)) return(var)
  mf <- model.frame(object)
  vars <- colnames(mf)
  if (var %in% vars) return(mf[, var])
  stop("Variable '", var, "' not found in model frame. ")
}

