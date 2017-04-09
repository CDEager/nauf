

## this file contains internal functions taken from lme4 v1.1-12
## which are altered as little as possible and preceded by 'lme4_'


lme4_safeDeparse <- function(x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}


lme4_barnames <- function(bars) {
  vapply(bars, function(x) lme4_safeDeparse(x[[3]]), "")
}


lme4_wmsg <- function(n, cmp.val, allow.n, msg1 = "", msg2 = "", msg3 = "") {
  if (allow.n) {
    unident <- n < cmp.val
    cmp <- "<"
    rstr <- ""
  } else {
    unident <- n <= cmp.val
    cmp <- "<="
    rstr <- " and the residual variance (or scale parameter)"
  }
  wstr <- sprintf("%s (=%d) %s %s (=%d)%s; the random-effects parameters%s are probably unidentifiable", 
    msg1, n, cmp, msg2, cmp.val, msg3, rstr)
    
  list(unident = unident, wstr = wstr)
}


lme4_RHSForm <- function(form, as.form = FALSE) {
  rhsf <- form[[length(form)]]
  if (as.form) return(stats::reformulate(deparse(rhsf)))
  return(rhsf)
}


lme4_`RHSForm<-` <- function(formula, value) {
  formula[[length(formula)]] <- value
  formula
}


lme4_getFixedFormula <- function(form) {
  lme4_RHSForm(form) <- lme4::nobars(lme4_RHSForm(form))
  form
}


lme4_reOnly <- function(f, response = FALSE) {
  if (response && length(f) == 3) {
    response <- f[[2]]
  } else {
    response <- NULL
  }
  stats::reformulate(paste0("(", vapply(lme4::findbars(f), lme4_safeDeparse,
    ""), ")"), response = response)
}


lme4_reFormHack <- function(re.form, ReForm, REForm, REform) {
  warnDeprec <- function(name) {
    warning(gettextf("'%s' is deprecated; use '%s' instead", name, "re.form"),
      call. = FALSE, domain = NA)
  }
  if (!missing(ReForm)) {
    warnDeprec("ReForm")
    return(ReForm)
  }
  if (!missing(REForm)) {
    warnDeprec("REForm")
    return(REForm)
  }
  if (!missing(REform)) {
    warnDeprec("REform")
    return(REform)
  }
  re.form
}


lme4_isGLMM <- function(x, ...) {
  as.logical(x@devcomp$dims[["GLMM"]])
}


lme4_isNLMM <- function(x, ...) {
  as.logical(x@devcomp$dims[["NLMM"]])
}


lme4_isLMM <- function(x, ...) {
  !lme4_isGLMM(x) && !lme4_isNLMM(x)
}


lme4_checkArgs <- function(type, ...) {
  l... <- list(...)
  if (isTRUE(l...[["sparseX"]])) {
    warning("sparseX = TRUE has no effect at present", call. = FALSE)
  }
  
  if (length(l... <- list(...))) {
    if (!is.null(l...[["family"]])) {
      warning("calling lmer with family() is deprecated: please use glmer() instead", 
        call. = FALSE)
      type <- "glmer"
    }
    
    if (!is.null(l...[["method"]])) {
      msg <- paste("Argument", sQuote("method"), "is deprecated.")
      if (type == "lmer") {
        msg <- paste(msg, "Use the REML argument to specify ML or REML estimation.")
      } else if (type == "glmer") {
        msg <- paste(msg, "Use the nAGQ argument to specify Laplace (nAGQ=1) or adaptive", 
          "Gauss-Hermite quadrature (nAGQ>1).  PQL is no longer available.")
      }
      warning(msg, call. = FALSE)
      l... <- l...[names(l...) != "method"]
    }
    
    if (length(l...)) {
      warning("extra argument(s) ", paste(sQuote(names(l...)),
        collapse = ", "), " disregarded", call. = FALSE)
    }
  }
}


lme4_checkCtrlLevels <- function(cstr, val, smallOK = FALSE) {
  bvals <- c("stop", "warning", "ignore")
  if (smallOK) bvals <- outer(bvals, c("", "Small"), paste0)
  if (!is.null(val) && !val %in% bvals) {
    stop("invalid control level ", sQuote(val), " in ", cstr, 
      ": valid options are {", paste(sapply(bvals, sQuote),
      collapse = ","), "}")
  }
  invisible(NULL)
}


lme4_checkFormulaData <- function(formula, data, checkLHS = TRUE,
                                  debug = FALSE) {
  nonexist.data <- missDataFun(data)
  wd <- tryCatch(eval(data), error = identity)
  if (wrong.data <- inherits(wd, "simpleError")) {
    wrong.data.msg <- wd$message
  }
  bad.data <- nonexist.data || wrong.data
  
  if (bad.data || debug) {
    varex <- function(v, env) exists(v, envir = env, inherits = FALSE)
    allvars <- all.vars(stats::as.formula(formula))
    allvarex <- function(env, vvec = allvars) all(vapply(vvec, varex, NA, env))
  }
  
  if (bad.data) {
    if (allvarex(environment(formula))) {
      stop("bad 'data', but variables found in environment of formula: ", 
        "try specifying 'formula' as a formula rather ", 
        "than a string in the original model", call. = FALSE)
    } else {
      if (nonexist.data) {
        stop("'data' not found, and some variables missing from formula environment", 
          call. = FALSE)
      }
      else {
        stop("bad 'data': ", wrong.data.msg)
      }
    }
  } else {
    if (is.null(data)) {
      if (!is.null(ee <- environment(formula))) {
        denv <- ee
      } else {
        denv <- parent.frame(2L)
      }
    } else {
      denv <- list2env(data)
    }
  }
  
  if (debug) {
    cat("Debugging parent frames in checkFormulaData:\n")
    glEnv <- 1L
    while (!identical(parent.frame(glEnv), .GlobalEnv)) {
      glEnv <- glEnv + 1L
    }
    for (i in 1:glEnv) {
      OK <- allvarex(parent.frame(i))
      cat("vars exist in parent frame ", i)
      if (i == glEnv) cat(" (global)")
      cat(" ", OK, "\n")
    }
    cat("vars exist in env of formula ", allvarex(denv), "\n")
  }
  
  stopifnot(!checkLHS || length(as.formula(formula, env = denv)) == 3)
  
  return(denv)
}


lme4_doCheck <- function(x) {
  is.character(x) && !any(x == "ignore")
}


lme4_checkNlevels <- function(flist, n, ctrl, allow.n = FALSE) {
  stopifnot(is.list(ctrl), is.numeric(n))
  nlevelVec <- unlist(lapply(flist, function(x) nlevels(droplevels(x))))
  
  cstr <- "check.nlev.gtr.1"
  lme4_checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (lme4_doCheck(cc) && any(nlevelVec < 2)) {
    wstr <- "grouping factors must have > 1 sampled level"
    switch(cc, warning = warning(wstr, call. = FALSE), stop = stop(wstr, 
      call. = FALSE), stop(gettextf("unknown check level for '%s'", 
      cstr), domain = NA))
  } else {
    wstr <- character()
  }
  
  cstr <- "check.nobs.vs.nlev"
  lme4_checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (lme4_doCheck(cc) && any(if (allow.n) nlevelVec > n else nlevelVec >= n)) {
    wst2 <- gettextf("number of levels of each grouping factor must be %s number of observations", 
      if (allow.n) "<=" else "<")
    switch(cc, warning = warning(wst2, call. = FALSE), stop = stop(wst2, 
      call. = FALSE))
  } else {
    wst2 <- character()
  }
  
  cstr <- "check.nlev.gtreq.5"
  lme4_checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (lme4_doCheck(cc) && any(nlevelVec < 5)) {
    wst3 <- "grouping factors with < 5 sampled levels may give unreliable estimates"
    switch(cc, warning = warning(wst3, call. = FALSE), stop = stop(wst3, 
      call. = FALSE), stop(gettextf("unknown check level for '%s'", 
      cstr), domain = NA))
  } else {
    wst3 <- character()
  }
  
  c(wstr, wst2, wst3)
}


lme4_checkZdims <- function(Ztlist, n, ctrl, allow.n = FALSE) {
  stopifnot(is.list(Ztlist), is.numeric(n))
  cstr <- "check.nobs.vs.nRE"
  lme4_checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  term.names <- names(Ztlist)
  rows <- vapply(Ztlist, nrow, 1L)
  cols <- vapply(Ztlist, ncol, 1L)
  stopifnot(all(cols == n))
  
  if (lme4_doCheck(cc)) {
    unique(unlist(lapply(seq_along(Ztlist), function(i) {
      ww <- lme4_wmsg(cols[i], rows[i], allow.n, "number of observations", 
        "number of random effects", sprintf(" for term (%s)", term.names[i]))
      if (ww$unident) {
        switch(cc, warning = warning(ww$wstr, call. = FALSE), 
          stop = stop(ww$wstr, call. = FALSE), stop(gettextf("unknown check level for '%s'", 
            cstr), domain = NA))
        ww$wstr
      } else {
        character()
      }
    })))
  } else {
    character()
  }
}


lme4_checkZrank <- function(Zt, n, ctrl, nonSmall = 1e+06, allow.n = FALSE) {
  stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
  cstr <- "check.nobs.vs.rankZ"
  if (lme4_doCheck(cc <- ctrl[[cstr]])) {
    lme4_checkCtrlLevels(cstr, cc, smallOK = TRUE)
    d <- dim(Zt)
    doTr <- d[1L] < d[2L]
    if (!(grepl("Small", cc) && prod(d) > nonSmall)) {
      rankZ <- Matrix::rankMatrix(if (doTr) Matrix::t(Zt) else Zt,
        method = "qr", sval = numeric(min(d)))
      ww <- lme4_wmsg(n, rankZ, allow.n, "number of observations", "rank(Z)")
      if (ww$unident) {
        switch(cc, warningSmall = , warning = warning(ww$wstr, 
          call. = FALSE), stopSmall = , stop = stop(ww$wstr, 
          call. = FALSE), stop(gettextf("unknown check level for '%s'", 
          cstr), domain = NA))
        ww$wstr
      } else {
        character()
      }
    } else {
      character()
    }
  } else {
    character()
  }
}


lme4_chkRank.drop.cols <- function(X, kind, tol = 1e-07, method = "qr.R") {
  stopifnot(is.matrix(X))
  kinds <- eval(formals(lme4::lmerControl)[["check.rankX"]])
  if (!kind %in% kinds) stop(gettextf("undefined option for 'kind': %s", kind))
  if (kind == "ignore") return(X)
  p <- ncol(X)
  if (kind == "stop.deficient") {
    if ((rX <- Matrix::rankMatrix(X, tol = tol, method = method)) < p) {
      stop(gettextf(sub("\n +", "\n", "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);\n                   the fixed effects will be jointly unidentifiable"), rX, p), call. = FALSE)
    }
  } else {
    qr.X <- qr(X, tol = tol, LAPACK = FALSE)
    rnkX <- qr.X$rank
    if (rnkX == p) return(X)
    
    msg <- sprintf(ngettext(p - rnkX, "fixed-effect model matrix is rank deficient so dropping %d column / coefficient", 
      "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"), 
      p - rnkX)
    if (kind != "silent.drop.cols") {
      if (kind == "warn+drop.cols") {
        warning(msg, domain = NA)
      } else {
        message(msg, domain = NA)
      }
    }
    
    contr <- attr(X, "contrasts")
    asgn <- attr(X, "assign")
    keep <- qr.X$pivot[seq_len(rnkX)]
    dropped.names <- colnames(X[, -keep, drop = FALSE])
    X <- X[, keep, drop = FALSE]
    
    if (Matrix::rankMatrix(X, tol = tol, method = method) < ncol(X)) {
      stop(gettextf("Dropping columns failed to produce full column rank design matrix"), 
        call. = FALSE)
    }
    if (!is.null(contr)) attr(X, "contrasts") <- contr
    if (!is.null(asgn)) attr(X, "assign") <- asgn[keep]
    attr(X, "msgRankdrop") <- msg
    attr(X, "col.dropped") <- stats::setNames(qr.X$pivot[(rnkX + 1L):p],
      dropped.names)
  }
  
  X
}


lme4_checkScaleX <- function (X, kind = "warning", tol = 1000) {
  kinds <- eval(formals(lme4::lmerControl)[["check.scaleX"]])
  if (!kind %in% kinds) stop(gettextf("unknown check-scale option: %s", kind))
  if (is.null(kind) || kind == "ignore") return(X)
  
  cont.cols <- apply(X, 2, function(z) !all(z %in% c(0, 1)))
  col.sd <- apply(X[, cont.cols, drop = FALSE], 2L, sd)
  sdcomp <- outer(col.sd, col.sd, "/")
  logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
  logsd <- abs(log(col.sd))
  
  if (any(c(logcomp, logsd) > log(tol))) {
    wmsg <- "Some predictor variables are on very different scales:"
    if (kind %in% c("warning", "stop")) {
      wmsg <- paste(wmsg, "consider rescaling")
      switch(kind, warning = warning(wmsg, call. = FALSE), 
        stop = stop(wmsg, call. = FALSE))
    } else {
      wmsg <- paste(wmsg, "auto-rescaled (results NOT adjusted)")
      X[, cont.cols] <- sweep(X[, cont.cols, drop = FALSE], 2, col.sd, "/")
      attr(X, "scaled:scale") <- stats::setNames(col.sd, colnames(X)[cont.cols])
      if (kind == "warn+rescale") warning(wmsg, call. = FALSE)
    }
  } else {
    wmsg <- character()
  }
  
  structure(X, msgScaleX = wmsg)
}


lme4_checkHess <- function(H, tol, hesstype = "") {
  res <- list(code = numeric(0), messages = list())
  
  evd <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, 
      error = function(e) e)
      
  if (inherits(evd, "error")) {
      res$code <- -6L
      res$messages <- gettextf("Problem with Hessian check (infinite or missing values?)")
      
  } else {
    negative <- sum(evd < -tol)
    if (negative) {
      res$code <- -3L
      res$messages <- gettextf(paste("Model failed to converge:", 
        "degenerate", hesstype, "Hessian with %d negative eigenvalues"),
        negative)
        
    } else {
      zero <- sum(abs(evd) < tol)
      if (zero || inherits(tryCatch(chol(H), error = function(e) e), "error")) {
        res$code <- -4L
        res$messages <- paste(hesstype, "Hessian is numerically singular: parameters are not uniquely determined")
        
      } else {
        res$cond.H <- max(evd) / min(evd)
        if (max(evd) * tol > 1) {
          res$code <- c(res$code, 2L)
          res$messages <- c(res$messages, paste("Model is nearly unidentifiable: ", 
            "very large eigenvalue", "\n - Rescale variables?", sep = ""))
        }
        
        if ((min(evd)/max(evd)) < tol) {
          res$code <- c(res$code, 3L)
          if (!5L %in% res$code) {
            res$messages <- c(res$messages, paste("Model is nearly unidentifiable: ", 
              "large eigenvalue ratio", "\n - Rescale variables?", sep = ""))
          }
        }
      }
    }
  }
  
  if (length(res$code) == 0) res$code <- 0
  
  res
}


lme4_checkConv <- function(derivs, coefs, ctrl, lbound, debug = FALSE) {
  if (is.null(derivs)) return(NULL)
  if (anyNA(derivs$gradient)) {
    return(list(code = -5L, messages = gettextf("Gradient contains NAs")))
  }
  
  ntheta <- length(lbound)
  res <- list()
  
  ccl <- ctrl[[cstr <- "check.conv.grad"]]
  lme4_checkCtrlLevels(cstr, cc <- ccl[["action"]])
  wstr <- NULL
  if (lme4_doCheck(cc)) {
    scgrad <- tryCatch(with(derivs, solve(chol(Hessian), 
      gradient)), error = function(e) e)
      
    if (inherits(scgrad, "error")) {
      wstr <- "unable to evaluate scaled gradient"
      res$code <- -1L
    } else {
      mingrad <- pmin(abs(scgrad), abs(derivs$gradient))
      maxmingrad <- max(mingrad)
      if (maxmingrad > ccl$tol) {
        w <- which.max(maxmingrad)
        res$code <- -1L
        wstr <- gettextf("Model failed to converge with max|grad| = %g (tol = %g, component %d)", 
          maxmingrad, ccl$tol, w)
      }
    }
    
    if (!is.null(wstr)) {
      res$messages <- wstr
      switch(cc, warning = warning(wstr), stop = stop(wstr), 
        stop(gettextf("unknown check level for '%s'", cstr), domain = NA))
    }
    
    if (!is.null(ccl$relTol) && (max.rel.grad <- max(abs(derivs$gradient/coefs))) > ccl$relTol) {
      res$code <- -2L
      wstr <- gettextf("Model failed to converge with max|relative grad| = %g (tol = %g)", 
        max.rel.grad, ccl$relTol)
      res$messages <- wstr
      switch(cc, warning = warning(wstr), stop = stop(wstr), 
        stop(gettextf("unknown check level for '%s'", cstr), domain = NA))
    }
  }
  
  ccl <- ctrl[[cstr <- "check.conv.singular"]]
  lme4_checkCtrlLevels(cstr, cc <- ccl[["action"]])
  if (lme4_doCheck(cc)) {
    bcoefs <- seq(ntheta)[lbound == 0]
    if (any(coefs[bcoefs] < ccl$tol)) {
      wstr <- "singular fit"
      res$messages <- c(res$messages, wstr)
      switch(cc, warning = warning(wstr), stop = stop(wstr), 
        stop(gettextf("unknown check level for '%s'", cstr), domain = NA))
    }
  }
  
  ccl <- ctrl[[cstr <- "check.conv.hess"]]
  lme4_checkCtrlLevels(cstr, cc <- ccl[["action"]])
  if (lme4_doCheck(cc)) {
    if (length(coefs) > ntheta) {
      H.beta <- derivs$Hessian[-seq(ntheta), -seq(ntheta)]
      resHess <- lme4_checkHess(H.beta, ccl$tol, "fixed-effect")
      if (any(resHess$code != 0)) {
        res$code <- resHess$code
        res$messages <- c(res$messages, resHess$messages)
        wstr <- paste(resHess$messages, collapse = ";")
        switch(cc, warning = warning(wstr), stop = stop(wstr), 
          stop(gettextf("unknown check level for '%s'", cstr), domain = NA))
      }
    }
    
    resHess <- lme4_checkHess(derivs$Hessian, ccl$tol)
    if (any(resHess$code != 0)) {
      res$code <- resHess$code
      res$messages <- c(res$messages, resHess$messages)
      wstr <- paste(resHess$messages, collapse = ";")
      switch(cc, warning = warning(wstr), stop = stop(wstr), 
        stop(gettextf("unknown check level for '%s'", cstr), domain = NA))
    }
  }
  
  if (debug && length(res$messages) > 0) print(res$messages)
  
  res
}


lme4_updateStart <- function(start, theta) {
  if (is.null(start)) return(NULL)
  if (is.numeric(start)) {
    start <- theta
  } else if (!is.null(start$theta)) { 
    start$theta <- theta
  }
  start
}


lme4_est_theta <- function(object, limit = 20, eps = .Machine$double.eps^0.25, 
                           trace = 0) {
  Y <- stats::model.response(model.frame(object))
  mu <- fitted(object)
  MASS::theta.ml(Y, mu, weights = object@resp$weights, limit = limit, 
    eps = eps, trace = trace)
}


# Not sure if this one will work
lme4_setNBdisp <- function(object, theta) {
  rr <- object@resp
  newresp <- do.call(lme4::glmResp$new, c(lapply(
    stats::setNames(nm = glmNB.to.change), rr$field),
    list(family = MASS::negative.binomial(theta = theta))))
  newresp$setOffset(rr$offset)
  newresp$updateMu(rr$eta - rr$offset)
  object@resp <- newresp
  object
}


lme4_refitNB <- function(object, theta, control = NULL) {
  lme4::refit(lme4_setNBdisp(object, theta), control = control)
}


lme4_getNBdisp <- function(object) {
  environment(object@resp$family$aic)[[".Theta"]]
}


lme4_optTheta <- function(object, interval = c(-5, 5),
                          tol = .Machine$double.eps^0.25, verbose = FALSE,
                          control = NULL) {
  lastfit <- object
  it <- 0L
  
  NBfun <- function(t) {
    dev <- -2 * logLik(lastfit <<- lme4_refitNB(lastfit, theta = exp(t),
      control = control))
    it <<- it + 1L
    if (verbose) {
      cat(sprintf("%2d: th=%#15.10g, dev=%#14.8f, beta[1]=%#14.8f\n", 
        it, exp(t), dev, lastfit@beta[1]))
    }
    dev
  }
  
  optval <- stats::optimize(NBfun, interval = interval, tol = tol)
  stopifnot(all.equal(optval$minimum, log(lme4_getNBdisp(lastfit))))
  attr(lastfit, "nevals") <- it
  lastfit@call$family[["theta"]] <- exp(optval$minimum)
  
  lastfit
}

