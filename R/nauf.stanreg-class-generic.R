

#' Class for fitted Bayesian models with \code{nauf} contrasts.
#'
#' ADD DESCRIPTION
#'
#' @seealso \code{\link{nauf_stan_glm}}, \code{\link{nauf_stan_glmer}},
#'   \code{\link{nauf_contrasts}}, and \code{\link[rstanarm]{stanreg-objects}}.
#'
#' @name nauf.stanreg
NULL


#' @export
formula.nauf.stanreg <- function(x, fixed.only = FALSE, random.only = FALSE,
                                 ...) {
  form <- x$formula
  
  if (is.nauf.stanmer(x)) {
    if (missing(fixed.only) && random.only) fixed.only <- FALSE
    if (fixed.only && random.only) {
      stop("can't specify 'only fixed' and 'only random' terms")
    }
    if (fixed.only) form <- lme4_getFixedFormula(form)
    if (random.only) form <- lme4_reOnly(form, response = TRUE)
  }
  
  class(form) <- c("nauf.formula", "formula")

  return(form)
}


#' @export
terms.nauf.stanreg <- function(x, fixed.only = TRUE, random.only = FALSE, ...) {
  if (is.nauf.stanfer(x)) return(x$terms)
  
  if (missing(fixed.only) && random.only) fixed.only <- FALSE
  if (fixed.only && random.only) {
    stop("can't specify 'only fixed' and 'only random' terms")
  }

  tt <- attr(x$glmod$fr, "terms")

  if (fixed.only) {
    tt <- terms.formula(formula(x, fixed.only = TRUE))
    attr(tt, "predvars") <- attr(terms(x$glmod$fr), "predvars.fixed")
  }
  if (random.only) {
    tt <- terms.formula(lme4::subbars(formula(x, random.only = TRUE)))
    attr(tt, "predvars") <- attr(terms(x$glmod$fr), "predvars.random")
  }

  attr(tt, "nauf.info") <- nauf.info(x$glmod$fr)
    # calling nauf.info(x) would create infinite recursion
  class(tt) <- c("nauf.terms", "terms", "formula")

  return(tt)
}


#' @export
print.nauf.stanreg <- function(x, ...) {
  NextMethod("print")
  cat("\n")
}


#' @export
predict.nauf.stanreg <- function(object, newdata = NULL, re.form = NULL,
                                 type = c("link", "response"), offset = NULL,
                                 na.action = na.pass, samples = FALSE, ...) {
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("'newdata', if specified, must be a data.frame.")
  }
  type <- match.arg(type)
  
  dat <- nauf_pp_data(object, newdata = newdata, re.form = re.form,
    offset = offset, ...)
    
  eta <- nauf_pp_eta(object, dat)
  if (type == "response") eta <- object$family$linkinv(eta)
  
  fit <- colMeans(eta)
  if (is.null(newdata)) {
    names(fit) <- nms <- rownames(newdata)
  } else {
    names(fit) <- nms <- rownames(model.frame(object))
  }
  
  if (!samples) return(fit)
  
  eta <- nauf.mcmc(eta, nchain(object), nms)
  
  return(list(fit = fit, samples = eta))
}


#' @export
posterior_predict.nauf.stanreg <- function(object, newdata = NULL, draws = NULL,
                                           re.form = NULL, fun = NULL,
                                           seed = NULL, offset = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(fun)) fun <- match.fun(fun)
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("'newdata', if specified, must be a data.frame.")
  }
  
  dat <- nauf_pp_data(object, newdata = newdata, re.form = re.form,
    offset = offset, ...)
    
  ppargs <- nauf_pp_args(object, data = nauf_pp_eta(object, dat, draws))
  if (object$family$family == "binomial") {
    ppargs$trials <- rep(1, nrow(dat$x))
  }
  
  ppfun <- rsa_pp_fun(object)
  ytilde <- do.call(ppfun, ppargs)
  
  if (!is.null(newdata) && nrow(newdata) == 1L) ytilde <- t(ytilde)
  if (!is.null(fun)) ytilde <- do.call(fun, list(ytilde))
  if (is.null(newdata)) {
    colnames(ytilde) <- rownames(model.frame(object))
  } else {
    colnames(ytilde) <- rownames(newdata)
  }
  
  return(structure(ytilde, class = c("ppd", class(ytilde))))
}


#' @export
log_lik.nauf.stanreg <- function(object, newdata = NULL, offset = NULL, ...) {
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("'newdata', if specified, must be a data.frame.")
  }
  
  calling_fun <- as.character(sys.call(-1))[1]
  args <- nauf_ll_args(object, newdata = newdata, offset = offset, 
    reloo_or_kfold = calling_fun %in% c("kfold", "reloo"), ...)
    
  fun <- rsa_ll_fun(object)
  
  out <- vapply(seq_len(args$N), FUN = function(i) {
    as.vector(fun(i = i, data = args$data[i, , drop = FALSE],
      draws = args$draws))
  }, FUN.VALUE = numeric(length = args$S))
  
  if (is.null(newdata)) {
    colnames(out) <- rownames(model.frame(object))
  } else {
    colnames(out) <- rownames(newdata)
  }
  
  return(out)
}


#' @export
fixef.nauf.stanreg <- function(object, method = c("median", "mean", "samples"),
                               permuted = TRUE, ...) {
  method <- match.arg(method)
  
  if (method == "median") return(NextMethod("fixef"))
  
  if (method == "mean") {
    return(object$stan_summary[1:ncol(model.matrix(object)), "mean"])
  }
  
  samp <- rstan::extract(object$stanfit,
    c(if (colnames(model.matrix(object)) == "(Intercept)") "alpha", "beta"),
    permuted)
  
  if (permuted) return(do.call(cbind, samp))
  
  return(nauf.mcmc(samp))
}


#' @export
ranef.nauf.stanreg <- function(object, method = c("median", "mean", "samples"),
                               permuted = TRUE, by_group = TRUE,
                               drop_NEW = TRUE, ...) {
  method <- match.arg(method)
  
  all_names <- object$stanfit@sim$fnames_oi
  sel <- rsa_b_names(all_names)
  
  fl <- as.list(object$glmod$reTrms$flist)
  levs <- lapply(fl, levels)
  
  if (drop_NEW) {
    not_NEW <- which(!grepl("_NEW_", all_names[sel], fixed = TRUE))
    sel <- setdiff(sel, which(grepl("_NEW_", all_names, fixed = TRUE)))
  } else {
    levs <- lapply(levs, c, "_NEW_")
  }
  
  bnms <- all_names[sel]
  bnms <- substr(bnms, 3, nchar(bnms) - 1)
  
  asgn <- attr(fl, "assign")
  attr(fl, "assign") <- NULL
  cnms <- object$glmod$reTrms$cnms
  mark <- !grepl("^Xr", names(cnms))
  cnms <- cnms[mark]
  asgn <- asgn[mark]
  
  nc <- vapply(cnms, length, 1L)
  nm <- vapply(levs[asgn], length, 1L)
  nb <- nc * nm
  nbseq <- rep.int(seq_along(nb), nb)
    
  if (method != "samples") {
    ans <- object$stan_summary[sel, if (method == "mean") "mean" else "50%"]
    
    if (!by_group) {
      names(ans) <- bnms
      return(ans)
    }
    
    ml <- split(ans, nbseq)
    for (i in seq_along(ml)) {
      ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
        dimnames = list(NULL, cnms[[i]]))
    }
    ans <- sapply(intersect(names(fl), names(cnms)), function(i) {
      data.frame(do.call(cbind, ml[names(ml) == i]), row.names = levs[[i]],
        check.names = FALSE)
    }, simplify = FALSE)
    
    return(structure(ans, class = "ranef.mer"))
  }
  
  if (permuted) {
    ans <- rstan::extract(object$stanfit, "b")[[1]]
    if (drop_NEW) ans <- ans[, not_NEW, drop = FALSE]
    
  } else if (!by_group) {
    return(nauf.mcmc(rstan::extract(object$stanfit, "b",
      FALSE)[, , sel, drop = FALSE], NULL, bnms))
      
  } else {
    ans <- as.matrix(object$stanfit)[, sel, drop = FALSE]
  }
  
  #  dimnames(ans) <- list(iterations = NULL, parameters = bnms)
  if (!by_group) return(ans)
  
  
}

function (object, ...) 
{
    all_names <- if (used.optimizing(object)) 
        rownames(object$stan_summary)
    else object$stanfit@sim$fnames_oi
    sel <- b_names(all_names)
    ans <- object$stan_summary[sel, select_median(object$algorithm)]
    ans <- ans[!grepl("_NEW_", names(ans), fixed = TRUE)]
    fl <- .flist(object)
    levs <- lapply(fl, levels)
    asgn <- attr(fl, "assign")
    cnms <- .cnms(object)
    mark <- !grepl("^Xr", names(cnms))
    fl <- fl[mark]
    asgn <- asgn[mark]
    levs <- levs[mark]
    cnms <- cnms[mark]
    nc <- vapply(cnms, length, 1L)
    nb <- nc * vapply(levs, length, 1L)
    
    nbseq <- rep.int(seq_along(nb), nb)
    ml <- split(ans, nbseq)
    for (i in seq_along(ml)) {
        ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE, 
            dimnames = list(NULL, cnms[[i]]))
    }
    ans <- lapply(seq_along(fl), function(i) {
        data.frame(do.call(cbind, ml[i]), row.names = levs[[i]], 
            check.names = FALSE)
    })
    names(ans) <- names(fl)
    structure(ans, class = "ranef.mer")
}



#' @export
ranef.nauf.stanreg <- function(object, samples = FALSE, permuted = TRUE,
                               arr = FALSE, drop_NEW = TRUE, ...) {
  if (!is.nauf.stanmer(object)) {
    stop("model contains no random effects")
  }
  if (!samples) return(NextMethod("ranef"))
  
  g <- rstan::extract(object$stanfit, "b", permuted = permuted)
  all_names <- object$stanfit@sim$fnames_oi
  b_names <- all_names[grep("^b\\[", all_names)]
  keep <- !grepl("_NEW_", b_names, fixed = TRUE)
  
  if (drop_NEW) {
    if (permuted) {
      g <- g[[1]][, keep, drop = FALSE]
    } else {
      g <- g[, , keep, drop = FALSE]
    }
  } else if (permuted) {
    g <- g[[1]]
  }
  
  re <- object$glmod$reTrms
  
  groups <- colnames(re$flist)
  members <- sapply(re$flist, levels, simplify = FALSE)
  if (!drop_NEW) {
    for (g in groups) {
      members[[g]] <- c(members[[g]], "_NEW_")
    }
  }
  
  if (arr) {
    gnums <- vector("list", length(groups))
    names(gnums) <- groups
    nms <- gnums
    
    asgn <- attr(re$flist, "assign")
    cnms <- re$cnms
    nc <- lengths(cnms)
    nm <- lengths(members)[asgn]
    nb <- nc * nm
    bnums <- lapply(nb, function(x) 1:x)
    if (length(bnums) > 1) {
      for (j in 2:length(bnums)) {
        bnums[[j]] <- bnums[[j]] + max(bnums[[j - 1]])
      }
    }
    mnums <- list()
    for (j in 1:length(bnums)) {
      mnums[[j]] <- rep(1:nm[j], each = nc[j])
    }
    
    for (j in 1:length(groups)) {
      w <- which(asgn == j)
      gnums[[j]] <- unname(unlist(bnums[w]))[order(unname(unlist(mnums[w])))]
      nms[[j]] <- list(m = members[[j]], parameters = unname(unlist(cnms[w])))
      names(nms[[j]])[1] <- groups[j]
    }
    
    matdims <- lapply(nms, lengths)
    
    if (permuted) {
      return(mapply(function(loc, md, g, n) {
        iter <- dim(g)[1]
        a <- array(dim = c(iter, md), dimnames = c(dimnames(g)[1], n))
        for (i in 1:iter) {
          a[i, , ] <- matrix(g[i, loc], md[1], md[2], byrow = TRUE)
        }
        return(a)
      }, gnums, matdims, nms, MoreArgs = list(g = g), SIMPLIFY = FALSE))
      
    } else {
      return(mapply(function(loc, md, g, n) {
        iter <- dim(g)[1]
        chains <- dim(g)[2]
        a <- array(dim = c(iter, chains, md), dimnames = c(dimnames(g)[1:2], n))
        for (ch in 1:chains) {
          for (i in 1:iter) {
            a[i, ch, , ] <- matrix(g[i, ch, loc], md[1], md[2], byrow = TRUE)
          }
        }
        return(a)
      }, gnums, matdims, nms, MoreArgs = list(g = g), SIMPLIFY = FALSE))
    }
  }
  
  nms <- sapply(groups, function(x) paste(x, members[[x]], sep = "_"),
    simplify = FALSE)
  nms <- unlist(sapply(1:length(re$cnms), function(x)
    as.vector(t(outer(nms[[names(re$cnms)[x]]], re$cnms[[x]], paste,
    sep = ":"))), simplify = FALSE))
    
  if (permuted) {
    dimnames(g) <- list(iterations = NULL, parameters = nms)
  } else {
    dimnames(g)[[3]] <- nms
  }
  
  return(g)
}


#' @export
coef.nauf.stanreg <- function(object, samples = FALSE, permuted = TRUE,
                              drop_NEW = TRUE, ...) {
  if (!samples) return(NextMethod("coef"))
  
  fef <- fixef(object, samples, permuted)
  if (!inherits(object, "lmerMod")) return(fe)
  ref <- ranef(object, samples, permuted, arr = TRUE, drop_NEW = drop_NEW)
  re <- object$glmod$reTrms
  fl <- re$flist
  groups <- colnames(fl)
  glist <- vector("list", ncol(fl))
  names(glist) <- groups
  
  if (permuted) {
    fnms <- colnames(fef)
    for (g in groups) {
      rnms <- dimnames(ref[[g]])[[3]]
      dnms <- dimnames(ref[[g]])
      dnms[[3]] <- union(fnms, rnms)
      
      glist[[g]] <- array(0, dim = c(dim(ref[[g]])[1:2], length(dnms[[3]])),
        dimnames = dnms)
      glist[[g]][, , rnms] <- ref[[g]][, , rnms]
      
      for (f in fnms) {
        glist[[g]][, , f] <- glist[[g]][, , f] + matrix(fef[, f],
          nrow(fef), dim(ref[[g]])[2])
      }
    }
    
    return(glist)
  }
  
  fnms <- dimnames(fef)[[3]]
  for (g in groups) {
    rnms <- dimnames(ref[[g]])[[4]]
    dnms <- dimnames(ref[[g]])
    dnms[[4]] <- union(fnms, rnms)
    
    glist[[g]] <- array(0, dim = c(dim(ref[[g]])[1:3], length(dnms[[4]])),
      dimnames = dnms)
    glist[[g]][, , , rnms] <- ref[[g]][, , , rnms]
    
    for (f in fnms) {
      for (ch in 1:dim(ref[[g]])[2]) {
        glist[[g]][, ch, , f] <- glist[[g]][, ch, , f] + matrix(fef[, ch, f],
          dim(fef)[1], dim(ref[[g]])[3])
      }
    }
  }
  
  return(glist)
}

