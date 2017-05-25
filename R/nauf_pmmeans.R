

#' @rdname nauf-pmmeans
#'
#' @importFrom pbkrtest Lb_ddf vcovAdj
#' @importFrom lsmeans ref.grid
#' @importFrom lmerTest calcSatterth
#'
#' @export
nauf_ref.grid <- function(mod, KR = FALSE, ...) {
  if (!is.nauf.model(mod)) stop("Must supply a nauf model")

  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quotes(names(dots)))

  info <- nauf.info(mod)
  fenr <- stats::delete.response(terms(mod))
  vars <- varnms(fenr)
  mf <- model.frame(mod)
  
  uf <- intersect(vars, names(info$uf))
  xlev <- c(
    mlapply(fac = info$uf[uf], hasna = info$hasna[uf],
      fun = function(fac, hasna) c(fac[[1]], if (hasna) NA)),
    sapply(info$of[intersect(vars, names(info$of))], `[[`, "levels",
      simplify = FALSE)
  )
  lvs <- c(xlev, info$num[intersect(vars, names(info$num))])
  mlvs <- info$mat[intersect(vars, names(info$mat))]
  
  # reorder for consistency
  lvs <- lvs[intersect(vars, names(lvs))]
  xlev <- xlev[intersect(vars, names(xlev))]

  g <- expand.grid(lvs)
  g[mames(mlvs)] <- mlapply(data = mlvs, ncol = lengths(mlvs),
    same = list(nrow = nrow(g), byrow = TRUE), fun = matrix)

  mm <- model.matrix(fenr, g)
  asgn <- attr(mm, "assign")

  g[[".wgt."]] <- data.frame(xtabs(~ ., mf[names(xlev)]))$Freq

  if (!(bayes <- is.nauf.stanreg(mod))) summ <- summary(mod)

  temp_dat <- data.frame(y = 1:5, x = c(2, 5, 3, 6, 1))
  temp_mod <- lm(y ~ x, temp_dat)
  rg <- lsmeans::ref.grid(temp_mod, data = temp_dat)

  rg@model.info <- list(
    call = if (bayes) mod$call else summ$call,
    terms = fenr,
    xlev = xlev)
  rg@roles <- list(
    predictors = vars,
    responses = character(),
    multresp = character())
    
  rg@grid <- g
  rg@levels <- lvs
  rg@matlevs <- mlvs
  rg@linfct <- mm
  rg@bhat <- if (bayes) numeric() else summ$coefficients[, 1]
  rg@nbasis <- matrix()
  
  if (is.nauf.merMod(mod)) {
    rg@V <- as.matrix(summ$vcov)
  } else if (bayes) {
    rg@V <- matrix()
  } else {
    rg@V <- vcov(mod)
  }
  
  if (is.nauf.lmerMod(mod)) {
    if (!lme4::isREML(mod)) KR <- FALSE
    if (KR) {
      rg@dffun <- function(k, dfargs) {
        pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
      }
      rg@dfargs <- list(
        unadjV = vcov(mod),
        adjV = pbkrtest::vcovAdj(mod))
    } else {
      rg@dffun <- function(k, dfargs) {
        lmerTest::calcSatterth(dfargs$object, k)$denom
      }
      rg@dfargs <- list(object = mod)
    }
  } else if (is.nauf.glmerMod(mod) || bayes) {
    rg@dffun <- function(k, dfargs) NA
    rg@dfargs <- list()
  } else {
    rg@dffun <- function(k, dfargs) dfargs$df
    rg@dfargs <- list(df = mod$df)
  }
  
  rg@misc <- list(
    estName = "prediction",
    estType = "prediction",
    infer = c(FALSE, FALSE),
    level = 0.95,
    adjust = "none",
    famSize = ncol(g) - 1,
    avgd.over = character(),
    assign = asgn,
    chains = nchain(mod))
    
  if (bayes) {
    rg@post.beta <- as.matrix(mod)[, 1:ncol(model.matrix(mod)), drop = FALSE]
  }
  
  if (!is.linear(family <- get_family(mod))) {
    rg@misc$tran <- family$link
    rg@misc$family <- family
    family <- family$family
    
    rate <- length(grep("poisson", family)) + length(grep("Negative", family))
    prob <- length(grep("binomial", family))
    if (rate) {
      rg@misc$inv.lbl <- "rate"
    } else if (prob) {
      rg@misc$inv.lbl <- "prob"
    } else {
      rg@misc$inv.lbl <- "response"
    }
  } else {
    rg@misc$family <- family
  }

  return(structure(list(ref.grid = rg), class = "nauf.ref.grid"))
}


#' @rdname nauf-pmmeans
#'
#' @importFrom lsmeans contrast lsmeans pmmeans
#'
#' @export
nauf_pmmeans <- function(object, specs, pairwise = FALSE, subset = NULL,
                         na_as_level = NULL, by = NULL, ...) {
  if (!is.nauf.ref.grid(object)) stop("must supply a nauf.ref.grid")

  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quotes(names(dots)))
  
  rg <- object$ref.grid

  specs <- check_specs(specs, pairwise, rg@model.info$xlev, rg@model.info$terms,
    na_as_level, by)
  subset <- grid_subset(subset, specs)

  if (!is.null(subset$keep)) {
    keep <- eval(subset$keep, envir = rg@grid)
    rg@grid <- rg@grid[keep, , drop = FALSE]
    rg@linfct <- rg@linfct[keep, , drop = FALSE]
  }

  if (length(specs$facs)) {
    est_grid <- unique(rg@grid[, specs$facs, drop = FALSE])
    if (!nrow(est_grid)) {
      stop("Factors are given in 'specs' but no levels exist in the specified",
        " subset")
    }
    if (!is.null(subset$est_keep)) {
      keep <- eval(subset$est_keep, envir = est_grid)
      est_grid <- est_grid[keep, , drop = FALSE]
    }
    if (!nrow(est_grid)) {
      stop("Factors are given in 'specs' but no levels exist in the specified",
        " subset")
    }
  }

  if (length(specs$nums)) {
    if (!length(specs$facs) || !nrow(est_grid)) {
      est_grid <- data.frame(t(rep("inc_1", length(specs$nums))))
      colnames(est_grid) <- specs$nums
    } else {
      est_grid[specs$nums] <- "inc_1"
    }
    for (v in specs$nums) {
      rg@grid[[v]] <- rg@grid[[v]] + 1
    }
    rg@linfct <- model.matrix(rg@model.info$terms, rg@grid) - rg@linfct
  }
  
  if (!nrow(est_grid)) stop("No grid rows satisfy subsetting conditions")
  
  if (length(specs$by)) {
    by_fac <- est_grid[specs$by]
    for (j in specs$by) {
      by_fac[[j]] <- droplevels(by_fac[[j]])
      if (anyNA(by_fac[[j]])) by_fac[[j]] <- addNA(by_fac[[j]])
    }
    by_fac <- interaction(by_fac, drop = TRUE)
    if (any(xtabs(~ by_fac) < 2)) {
      cat("Not all combinations of ", paste(specs$by, collapse = ":"),
        " have more than 1 row in the subsetted estimate grid:\n\n", sep = "")
      print(est_grid)
      stop("Cannot condition pairwise comparisons on 'by'")
    }
  } else if (nrow(est_grid) == 1 && specs$pw) {
    cat("Specified 'pairwise' but subsetted estimate grid only has one row:\n\n")
    print(est_grid)
    stop("Cannot perform pairwise comparisons.")
  }

  if (length(specs$facs)) {
    est <- estimate_contrasts(est_grid[specs$facs])
    est <- eval(est, envir = rg@grid)
    names(est) <- paste(1:length(est))
  } else {
    est <- list("1" = rep(1 / nrow(rg@grid), nrow(rg@grid)))
  }
  
  freq <- !rg@misc$chains

  pmms <- list()

  if (freq) {
    pmms$pmmeans <- pmm_freq(rg, est, est_grid, estName = "pmmean",
      infer = c(TRUE, FALSE), famSize = nrow(est_grid))
  } else {
    pmms$pmmeans <- pmm_bayes(rg, est, est_grid, estName = "pmmean")
  }

  if (length(specs$by)) {
    by_num <- as.numeric(by_fac)
    est_grid <- est_grid[c(specs$by, setdiff(specs$vars, specs$by))]
    est.contr <- split(est, by_num)
    est_grid.contr <- split(est_grid, by_num)
    contr <- mlapply(est = est.contr, eg = est_grid.contr, fun = pairwise_contr)
    
    if (freq) {
      pmms[paste0("contrasts_", seq_along(contr))] <- mlapply(contr = contr,
        famSize = lengths(est.contr), same = list(rg = rg, estType = "pairs",
        adjust = "tukey", methDesc = "pairwise differences"), fun = pmm_freq)
    } else {
      pmms[paste0("contrasts_", seq_along(contr))] <- mlapply(contr = contr,
        same = list(rg = rg, estName = "pairs"), fun = pmm_bayes)
    }
    
  } else if (specs$pw) {
    contr <- pairwise_contr(est, est_grid)
    
    if (freq) {
      pmms$contrasts <- pmm_freq(rg, contr, estType = "pairs", adjust = "tukey",
        methDesc = "pairwise differences", famSize = nrow(est_grid))
    } else {
      pmms$contrasts <- pmm_bayes(rg, contr, estName = "pairs")
    }
  }

  if (any(lengths(rg@matlevs) > 1)) {
    warning("Some predictors are matrices with more than one column.",
      "\n  If these are from a call to poly(), results may be misleading.")
  }
  
  return(structure(pmms, class = "nauf.pmm.list",
    specs = list(variables = specs$vars, pairwise = specs$pw,
    averaged_over = setdiff(specs$avgfac, subset$cond),
    held_at_mean = specs$avgcov, conditioned_on = subset$cond,
    keep_NA = specs$keepNA, drop_NA = specs$dropNA, subset = subset$subset,
    by = specs$by, bayes = !freq, note = specs$note)))
}


pmm_freq <- function(rg, contr, eg = NULL, ...) {
  pmm <- lsmeans::contrast(rg, contr)
  
  if (!is.null(eg)) {
    pmm@grid <- as.data.frame(sapply(eg, paste))
    pmm@roles$predictors <- colnames(eg)
  }
  rownames(pmm@grid) <- NULL
  
  if (length(dots <- list(...))) pmm@misc[names(dots)] <- dots
  if (length(tran <- pmm@misc$orig.tran)) pmm@misc$tran <- tran
  
  return(pmm)
}


pmm_bayes <- function(rg, contr, eg = NULL, ...) {
  if (is.null(eg)) {
    nms <- names(contr)
    g <- data.frame(contrast = nms)
  } else {
    nms <- sapply(eg, paste)
    g <- as.data.frame(nms)
    for (j in 1:ncol(nms)) {
      nms[, j] <- paste(colnames(nms)[j], nms[, j], sep = ":")
    }
    nms <- apply(nms, 1, paste, collapse = ",")
  }
  rownames(g) <- NULL
  
  m <- do.call(rbind, contr) %*% rg@linfct
  
  samp <- make_stan_mcmc_array(samples = rg@post.beta %*% t(m),
    chains = rg@misc$chains, nms = nms)
    
  return(structure(list(names = g, contrasts = m, samples = samp,
    family = rg@misc$family, inv.lbl = rg@misc$inv.lbl, misc = list(...)),
    class = "nauf.pmm.stan"))
}


make_stan_mcmc_array <- function(samples, chains, nms = NULL) {
  if (!is.matrix(samples) || nrow(samples) %% chains) {
    stop("'samples' should be a matrix with number of rows divisible by 'chains'")
  }
  pars <- ncol(samples)
  if (is.null(nms)) {
    if (is.null(colnames(samples))) {
      nms <- paste0("V", 1:ncol(samples))
    } else {
      nms <- colnames(samples)
    }
  } else if (!is.vector(nms) || length(nms) != pars) {
    stop("'nms', if specified, must be a vector with an element for each",
      " column in 'samples'.")
  }
  iter <- nrow(samples) / chains
  
  a <- array(dim = c(iter, chains, pars), dimnames = list(iterations = NULL,
    chains = paste0("chain:", 1:chains), parameters = nms))
  for (ch in 1:chains) {
    a[, ch, ] <- samples[(iter * (ch - 1) + 1):(iter * ch), , drop = FALSE]
  }
  
  return(a)
}


#' @importFrom utils combn
pairwise_contr <- function(est, eg) {
  eg <- sapply(eg, as.character)
  
  combs <- utils::combn(length(est), 2)
  
  contr <- list_mat_cols(apply(combs, 2,
    function(x) est[[x[1]]] - est[[x[2]]]))
    
  names(contr) <- apply(combs, 2, function(x) paste0(
    paste(eg[x[1], ], collapse = ","),
    " - ",
    paste(eg[x[2], ], collapse = ",")
  ))
  
  return(contr)
}


check_specs <- function(specs, pw, xlev, mt, keepNA, by) {
  vv <- varnms(mt)

  if (inherits(specs, "formula")) {
    resp <- has_resp(specs <- stats::terms(specs))
    specs <- varnms(specs)
    if (is.null(specs)) stop("'specs' has no variables")

    if (pw <- (resp && specs[1] == "pairwise")) {
      specs <- specs[-1]
    } else if (resp) {
      stop("If 'specs' has a left hand side, it must be 'pairwise'")
    }
  }
  if (!is.character(specs)) {
    stop("'specs' should be a character vector or formula specifying variables")
  }
  if (!length(specs)) stop("'specs' has no variables")
  if (length(not_in <- setdiff(specs, vv))) {
    stop("The following variables are not valid:\n  ", add_quotes(not_in),
      "\nValid variables are:\n  ", add_quotes(vv))
  }

  all_facs <- names(xlev)
  # intersect since info$uf may contain facs in ranef but not fixef
  all_uf <- names(uflev <- xlev[intersect(names(xlev), names(attr(mt,
    "nauf.info")$uf))])
  all_nauf <- all_facs[sapply(xlev, anyNA)]
  
  if (!is.null(by)) {
    if (!is.character(by)) {
      stop("'by', if specified, must be a character vector")
    }
    if (length(setdiff(by, all_uf))) {
      stop("'by' should only contain names of unordered factors")
    }
    by <- unique(by)
    specs <- union(specs, by)
    if (!length(intersect(setdiff(specs, by), all_facs))) {
      stop("If 'by' is specified, 'specs' must contain at least one factor",
        " that is not in 'by'")
    }
    pw <- TRUE
  }
  
  # reorder for consistency; gets rid of possible duplicates
  specs <- intersect(vv, specs)

  facs <- intersect(specs, all_facs)
  avgfac <- setdiff(all_facs, facs)
  uf <- intersect(facs, all_uf)
  nauf <- intersect(uf, all_nauf)
  nums <- setdiff(specs, all_facs)
  avgcov <- setdiff(vv, c(all_facs, nums))

  note <- character()
  fmat <- attr(mt, "factors")
  specs_which <- which(colnames(fmat) == (specs_term <- paste(specs,
    collapse = ":")))
  if (!length(specs_which)) {
    note <- paste0("The interaction term '", specs_term,
      "' is not in the model.")

  } else {
    ord <- attr(mt, "order")
    specs_order <- ord[specs_which]
    has_specs <- sapply(fmat[specs, , drop = FALSE], all)
    ord[!has_specs] <- 0
    if (length(higher <- which(ord > specs_order))) {
      note <- paste(add_quotes(specs_term), "is included in higher order",
        "interaction(s)", add_quotes(colnames(fmat)[higher]))
    }
  }

  if (is.null(keepNA)) keepNA <- character()
  arg_keepNA <- keepNA
  keepNA <- intersect(keepNA, nauf)
  if (length(dropped <- setdiff(arg_keepNA, keepNA))) {
    warning("Dropped from 'keepNA' (not unordered factors with NA values ",
      "included in 'specs' or 'by':\n  ", add_quotes(dropped))
  }
  dropNA <- setdiff(nauf, keepNA)


  return(list(vars = specs, facs = facs, nums = nums, uf = uf, uflev = uflev,
    pw = pw, keepNA = keepNA, avgfac = avgfac, avgcov = avgcov, nauf = nauf,
    dropNA = dropNA, note = note, by = by))
}


grid_subset <- function(subset, specs) {
  if (length(specs$dropNA)) {
    est_keep <- drop_nas(specs$dropNA)
  } else {
    est_keep <- NULL
  }

  if (is.null(subset)) {
    return(list(keep = NULL, est_keep = est_keep, cond = character(),
      subset = NULL))
  }

  if (!is.list(subset) || !length(subset)) {
    stop("'subset' must be a list specifying groups")
  }
  if (!any(sapply(subset, is.list))) subset <- list(subset)

  msg <- mapply(check_named_list, subset, MoreArgs = list(nms = specs$uflev),
    SIMPLIFY = FALSE)
  msg <- msg[sapply(msg, function(u) !isTRUE(u))]
  if (length(msg)) {
    stop("'subset' not properly specified:\n  ", do.call(paste, c(msg,
      list(collapse = "\n  "))))
  }

  cond <- unname(unique(unlist(lapply(subset, names))))
  keep <- any_subset(subset)

  return(list(keep = keep, est_keep = est_keep, cond = cond, subset = subset))
}


estimate_contrasts <- function(eg) {
  est <- call("list")

  # otherwise level labels are replaced with integers
  eg <- as.data.frame(lapply(eg, as.character), stringsAsFactors = FALSE)

  for (i in 1:nrow(eg)) {
    est[[i + 1]] <- call("as_simplex")
    est[[i + 1]][[2]] <- list_to_subset(eg[i, , drop = FALSE])
  }

  return(est)
}


check_named_list <- function(x, nms) {
  if (!is.list(x)) return("Not a list")

  if (!length(x)) return("Empty list")

  if (is.null(nx <- names(x))) return("Not a named list")

  if (any(sapply(x, function(u) !is.character(u) && !all(is.na(u))))) {
    return("Not all list elements are character vectors")
  }

  if (is.list(nms)) {
    if (length(not_in <- setdiff(nx, names(nms)))) {
      return(paste(add_quotes(not_in),
        "is/are not unordered factors in the model"))
    }

    not_in <- mapply(setdiff, x, nms[nx], SIMPLIFY = FALSE)
    if (length(wrong <- which(lengths(not_in) > 0))) {
      nx <- nx[wrong[1]]
      wrong <- not_in[[wrong[1]]]
      return(paste(add_quotes(wrong), "is/are not valid levels for the factor",
        add_quotes(nx)))
    }

  } else {
    if (length(not_in <- setdiff(nx, nms))) {
      return(paste(add_quotes(not_in),
        "is/are not unordered factors in the model"))
    }
  }

  if (length(nx) != length(unique(nx))) {
    return("list names are not unique")
  }

  return(TRUE)
}


join_and <- function(conds) {
  if (!length(conds)) return(NULL)
  conds <- conds[sapply(conds, function(x) !is.null(x))]
  if (!length(conds)) return(NULL)

  if (length(conds) > 1) {
    cl <- call("&")
    cl[[3]] <- conds[[1]]
    cl[[2]] <- conds[[2]]
    if (length(conds) > 2) {
      for (k in 3:length(conds)) {
        cl <- substitute(a & b, list(a = cl, b = conds[[k]]))
      }
    }
  } else {
    cl <- conds[[1]]
  }

  return(cl)
}


join_or <- function(conds) {
  if (!length(conds)) return(NULL)
  conds <- conds[sapply(conds, function(x) !is.null(x))]
  if (!length(conds)) return(NULL)

  if (length(conds) > 1) {
    cl <- call("|")
    cl[[3]] <- conds[[1]]
    cl[[2]] <- conds[[2]]
    if (length(conds) > 2) {
      for (k in 3:length(conds)) {
        cl <- substitute(a | b, list(a = cl, b = conds[[k]]))
      }
    }
  } else {
    cl <- conds[[1]]
  }

  return(cl)
}


fac_in_levs <- function(fac, levs) {
  return(substitute(F %in% L, list(F = as.name(fac), L = levs)))
}


list_to_subset <- function(subspecs) {
  return(join_and(mapply(fac_in_levs, names(subspecs), subspecs,
    SIMPLIFY = FALSE)))
}


any_subset <- function(subspecs) {
  return(join_or(lapply(subspecs, list_to_subset)))
}


drop_nas <- function(facs) {
  return(substitute(!X, list(X = join_or(lapply(facs,
    function(fac) substitute(is.na(F), list(F = as.name(fac))))))))
}

