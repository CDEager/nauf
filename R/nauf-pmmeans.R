

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

  xlev <- lvs <- mlvs <- list()
  for (v in intersect(vars, names(info$uf))) {
    lvs[[v]] <- info$uf[[v]][[1]]
    if (info$hasna[v]) {
      lvs[[v]][length(lvs[[v]]) + 1] <- NA
      mf[[v]] <- addNA(mf[[v]])
    }
    xlev[[v]] <- lvs[[v]]
  }
  for (v in intersect(vars, names(info$of))) {
    xlev[[v]] <- lvs[[v]] <- info$of[[v]]$levels
  }
  for (v in intersect(vars, names(info$num))) {
    lvs[[v]] <- info$num[[v]]
  }
  for (v in intersect(vars, names(info$mat))) {
    mlvs[[v]] <- info$mat[[v]]
  }

  # reorder for consistency
  lvs <- lvs[intersect(vars, names(lvs))]
  xlev <- xlev[intersect(vars, names(xlev))]

  g <- expand.grid(lvs)
  g[names(mlvs)] <- lapply(mlvs, function(m) {
    matrix(m, nrow(g), length(m), byrow = TRUE)
  })

  mm <- model.matrix(fenr, g)
  asgn <- attr(mm, "assign")

  g[[".wgt."]] <- data.frame(xtabs(~ ., mf[names(xlev)]))$Freq

  if (!(bayes <- is.nauf.stanreg(mod))) summ <- summary(mod)

  rg <- list()
  class(rg) <- "nauf.ref.grid"
  temp_dat <- data.frame(y = 1:5, x = c(2, 5, 3, 6, 1))
  temp_mod <- lm(y ~ x, temp_dat)
  rg$ref.grid <- lsmeans::ref.grid(temp_mod, data = temp_dat)

  rg$ref.grid@model.info <- list(
    call = if (bayes) mod$call else summ$call,
    terms = fenr,
    xlev = xlev)
  rg$ref.grid@roles <- list(
    predictors = vars,
    responses = character(),
    multresp = character())
  rg$ref.grid@grid <- g
  rg$ref.grid@levels <- lvs
  rg$ref.grid@matlevs <- mlvs
  rg$ref.grid@linfct <- mm
  rg$ref.grid@bhat <- if (bayes) numeric() else summ$coefficients[, 1]
  rg$ref.grid@nbasis <- matrix()
  if (is.nauf.merMod(mod)) {
    rg$ref.grid@V <- as.matrix(summ$vcov)
  } else if (bayes) {
    rg$ref.grid@V <- matrix()
  } else {
    rg$ref.grid@V <- vcov(mod)
  }
  if (is.nauf.lmerMod(mod)) {
    if (!lme4::isREML(mod)) KR <- FALSE
    if (KR) {
      rg$ref.grid@dffun <- function(k, dfargs) {
        pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
      }
      rg$ref.grid@dfargs <- list(
        unadjV = vcov(mod),
        adjV = pbkrtest::vcovAdj(mod))
    } else {
      rg$ref.grid@dffun <- function(k, dfargs) {
        lmerTest::calcSatterth(dfargs$object, k)$denom
      }
      rg$ref.grid@dfargs <- list(object = mod)
    }
  } else if (is.nauf.glmerMod(mod) || bayes) {
    rg$ref.grid@dffun <- function(k, dfargs) NA
    rg$ref.grid@dfargs <- list()
  } else {
    rg$ref.grid@dffun <- function(k, dfargs) dfargs$df
    rg$ref.grid@dfargs <- list(df = mod$df)
  }
  rg$ref.grid@misc <- list(
    estName = "prediction",
    estType = "prediction",
    infer = c(FALSE, FALSE),
    level = 0.95,
    adjust = "none",
    famSize = ncol(g) - 1,
    avgd.over = character(),
    assign = asgn,
    post.beta = if (bayes) fixef(mod, TRUE, FALSE) else NULL)
  rg$ref.grid@post.beta <- matrix()
  if (!is.linear(family <- get_family(mod))) {
    rg$ref.grid@misc$tran <- family$link
    rg$ref.grid@misc$family <- family
    family <- family$family
    
    rate <- length(grep("poisson", family)) + length(grep("Negative", family))
    prob <- length(grep("binomial", family))
    if (rate) {
      rg$ref.grid@misc$inv.lbl <- "rate"
    } else if (prob) {
      rg$ref.grid@misc$inv.lbl <- "prob"
    } else {
      rg$ref.grid@misc$inv.lbl <- "response"
    }
  }

  return(rg)
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

  specs <- check_specs(specs, pairwise, object$ref.grid@model.info$xlev,
    object$ref.grid@model.info$terms, na_as_level, by)
  subset <- grid_subset(subset, specs)

  if (!is.null(subset$keep)) {
    keep <- eval(subset$keep, envir = object$ref.grid@grid)
    object$ref.grid@grid <- object$ref.grid@grid[keep, , drop = FALSE]
    object$ref.grid@linfct <- object$ref.grid@linfct[keep, , drop = FALSE]
  }

  if (length(specs$facs)) {
    est_grid <- unique(object$ref.grid@grid[, specs$facs, drop = FALSE])
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
      object$ref.grid@grid[[v]] <- object$ref.grid@grid[[v]] + 1
    }
    object$ref.grid@linfct <- model.matrix(object$ref.grid@model.info$terms,
      object$ref.grid@grid) - object$ref.grid@linfct
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
    est <- estimate_contrasts(est_grid[, specs$facs, drop = FALSE])
    est <- eval(est, envir = object$ref.grid@grid)
    names(est) <- paste(1:length(est))
  } else {
    est <- list("1" = rep(1 / nrow(object$ref.grid@grid),
      nrow(object$ref.grid@grid)))
  }
  
  freq <- is.null(object$ref.grid@misc$post.beta)

  pmms <- list()

  if (freq) {
    pmms$pmmeans <- pmm_freq(object, est, est_grid, list(estName = "pmmean",
      infer = c(TRUE, FALSE), famSize = nrow(est_grid)))
  } else {
    pmms$pmmeans <- pmm_bayes(object, est, est_grid, list(estName = "pmmean"))
  }

  if (length(specs$by)) {
    by_num <- as.numeric(by_fac)
    est_grid <- est_grid[c(specs$by, setdiff(specs$vars, specs$by))]
    contr <- lapply(1:max(by_num), function(x) {
      w <- which(by_num == x)
      return(pairwise_contr(est[w], est_grid[w, , drop = FALSE]))
    })
    
    if (freq) {
      for (i in 1:length(contr)) {
        pmms[[paste("contrasts", i, sep = "_")]] <- pmm_freq(object,
          contr[[i]], NULL, list(estType = "pairs", adjust = "tukey",
          methDesc = "pairwise differences", famSize = sum(by_num == i)))
      }
    } else {
      for (i in 1:length(contr)) {
        pmms[[paste("contrasts", i, sep = "_")]] <- pmm_bayes(object,
          contr[[i]], NULL, list(estName = "pairs"))
      }
    }
    
  } else if (specs$pw) {
    contr <- pairwise_contr(est, est_grid)
    
    if (freq) {
      pmms$contrasts <- pmm_freq(object, contr, NULL, list(estType = "pairs",
        adjust = "tukey", methDesc = "pairwise differences",
        famSize = nrow(est_grid)))
    } else {
      pmms$contrasts <- pmm_bayes(object, contr, NULL, list(estName = "pairs"))
    }
  }

  if (length(object$ref.grid@matlevs) && any(sapply(object$ref.grid@matlevs,
  length) > 1)) {
    warning("Some matrix elements have more than one column.  If these are\n",
      "  from a call to poly(), results may be misleading.")
  }

  class(pmms) <- c("nauf.pmm", if (freq) "lsm.list", "list")
  attr(pmms, "nauf.specs") <- list(
    variables = specs$vars,
    pairwise = specs$pw,
    averaged_over = setdiff(specs$avgfac, subset$cond),
    held_at_mean = specs$avgcov,
    conditioned_on = subset$cond,
    keep_NA = specs$keepNA,
    drop_NA = specs$dropNA,
    subset = subset$subset,
    by = specs$by,
    bayes = !freq,
    note = specs$note)

  return(pmms)
}


pmm_freq <- function(rg, contr, eg = NULL, misc = NULL) {
  pmm <- lsmeans::contrast(rg$ref.grid, contr)
  
  if (!is.null(eg)) {
    pmm@grid <- as.data.frame(sapply(eg, paste))
    pmm@roles$predictors <- colnames(eg)
  }
  rownames(pmm@grid) <- NULL
  
  pmm@misc[names(misc)] <- misc
  if (length(tran <- pmm@misc$orig.tran)) pmm@misc$tran <- tran
  
  return(pmm)
}


pmm_bayes <- function(rg, contr, eg = NULL, misc = NULL) {
  pmm <- list()
  
  pmm$model.info <- rg$ref.grid@model.info
  
  if (is.null(eg)) {
    pmm$grid <- data.frame(contrast = names(contr))
  } else {
    pmm$grid <- as.data.frame(sapply(eg, paste))
  }
  rownames(pmm$grid) <- NULL
  
  pmm$linfct <- do.call(rbind, contr) %*% rg$ref.grid@linfct
  pmm$mcmc <- pmm$linfct %*a% rg$ref.grid@misc$post.beta
  
  pmm$misc <- list(tran = rg$ref.grid@misc$tran,
    inv.lbl = rg$ref.grid@misc$inv.lbl)
  pmm$misc[names(misc)] <- misc
  pmm$family <- rg$ref.grid@misc$family
    
  class(pmm) <- "nauf.postmm"
  return(pmm)
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
    stop("'subset' not properly specified:\n  ", msg[[1]])
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

