
#' @export
nauf_grid <- function(mod) {
  if (!is.nauf(mod) || !inherits(mod, "lm")) {
    stop("must supply a nauf model")
  }

  mtnr <- delete.response(mod$terms)
  dc <- attr(mtnr, "dataClasses")
  cn <- attr(mod$terms, "factors")
  cn <- rownames(cn)[rowSums(cn) > 0]
  dc <- dc[names(dc) %in% cn]
  lvs <- xlvs <- mod$xlevels
  hasna <- attr(mod$terms, "hasna")
  dat <- mod$model
  for (n in names(hasna)[hasna]) {
    lvs[[n]] <- c(lvs[[n]], NA)
    dat[, n] <- addNA(dat[, n])
  }
  for (n in names(dc)[dc == "numeric"]) {
    lvs[[n]] <- mean(dat[, n])
  }
  g <- expand.grid(lvs)
  counts <- data.frame(xtabs(~ ., dat[, names(lvs)[sapply(lvs, length) > 1]]))  
  g[, ncol(g) + 1] <- counts$Freq
  colnames(g)[ncol(g)] <- ".wgt."
  
  mf <- nauf_on(g)
  attr(mf, "terms") <- mtnr
  mefc <- attr(mtnr, "mefc")
  for (f in names(mefc)) {
    contrasts(mf[, f]) <- mefc[[f]]$contrasts
  }
  mm <- model.matrix(mf)
  summ <- summary(mod)
  
  rgargs <- list(
    model.info = list(
      call = mod$call,
      terms = mod$terms,
      xlev = mod$xlevels),
    roles = list(
      predictors = all.vars(mtnr),
      responses = character(),
      multresp = character()),
    grid = g,
    levels = lvs,
    matlevs = list(),
    linfct = mm,
    bhat = summ$coefficients[, 1],
    nbasis = matrix(),
    V = summ$cov,
    dffun = function(k, dfargs) dfargs$df,
    dfargs = list(df = mod$df),
    misc = list(
      estName = "prediction",
      estType = "prediction",
      infer = c(FALSE, FALSE),
      level = 0.95,
      adjust = "none",
      famSize = ncol(g) - 1,
      avgd.over = character()),
    post.beta = matrix())

  d <- data.frame(y = 1:5, x = c(2,5,3,6,1))
  rg <- lm(y ~ x, d)
  rg <- lsmeans::ref.grid(rg, data = d)
  for (i in names(rgargs)) {
    slot(rg, i) <- rgargs[[i]]
  }
  
  return(rg)
}


#' @export
nauf_pmmeans <- function(object, facs, keep_level = NULL, drop_level = NULL,
                         keep_group = NULL, drop_group = NULL, pairwise = TRUE,
                         adjust = "tukey") {

  if (inherits(object, "ref.grid")) {
    if (is.nauf(object@model.info$terms)) {
      rg <- object
    } else {
      stop("must supply a ref.grid or model object fit with nauf functions")
    }
  } else if (is.nauf(object) && inherits(object, "lm")) {
    rg <- nauf_grid(object)
  } else {
    stop("must supply a ref.grid or model object fit with nauf functions")
  }
  if (!length(facs) || !is.character(facs) ||
  !all(facs %in% names(rg@levels))) {
    stop("'facs' must be a character vector specifying factor(s) in the model")
  }

  if (!is.null(keep_level)) {
    if (!is.list(keep_level) || !all(names(keep_level) %in% names(rg@levels))) {
      stop("'keep_level' must be a named list of factor levels to keep_level")
    }
    k <- rep(TRUE, nrow(rg@grid))
    for (n in names(keep_level)) {
      if(!all(keep_level[[n]] %in% rg@levels[[n]])) {
        stop("'keep_level' contains invalid factor levels")
      }
      k <- k & rg@grid[, n] %in% keep_level[[n]]
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }
  
  if (!is.null(drop_level)) {
    if (!is.list(drop_level) || !all(names(drop_level) %in% names(rg@levels))) {
      stop("'drop_level' must be a named list of factor levels to drop_level")
    }
    k <- rep(TRUE, nrow(rg@grid))
    for (n in names(drop_level)) {
      if(!all(drop_level[[n]] %in% rg@levels[[n]])) {
        stop("'drop_level' contains invalid factor levels")
      }
      k <- k & !(rg@grid[, n] %in% drop_level[[n]])
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }
  
  if (!is.null(keep_group)) {
    if (!is.list(keep_group) || is.data.frame(keep_group) ||
    !all(unlist(lapply(keep_group, function(x) is.list(x) &
    !is.data.frame(x))))) {
      stop("'keep_group' must be a list of groups expressed as ",
        "named lists of factor levels")
    }
    k <- rep(FALSE, nrow(rg@grid))
    for (g in keep_group) {
      k <- k | apply(sapply(names(g),
        function(x) rg@grid[, x] %in% g[[x]]), 1, all)
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  if (!is.null(drop_group)) {
    if (!is.list(drop_group) || is.data.frame(drop_group) ||
    !all(unlist(lapply(drop_group, function(x) is.list(x) &
    !is.data.frame(x))))) {
      stop("'drop_group' must be a list of groups expressed as ",
        "named lists of factor levels")
    }
    k <- rep(FALSE, nrow(rg@grid))
    for (g in drop_group) {
      k <- k | apply(sapply(names(g),
        function(x) !(rg@grid[, x] %in% g[[x]])), 1, any)
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  estgrid <- expand.grid(rg@levels[facs])
  k <- rep(TRUE, nrow(estgrid))
  est <- list()
  for (i in 1:nrow(estgrid)) {
    est[[i]] <- rep(TRUE, nrow(rg@grid))
    for (j in facs) {
      est[[i]] <- est[[i]] & (rg@grid[, j] %in% estgrid[i, j])
    }
    if ((n <- sum(est[[i]]))) {
      est[[i]] <- est[[i]] / n
    } else {
      k[i] <- FALSE
    }
  }
  if (sum(k) < 1) {
    stop("invalid specification; must be at least one group")
  }
  est <- est[k]
  estgrid <- estgrid[k, , drop = FALSE]
  names(est) <- apply(estgrid, 1, function(x) paste(x, collapse = ","))
  est.lsm <- lsmeans::contrast(rg, est)
  est.lsm@roles$predictors <- "estimate"
  colnames(est.lsm@grid) <- "estimate"
  est.lsm@misc$estName <- "pmmean"
  est.lsm@misc$infer <- c(TRUE, FALSE)
  est.lsm@misc$famSize <- nrow(est.lsm@grid)
  
  res <- list(pmmeans = est.lsm)
  class(res) <- c("lsm.list", "list")
  
  if (pairwise & sum(k) > 1) {
    combs <- combn(length(est), 2)
    contr <- as.list(as.data.frame(apply(combs, 2,
      function(x) est[[x[1]]] - est[[x[2]]])))
    names(contr) <- apply(combs, 2,
      function(x) paste(names(est)[x], collapse = " - "))
    contr <- lsmeans::contrast(rg, contr)
    contr@misc$estType <- "pairs"
    contr@misc$adjust <- adjust
    contr@misc$methDesc <- "pairwise differences"
    contr@misc$famSize <- est.lsm@misc$famSize
    res$contrasts <- contr
  }
  
  attr(res, "specs") <- list(
    facs = facs,
    keep_level = keep_level,
    keep_group = keep_group,
    drop_level = drop_level,
    drop_group = drop_group)

  return(nauf_on(res))
}

