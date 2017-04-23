

#' Create a reference grid for a nauf model.
#'
#' Serves the same purpose as \code{\link[lsmeans]{ref.grid}}, but only takes
#' one argument (a \code{nauf} model).  See 'Details'.
#'
#' ADD DETAILS
#'
#' @param mod Any \code{nauf} model (see \code{\link{is.nauf.model}}).
#'
#' @importFrom pbkrtest Lb_ddf vcovAdj
#' @importFrom lsmeans ref.grid
#'
#' @export
nauf_ref.grid <- function(mod) {
  if (!is.nauf.model(mod)) stop("Must supply a nauf model")

  fenr <- stats::delete.response(terms(mod))
  info <- attr(fenr, "nauf.info")
  vars <- varnms(fenr)
  mf <- model.frame(mod)
  
  xlev <- lvs <- mlvs <- list()
  for (v in vars[vars %in% names(info$uf)]) {
    lvs[[v]] <- info$uf[[v]][[1]]
    if (info$hasna[v]) {
      lvs[[v]][length(lvs[[v]]) + 1] <- NA
      mf[[v]] <- addNA(mf[[v]])
    }
    xlev[[v]] <- lvs[[v]]
  }
  for (v in vars[vars %in% names(info$of)]) {
    xlev[[v]] <- lvs[[v]] <- info$of[[v]]$levels
  }
  for (v in vars[vars %in% names(info$num)]) {
    lvs[[v]] <- info$num[[v]]
  }
  for (v in vars[vars %in% names(info$mat)]) {
    mlvs[[v]] <- info$mat[[v]]
  }
  
  g <- expand.grid(lvs)
  g[names(mlvs)] <- lapply(mlvs, function(m) {
    matrix(m, nrow(g), length(m), byrow = TRUE)
  })
  
  mm <- model.matrix(fenr, g)
  asgn <- attr(mm, "assign")
  
  g[[".wgt."]] <- data.frame(xtabs(~ ., mf[names(xlev)]))$Freq
  
  summ <- summary(mod)
  
  rg <- list()
  class(rg) <- c("nauf.ref.grid", "list")
  temp_dat <- data.frame(y = 1:5, x = c(2, 5, 3, 6, 1))
  temp_mod <- lm(y ~ x, temp_dat)
  rg$ref.grid <- lsmeans::ref.grid(temp_mod, data = temp_dat)

  rg$ref.grid@model.info <- list(
    call = summ$call,
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
  rg$ref.grid@bhat <- summ$coefficients[, 1]
  rg$ref.grid@nbasis <- matrix()
  rg$ref.grid@V <- as.matrix(summ$vcov)
  if (is.nauf.lmerMod(mod)) {
    rg$ref.grid@dffun <- function(k, dfargs) {
      pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
    }
    rg$ref.grid@dfargs <- list(
      unadjV = vcov(mod),
      adjV = pbkrtest::vcovAdj(mod))
  } else if (is.nauf.glmerMod(mod)) {
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
    assign = asgn)
  rg$ref.grid@post.beta <- matrix()
  family <- get_family(mod)
  if (!isTRUE(all.equal(family, gaussian()))) {
    rg$ref.grid@misc$tran <- family$link
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


#' Predicted marginal means for nauf models.
#' 
#' ADD DESCRIPTION
#'
#' ADD DETAILS
#'
#' @param object A \code{\linkS4class{nauf.ref.grid}} or any \code{nauf} model
#'   (see \code{\link{is.nauf.model}}).
#'
#' @importFrom lsmeans contrast lsmeans pmmeans
#'
#' @export
nauf_pmmeans <- function(object, specs, keep_level = list(),
                         drop_level = list(), keep_group = list(),
                         drop_group = list(), ...) {
  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quote(names(dots)))
  
  if (!is.nauf.ref.grid(object)) object <- nauf_ref.grid(object)
  rg <- object$ref.grid
  
  if (!inherits(specs, "formula")) stop("'specs' must be a formula")
  b0 <- has_b0(specs <- stats::terms(specs))
  specs <- varnms(specs)
  if (is.null(specs)) stop("'specs' has no variables")
  if (pw <- (b0 && specs[1] == "pairwise")) {
    specs <- specs[-1]
  } else if (b0) {
    stop("If 'specs' has a left hand side, it must be 'pairwise'")
  }
  if (!length(specs)) stop("'specs' has no variables")
  if (length(not_in <- specs[!(specs %in% rg@roles$predictors)])) {
    stop("The following variables are not valid:\n  ", add_quotes(not_in),
      "Valid variables are:\n  ", add_quotes(rg@roles$predictors))
  }
  
  all_facs <- names(xlev <- rg@model.info$xlev)
  
  check_subset(xlev, keep_level, drop_level, keep_group, drop_group)
  keep <- grid_subset(rg@grid, keep_level, drop_level, keep_group,
    drop_group)
  rg@grid <- rg@grid[keep, , drop = FALSE]
  rg@linfct <- rg@linfct[keep, , drop = FALSE]

  nums <- specs[!(specs %in% all_facs)]
  facs <- specs[specs %in% all_facs]
  
  est_lvs <- xlev[facs]
  if (length(nums)) {
    g <- rg@grid[, -ncol(rg@grid), drop = FALSE]
    for (v in nums) {
      g[[v]] <- g[[v]] + 1
      est_lvs[[v]] <- "_inc_1"
    }
    rg@linfct <- model.matrix(rg@model.info$terms, g) - rg@linfct
  }
  est_grid <- expand.grid(est_lvs)
  if (length(facs)) {
    u <- unique(rg@grid[, facs, drop = FALSE])
    keep <- apply(apply(u, 1, function(i) match_row(est_grid, i)), 1, any)
    est_grid <- est_grid[keep, , drop = FALSE]
  }
  if (!nrow(est_grid)) {
    stop("No possible combinations satisfy subsetting conditions")
  }
  
  pmms <- list()
  est <- estimate_contrasts(rg@grid, est_grid, facs)
  pmms$pmmeans <- marginal_means(rg, est, est_grid)
  if (pw) {
    if (length(est) == 1) {
      pw <- paste("Specified 'pairwise', but there is only one combination",
        "that satisfies subsetting conditions.  Contrasts not computed")
      warning(pw)
    } else {
      pmms$contrasts <- pairwise_contrasts(rg, est, est_grid)
    }
  }
  
  if (length(rg@matlevs) && any(sapply(rg@matlevs, length) > 1)) {
    warning("Some matrix elements have more than one column.  If these are\n",
      "  from a call to poly(), results may be misleading.")
  }

  class(pmms) <- c("lsm.list", "list")
  attr(pmms, "nauf.specs") <- list(specs = specs, keep_level = keep_level,
    drop_level = drop_level, drop_group = drop_group, pairwise = pw)

  return(pmms)
}


marginal_means <- function(rg, est, est_grid) {
  est <- lsmeans::contrast(rg, est)
  est@roles$predictors <- colnames(est_grid)
  for (j in 1:ncol(est_grid)) est_grid[[j]] <- paste(est_grid[[j]])
  est@grid <- est_grid
  est@misc$estName <- "pmmean"
  est@misc$infer <- c(TRUE, FALSE)
  est@misc$famSize <- nrow(est_grid)
  if (length(tran <- est@misc$orig.tran)) est@misc$tran <- tran
  return(est)
}


#' @importFrom utils combn
pairwise_contrasts <- function(rg, est, est_grid) {
  est_grid <- sapply(est_grid, as.character)
  combs <- utils::combn(length(est), 2)
  contr <- list_mat_cols(apply(combs, 2, function(x) est[[x[1]]] - est[[x[2]]]))
  names(contr) <- apply(combs, 2, function(x) paste0(
    paste(est_grid[x[1], ], collapse = ","),
    " - ",
    paste(est_grid[x[2], ], collapse = ",")
  ))
  contr <- lsmeans::contrast(rg, contr)
  contr@misc$estType <- "pairs"
  contr@misc$adjust <- "tukey"
  contr@misc$methDesc <- "pairwise differences"
  contr@misc$famSize <- nrow(est_grid)
  if (length(tran <- contr@misc$orig.tran)) contr@misc$tran <- tran
  return(contr)
}


check_subset <- function(xlev, keep_level, drop_level, keep_group,
                              drop_group) {
  all_facs <- names(xlev)
  
  validate_level <- function(u) {
    if (!is.list(u)) return("Not a list")
    if (!length(u)) return(TRUE)
    nn <- names(u)
    not_in <- names(u)[!(names(u) %in% all_facs)]
    if (length(not_in)) {
      return(paste0("The following are not factors in the model:\n  ",
        add_quotes(not_in)))
    }
    u <- lapply(nn, function(k) u[[k]][!(u[[k]] %in% xlev[[k]])])
    u <- u[sapply(u, function(k) length(k) > 0)]
    if (!length(u)) return(TRUE)
    return(paste0("The following are not valid levels for '", names(u)[1],
      "':\n  ", add_quotes(u[[1]])))
  }
  
  validate_group <- function(g) {
    if (!is.list(g)) return("Not a list")
    if (!length(g)) return(TRUE)
    check <- sapply(g, is.list)
    if (any(!check)) return("Contains elements that are not lists")
    check <- sapply(g, length)
    if (any(!check)) return("Contains empty list elments")
    check <- lapply(g, validate_level)
    check <- check[sapply(check, function(k) !isTRUE(k))]
    if (!length(check)) return(TRUE)
    return(check[[1]])
  }
  
  if (!isTRUE(msg <- validate_level(keep_level))) {
    stop("Invalid 'keep_level'.  ", msg)
  }
  if (!isTRUE(msg <- validate_level(drop_level))) {
    stop("Invalid 'drop_level'.  ", msg)
  }
  if (!isTRUE(msg <- validate_group(keep_group))) {
    stop("Invalid 'keep_group'.  ", msg)
  }
  if (!isTRUE(msg <- validate_group(drop_group))) {
    stop("Invalid 'drop_group'.  ", msg)
  }
  
  invisible(NULL)
}


grid_subset <- function(grid, keep_level, drop_level, keep_group,
                            drop_group) {
  keep <- matrix(TRUE, nrow(grid), 4)
  
  if (length(keep_level)) {
    keep[, 1] <- match_row(grid, keep_level)
  }
  
  if (length(drop_level)) {
    keep[, 2] <- !match_row(grid, drop_level, f = any)
  }
  
  if (length(keep_group)) {
    keep[, 3] <- apply(sapply(keep_group, function(g) {
      match_row(grid, g)
    }), 1, any)
  }
  
  if (length(drop_group)) {
    keep[, 4] <- !apply(sapply(keep_group, function(g) {
      match_row(grid, g)
    }), 1, any)
  }
  
  return(apply(keep, 1, all))
}


estimate_contrasts <- function(grid, est_grid, facs) {
  if (!length(facs)) {
    return(list("1" = rep(1 / nrow(grid), nrow(grid))))
  }
  return(list_mat_cols(apply(apply(est_grid, 1,
    function(i) match_row(grid, i, facs)), 2, as_simplex)))
}

