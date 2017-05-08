

#' Predicted marginal means for \code{nauf} models.
#'
#' Create a reference grid for a \code{nauf} model with \code{nauf_ref.grid},
#' and use the resultnig \code{nauf.ref.grid} as the \code{object} argument to
#' \code{nauf_pmmeans} to obtain predicted marginal means and pairwise
#' comparisons, optionally conditioning these predictions on certain subsets of
#' the data via the \code{subset} argument.
#'
#' A reference grid creates a data frame which contains all possible
#' combinations of the factors in a regression model, holding all covariates
#' at their mean values.  There are many options for
#' \code{\link[lsmeans]{ref.grid}} which are not currently supported for
#' \code{nauf} models.  The main functionality which is not currently supported
#' is that the reference grid cannot be created specifying certain levels for
#' variables (i.e. the \code{at} argument; this is handled through the
#' \code{subset} argument to \code{nauf_pmmeans}).  A direct call to
#' \code{\link[lsmeans]{ref.grid}} will result in warnings (or possibly errors),
#' and inference made with the resulting object will be misleading and/or
#' incorrect.  Only \code{nauf_ref.grid} should be used. The \code{nauf.ref.grid}
#' returned by \code{nauf_ref.grid} can then be used as the \code{object}
#' argument to \code{nauf_pmmeans} to obtain predicted marginal means and
#' pairwise comparisons with p-values that adjust for familywise error rate.
#'
#' The \code{specs} and \code{pairwise} arguments to \code{nauf_pmmeans}
#' indicate what variables marginal means should be calculated for and wheter
#' pairwise comparisons of these means should be made.  If \code{specs} is a
#' character vector, then \code{pairwise} is used; if \code{specs} is a formula,
#' then the full iteraction of the terms on the right hand side of the formula
#' is considered, and the left hand side is used to indicate pairwise
#' comparisons.  For example (where \code{rg} is a \code{nauf.ref.grid}):
#'
#' \preformatted{
#' # all of these calculate pmm's for each combination of the factors f1 and f2
#' # but not pairwise comparisons
#' nauf_pmmeans(rg, c("f1", "f2"))
#' nauf_pmmeans(rg, ~ f1 + f2)
#' nauf_pmmeans(rg, ~ f1 * f2)
#' nauf_pmmeans(rg, ~ f1:f2)
#'
#' # all of these calculate the same pmm's, and additionally pairwise comparions
#' nauf_pmmeans(rg, c("f1", "f2"), pairwise = TRUE)
#' nauf_pmmeans(rg, pairwise ~ f1 + f2)
#' nauf_pmmeans(rg, pairwise ~ f1 * f2)
#' nauf_pmmeans(rg, pairwise ~ f1:f2)
#' }
#'
#' If \code{specs} indicates a single covariate, the effect of an increase of
#' \code{1} in the covariate is computed.  If \code{specs} indicates multiple
#' covariates, the effect of a simultaneous increase of \code{1} in all of the
#' covariates is computed.  If \code{specs} indicates a combination of factors
#' and covariate(s), the the effect of an increase of \code{1} for the
#' covariates is calcualted for each level of the full interaction of the
#' factors.
#'
#' The reference grid returned by \code{nauf_ref.grid} contains combinations of
#' factors which are not actually possible in the data set.  For example,
#' if factor \code{f1} has levels \code{A} and \code{B}, and factor \code{f2}
#' is \code{NA} when \code{f1 = A}, and takes values \code{C} and \code{D}
#' when \code{f1 = B}, the reference grid will still contain the combinations
#' \code{f1 = A, f2 = C}; \code{f1 = A, f2 = D}; and \code{f1 = B, f2 = NA},
#' even though these combinations are not possible.  This is because it is
#' impossible to know without the user's knowledge which combinations
#' make sense.  In many cases, this is inconsequential for the computation
#' of predicted marginal means, since the coding of unordered factors in
#' \code{nauf} regressions will average over the effects.  In cases where
#' these rows in the reference grid will cause invalid estimates and pairwise
#' comparisons, the \code{subset} argument can be used in the call to
#' \code{nauf_pmmeans} to ensure only the correct subsets are considered.
#' The default for the \code{subset} argument is \code{NULL}, indicating that
#' the the entire reference grid should be considered.  If not \code{NULL}, then
#' \code{subset} must be a list which defines the valid subsets as lists of
#' named character vectors, where the name of the character vector is an
#' unordered factor in the model, and the vector itself contains the levels
#' which define the subset (including \code{NA} in the case of factors which
#' have \code{NA} values; when \code{NA} is specified as a level, there should
#' be no quotes around it).  Any row in the reference grid which matches the
#' definition of at least one of the groups defined in \code{subset} is kept,
#' and all others are dropped.  So, continuing with the \code{f1} and \code{f2}
#' example, if \code{f2 = NA} corresponds to \code{f2 = D} in meaning, and is
#' coded as \code{NA} because all \code{f1 = A} observations are by necessity
#' \code{f2 = D}, then to analyze the effect of \code{f1}, we want to compare
#' the groups \code{f1 = A, f2 = NA} and \code{f1 = B, f2 = D}, which we could
#' do with the following call:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = "D")))
#' }
#'
#' This would produce an estimate for \code{f1 = A} and \code{f1 = B}, but
#' conditioning on the subset where \code{f1} is truly contrastive based on
#' \code{f2}.  If, on the other hand, \code{f2 = NA} does not correspond in
#' interpretation to either \code{f2 = C} or \code{f2 = D}, but rather indicates
#' that \code{f2} is simply not meaningful when \code{f1 = A}, we would want to
#' average over the effect of \code{f2} within \code{f1 = B}, and compare this
#' result to \code{f1 = A, f2 = NA}, which we could do with the following call:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))))
#' }
#'
#' In this case, the second sub-list in the \code{subset} list indicates that if
#' \code{f1 = B} and either \code{f2 = C} or \code{f2 = D}, then it belongs to
#' the second subset.  In this case, the \code{subset} argument is actually not
#' necessary, since for \code{f1 = A}, we want to \emph{not consider} the effect
#' \code{f2}, and for \code{f2 = B}, we want to \emph{average over all possible
#' levels} of \code{f2}, and these are actually the same thing computationally
#' for unordered factors in \code{nauf} models.  That is, we would get the same
#' result with:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f1")
#' }
#'
#' Generally speaking, if all of the factors in \code{specs} do \emph{not}
#' contain \code{NA} values, then the \code{subset} argument is unnecessary.
#' If any of the factors in \code{specs} \emph{do} contain \code{NA} values,
#' then you will almost always want to use the \code{subset} argument.  Now
#' consider that we are interested now in \code{f2}.  Because \code{f2} is only
#' contrastive when \code{f1 = B}, we probably want to call:
#'
#' \preformatted{
#' # note that because there are not multiple subsets being specified, you
#' # don't have to specify subset = list(list(f1 = "B")); nauf_pmmeans will
#' # assume list(f1 = "B") means list(list(f1 = "B"))
#' nauf_pmmeans(rg, "f2", subset = list(f1 = "B"))
#' }
#'
#' This call will produce two estimates, one for \code{f2 = C} and one for
#' \code{f2 = D}, conditioning on \code{f1 = B}.  There will be no estimate for
#' \code{f2 = NA} because, by default, no estiamtes are produced for
#' combinations of factors where one factor is \code{NA}.  If we wanted to
#' compare the three possible groups (i.e. \code{f1 = A, f2 = NA};
#' \code{f1 = B, f2 = C}; and \code{f1 = B, f2 = D}), then we could additionally
#' use the \code{na_as_level} argument and change our \code{subset}:
#'
#' \preformatted{
#' nauf_pmmeans(rg, "f2", subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))),
#'   na_as_level = "f2")
#'
#' # this gives the same estimates, but the output will also show the
#' # corresponding level of f1, which is more transparent
#' nauf_pmmeans(rg, c("f1", "f2"), subset = list(
#'   list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D"))),
#'   na_as_level = "f2")
#' }
#'
#' The easiest way to use the \code{subset} argument is to create a list that
#' defines valid subsets for different regression terms of interest outside of
#' \code{nauf_pmmeans}, and then using the relevant element of the list in the
#' \code{nauf_pmmeans} call.  For example:
#'
#' \preformatted{
#' pmmsubs <- list()
#' pmmsubs$f1 <- list(list(f1 = "A", f2 = NA), list(f1 = "B", f2 = c("C", "D")))
#' pmmsubs$f2 <- list(f1 = "B")
#'
#' nauf_pmmeans(rg, "f1", subset = pmmsubs$f1)
#'
#' nauf_pmmeans(rg, "f2", subset = pmmsubs$f2)
#'
#' nauf_pmmeans(rg, c("f1", "f2"), subset = pmmsubs$f1, na_as_level = "f2")
#' }
#'
#' This way you can just define the different subsets of the data once and not
#' have to think about it at every \code{nauf_pmmeans} call.
#'
#' @param mod A regression model fit with \code{nauf} contrasts.
#' @param KR Only applies when \code{mod} is a \code{\linkS4class{nauf.lmerMod}}
#'   fit with \code{REML = TRUE}. If \code{KR = TRUE}, then the Kenward-Roger
#'   approximation is used to calculate degrees of freedom.  If
#'   \code{KR = FALSE} (the default), then the Satterthwaite approximation is
#'   used.  When \code{mod} is a
#'   \code{\linkS4class{nauf.lmerMod}} fit with \code{REML = FALSE}, then the
#'   Satterthwaite approximation is always used.  The Kenward-Roger method
#'   is implemented with \code{\link[pbkrtest]{Lb_ddf}} and the Satterthwaite
#'   method is implemented with \code{\link[lmerTest]{calcSatterth}}.
#' @param object A \code{nauf.ref.grid} object created with
#'   \code{nauf_ref.grid}.
#' @param specs The fixed effects for which the full interaction term should be
#'   considered in the calculation of predicted marginal means. The preferred
#'   method is to specify the variables as a character vector.  However, they
#'   can also be specified on the right hand side of a formula, optionally with
#'   the keyword \code{pairwise} on the left hand side to indicate that pairwise
#'   comparions should be performed.
#' @param pairwise A logical (default \code{FALSE}) indicating whether pairwise
#'   comparisons of the predicted marginal means should be performed.  If
#'   \code{specs} is a formula, then the \code{pairwise} argument is ignored and
#'   the left hand side of the formula is used to determined whether pairwise
#'   comparsions should be made.
#' @param subset A list indicating which subsets of the reference grid should be
#'   considered in the calculation of the predicted marginal means. See
#'   'Details'.
#' @param na_as_level A character vector of unordered factors in \code{specs}
#'   that have \code{NA} values that should be considered as levels. The default
#'   \code{NULL} indicates that \code{NA} should not be considered as a level
#'   for any unordered factors in \code{specs}. See 'Details'.
#' @param ... Additional arguments are ignored with a warning.
#'
#' @return
#' \code{nauf_ref.grid} returns a \code{nauf.ref.grid} object, which is just a
#' list with one element \code{ref.grid} of class
#' \code{\link[lsmeans]{ref.grid-class}}.  This reference grid should not be
#' used directly with \code{\link[lsmeans]{lsmeans}}, but rather only with
#' \code{nauf_pmmeans}.
#'
#' \code{nauf_pmmeans} returns a \code{nauf.pmm} object, which is a list
#' inheriting from \code{lsm.list} with an additional attribute
#' \code{nauf.specs} containing information about the variables and subsets from
#' the call to the function.  The \code{nauf.pmm} list contains an element
#' \code{pmmeans}, and, if pairwise comparisons were made, a second element
#' \code{contrasts}, both of which are \code{\link[lsmeans]{lsmobj-class}}
#' objects.  The \code{nauf.pmm} object has \code{summary} and \code{print}
#' methods which print information from the \code{nauf.specs} attribute, and
#' then call the \code{\link[lsmeans]{summary.ref.grid}} methods (to which
#' arguments \code{infer}, \code{type}, \code{adjust}, etc. can be passed).
#'
#' @seealso \code{\link{nauf_contrasts}}, \code{\link{nauf_glm}}, and
#'   \code{\link{nauf_glmer}}.
#'
#' @name nauf-pmmeans
NULL


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

  # reorder for consistency
  lvs <- lvs[vars[vars %in% names(lvs)]]
  xlev <- xlev[vars[vars %in% names(xlev)]]

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
  if (is.nauf.merMod(mod)) {
    rg$ref.grid@V <- as.matrix(summ$vcov)
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
  if (!is.linear(family <- get_family(mod))) {
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


#' @rdname nauf-pmmeans
#'
#' @importFrom lsmeans contrast lsmeans pmmeans
#' @importFrom utils combn
#'
#' @export
nauf_pmmeans <- function(object, specs, pairwise = FALSE, subset = NULL,
                         na_as_level = NULL, ...) {
  if (!is.nauf.ref.grid(object)) stop("must supply a nauf.ref.grid")
  rg <- object$ref.grid

  dots <- list(...)
  if (length(dots)) warning("Ignoring arguments: ", add_quotes(names(dots)))

  specs <- check_specs(specs, pairwise, object$ref.grid@model.info$xlev,
    object$ref.grid@model.info$terms, na_as_level)
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

  if (length(specs$facs)) {
    est <- estimate_contrasts(est_grid[, specs$facs, drop = FALSE])
    est <- eval(est, envir = object$ref.grid@grid)
    names(est) <- paste(1:length(est))
  } else {
    est <- list("1" = rep(1 / nrow(object$ref.grid@grid),
      nrow(object$ref.grid@grid)))
  }

  pmms <- list()

  pmms$pmmeans <- lsmeans::contrast(object$ref.grid, est)
  pmms$pmmeans@roles$predictors <- colnames(est_grid)
  pmms$pmmeans@grid <- est_grid
  for (j in 1:ncol(est_grid)) {
    pmms$pmmeans@grid[[j]] <- paste(est_grid[[j]])
  }
  pmms$pmmeans@misc$estName <- "pmmean"
  pmms$pmmeans@misc$infer <- c(TRUE, FALSE)
  pmms$pmmeans@misc$famSize <- nrow(est_grid)
  if (length(tran <- pmms$pmmeans@misc$orig.tran)) {
    pmms$pmmeans@misc$tran <- tran
  }

  if (specs$pw) {
    if (length(est) == 1) {
      specs$pw <- paste("Specified 'pairwise', but there is only one",
        "combination that satisfies subsetting conditions.")
      warning(specs$pw)

    } else {
      est_grid <- sapply(est_grid, as.character)
      combs <- utils::combn(length(est), 2)
      est <- list_mat_cols(apply(combs, 2,
        function(x) est[[x[1]]] - est[[x[2]]]))
      names(est) <- apply(combs, 2, function(x) paste0(
        paste(est_grid[x[1], ], collapse = ","),
        " - ",
        paste(est_grid[x[2], ], collapse = ",")
      ))

      pmms$contrasts <- lsmeans::contrast(object$ref.grid, est)
      pmms$contrasts@misc$estType <- "pairs"
      pmms$contrasts@misc$adjust <- "tukey"
      pmms$contrasts@misc$methDesc <- "pairwise differences"
      pmms$contrasts@misc$famSize <- nrow(est_grid)
      if (length(tran <- pmms$contrasts@misc$orig.tran)) {
        pmms$contrasts@misc$tran <- tran
      }
    }
  }

  if (length(object$ref.grid@matlevs) && any(sapply(object$ref.grid@matlevs,
  length) > 1)) {
    warning("Some matrix elements have more than one column.  If these are\n",
      "  from a call to poly(), results may be misleading.")
  }

  class(pmms) <- c("nauf.pmm", "lsm.list", "list")
  attr(pmms, "nauf.specs") <- list(
    variables = specs$vars,
    pairwise = specs$pw,
    averaged_over = setdiff(specs$avgfac, subset$cond),
    held_at_mean = specs$avgcov,
    conditioned_on = subset$cond,
    keep_NA = specs$keepNA,
    drop_NA = specs$dropNA,
    subset = subset$subset,
    note = specs$note)

  return(pmms)
}


check_specs <- function(specs, pw, xlev, mt, keepNA) {
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

  # reorder for consistency; gets rid of possible duplicates
  specs <- vv[vv %in% specs]

  all_facs <- names(xlev)
  # intersect since info$uf may contain facs in ranef but not fixef
  all_uf <- names(uflev <- xlev[intersect(names(xlev), names(attr(mt,
    "nauf.info")$uf))])
  all_nauf <- all_facs[sapply(xlev, anyNA)]

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
      "included in 'specs':\n  ", add_quotes(dropped))
  }
  dropNA <- setdiff(nauf, keepNA)


  return(list(vars = specs, facs = facs, nums = nums, uf = uf, uflev = uflev,
    pw = pw, keepNA = keepNA, avgfac = avgfac, avgcov = avgcov, nauf = nauf,
    dropNA = dropNA, note = note))
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

